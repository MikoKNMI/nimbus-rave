'''
Copyright (C) 2014- Swedish Meteorological and Hydrological Institute (SMHI)

This file is part of RAVE.

RAVE is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RAVE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with RAVE.  If not, see <http://www.gnu.org/licenses/>.

'''
## Python interface to be able to perform tiled compositing. This functionallity relies on the use of multiprocessing and
# storing intermediate files representing tiles. When the other processes are finishied with the processing, the result is
# put back into the full area again.
#

## @file
## @author Anders Henja, SMHI (original odc version Daniel Michelson, SMHI) 
## @date 2014-09-09
## @author Michal Koutek, KNMI (Fixes needed for OPERA NIMBUS)
## @date 2022-01-06

import rave_tile_registry
import compositing
import multiprocessing, rave_mppool
import math, os
import rave_pgf_logger
import area_registry
import numpy
import _rave, _area, _pycomposite, _projection, _raveio, _polarscan, _polarvolume, _transform, _cartesianvolume
from rave_defines import CENTER_ID, GAIN, OFFSET
from rave_defines import RAVE_TILE_COMPOSITING_PROCESSES, RAVE_QUALITY_CONTROL_PROCESSES
import rave_tempfile
import odim_source
import rave_pgf_quality_registry

logger = rave_pgf_logger.create_logger()

##
# The area registry to be used by this composite generator.
my_area_registry = area_registry.area_registry()

##
# The basic area definition that should be transfered to the tiling compositing instance.
# This definition will be pickled and sent to the receiving product generator.
#
class tiled_area_definition(object):
  def __init__(self, id, pcsdef, xscale, yscale, xsize, ysize, extent):
    self.id = id
    self.pcsdef = pcsdef
    self.xscale = xscale
    self.yscale = yscale
    self.xsize = xsize
    self.ysize = ysize
    self.extent = extent
  
  def __repr__(self):
    return "<%s, scale=%f * %f, size=%d * %d, extent=%s />"%(self.pcsdef, self.xscale, self.yscale, self.xsize, self.ysize, str(self.extent))


##
# Stores the objects as uncompressed temporary files on disc
# @param objects: a disctionary with filenames as keys and polar objects as values
# @return a dictionary with temporary filenames as keys and polar objects as values
# @throws Exception on error, for example if we run out of disc space
def store_temporary_files(objects):
  tempobjects={}
  try:
    for k in objects.keys():
      rio = _raveio.new()
      rio.object = objects[k]
      rio.compression_level = 0
      rio.fcp_istorek = 1
      rio.fcp_metablocksize = 0
      rio.fcp_sizes = (4,4)
      rio.fcp_symk = (1,1)
      rio.fcp_userblock = 0
      fileno, rio.filename = rave_tempfile.mktemp(suffix='.h5', close="True")
      tempobjects[rio.filename] = rio.object
      rio.save()
  except Exception as e:
    for tmpfile in tempobjects.keys():
      try:
        os.unlink(tmpfile)
      except:
        pass
    raise e
  return tempobjects

##
# The argument wrapper so that the arguments can be transfered to the composite generator taking care of the tile.
#
class multi_composite_arguments(object):
  ##
  # Constructor
  def __init__(self):
    self.xscale = 2000.0
    self.yscale = 2000.0
    self.detectors = []
    self.filenames = []
    self.ignore_malfunc = False
    self.prodpar = None
    self.product = _rave.Rave_ProductType_PCAPPI
    self.height = 1000.0
    self.elangle = 0.0
    self.range = 200000.0
    self.selection_method = _pycomposite.SelectionMethod_NEAREST
    self.interpolation_method = _pycomposite.InterpolationMethod_NEAREST
    self.qitotal_field = None
    self.applygra = False
    self.zr_A = 200.0
    self.zr_b = 1.6    
    self.applygapfilling = False
    self.applyctfilter = False
    self.quantity = "DBZH"
    self.gain = GAIN
    self.offset = OFFSET
    self.minvalue = -30.0
    self.reprocess_quality_field = False 
    self.area_definition = None
    self.verbose = False
    self.dump = False
    self.dumppath = None
    self.radar_index_mapping = {}
  
  ##
  # Generate function. Basically same as calling compositing.generate but the pyarea is created from the
  # area definition.
  # @param dd: date
  # @param dt: time
  # @param tid: the area identifier (only used for identification purpose, actual area is taken from the area_definition)
  # @return a filename pointing to the tile
  def generate(self, dd, dt, tid):
    comp = compositing.compositing()
    comp.xscale = self.xscale
    comp.yscale = self.yscale
    comp.detectors = self.detectors
    comp.ignore_malfunc = self.ignore_malfunc
    comp.prodpar = self.prodpar
    comp.product = self.product
    comp.height = self.height
    comp.elangle = self.elangle
    comp.range = self.range
    comp.selection_method = self.selection_method
    comp.interpolation_method = self.interpolation_method
    comp.qitotal_field = self.qitotal_field
    comp.applygra = self.applygra
    comp.zr_A = self.zr_A
    comp.zr_b = self.zr_b
    comp.applygapfilling = self.applygapfilling
    comp.applyctfilter = self.applyctfilter
    comp.quantity = self.quantity
    comp.gain = self.gain
    comp.offset = self.offset    
    comp.minvalue = self.minvalue
    comp.filenames = self.filenames
    comp.verbose = self.verbose
    comp.reprocess_quality_field = self.reprocess_quality_field
    comp.dump = self.dump
    comp.dumppath = self.dumppath
    comp.radar_index_mapping = self.radar_index_mapping

    pyarea = _area.new()
    pyarea.id = "tiled area subset %s"%tid
    pyarea.xsize = self.area_definition.xsize
    pyarea.ysize = self.area_definition.ysize
    pyarea.xscale = self.area_definition.xscale
    pyarea.yscale = self.area_definition.yscale
    pyarea.extent = self.area_definition.extent
    pyarea.projection = _projection.new("dynamic pcsid", "dynamic pcs name", self.area_definition.pcsdef)    
    
    comp.store_dt_area_info(dd, dt, pyarea)
    
    logger.debug("tiled_compositing.py [%s]:generate(): Generating composite for tile_area=%s"%(comp.info_str, self.area_definition.id))
    result = comp.generate(dd, dt, pyarea)
    
    if result == None:
      logger.debug("tiled_compositing.py [%s]:generate(): No composite for tile_area=%s could be generated.", comp.info_str, self.area_definition.id)
      return (tid, None)

    fileno, outfile = rave_tempfile.mktemp(suffix='.h5', close="True")
  
    rio = _raveio.new()
    rio.object = result
    rio.filename = outfile
    rio.save()

    logger.info("tiled_compositing.py [%s]:generate(): Finished generating composite for tile_area=%s; file =%s", comp.info_str, self.area_definition.id, comp.info_str,outfile)
    return (tid, rio.filename)

##
# The actual compositing instance forwarding requests to the tilers.
class tiled_compositing(object):
  ##
  # Constructor
  # @param c: the compositing instance
  # @param preprocess_qc: If the ingoing files should be preprocessed or not before
  # they are sent to the tile generators (might improve performance in some cases)
  # @param mp_process_qc: If preprocess_qc is True, then this flag will send the
  # files to separate processes for quality control execution
  # @param mp_process_qc_split_evenly: Of mp_process_qc is True, then this flag will
  # indicate if the incomming files should be splitted evenly between the
  # processes. If false, then one file at a time are handled.
  #
  def __init__(self, c, preprocess_qc=False, mp_process_qc=False, mp_process_qc_split_evenly=False):
    self.compositing = c
    # If preprocess_qc = False, then the tile generators will take care of the preprocessing of the tiles
    # otherwise, the files will be qc-processed, written to disk and these filepaths will be sent to the
    # tile generators instead. Might or might not improve performance depending on file I/O etc..
    self.preprocess_qc = preprocess_qc
    self.mp_process_qc = mp_process_qc
    self.mp_process_qc_split_evenly = mp_process_qc_split_evenly
    self.verbose = c.verbose
    self.logger = logger
    self.file_objects = {}
    self.nodes = ""
    self.how_tasks = ""
    self.number_of_quality_control_processes = RAVE_QUALITY_CONTROL_PROCESSES
    self._do_remove_temporary_files=False
    self.compositing.area_id = 'None'
    self.compositing.dt_str = "_dt_unset_"
    self.BALTRAD_DEV_TEST_LOGGING = False # Disabled during production
    if os.getenv('BLT_RAVE_NIMBUS')=='1':
      # Executed only when the NIMBUS RAVE extension is enabled. Regular BALTRAD is NOT affected.
      # For testing activate by: export BLT_RAVE_NIMBUS_DEV_TEST_LOGGING=1
      # NOTE: self.compositing is the overall "tiled_compositing" object 
      #       instantiated from "rave_pgf_composite_plugin.py" as:
      #       comp = compositing(ravebdb); comp.filenames = files 
      self.BALTRAD_DEV_TEST_LOGGING = (os.getenv('BLT_RAVE_NIMBUS_DEV_TEST_LOGGING')=='1')
      #self.BALTRAD_DEV_TEST_LOGGING = True  # Activated only for testing
      self.nimbusQc_detected = False
      self.whichNimbus_qc = ""
      self.nimbusQc_tasks_main = []
      self.how_task_args_value_stored = ""
      self.nimbusQc_dict_qc_nodes = {}
      self.nimbusQc_dict_qc_nodes_str = ""
      self.nimbusQc_short_tasks_str = ""
      self.nimbusQc_satfilter_file = ""
      self.nimbusQc_tasks = []
      self.nimbusQc_tasks_str = "" # unused
      # self.nimbusQc_tasks_str = ','.join(self.nimbusQc_tasks)



  ##
  # Fetches the file objects and if self.preprocess_qc is True performs the quality controls.
  # If quality control processing is performed successfully, then temporary files are created
  # and their removal is required if self._do_remove_temporary_files = True 
  # @return (a dictionary with filenames as keys and objects as values, and a string with all included nodes names) 
  #
  def _fetch_file_objects(self):
    self.logger.info("tiled_compositing.py [%s]:_fetch_file_objects(): Fetching (and processing) %d files for tiled compositing"%(self.compositing.info_str, len(self.compositing.filenames)))

    result, nodes, how_tasks, all_files_malfunc = self.compositing.fetch_objects()
    if self.preprocess_qc:
      self._do_remove_temporary_files=False
      result, algorithm, qfields = self.compositing.quality_control_objects(result)
      try:
        result = store_temporary_files(result)
        self.compositing.filenames = result.keys()
        self._do_remove_temporary_files=True
      except Exception:
        self.logger.exception("tiled_compositing.py [%s]:_fetch_file_objects(): Failed to create temporary files. will not preprocess qc."%(self.compositing.info_str))

    self.logger.info("tiled_compositing.py [%s]:_fetch_file_objects(): Finished fetching (and processing) %d files for tiled compositing"%(self.compositing.info_str, len(self.compositing.filenames)))

    return (result, nodes, how_tasks, all_files_malfunc)

  ##
  # Fetches the file objects including the quality control by utilizing the multiprocessing
  # capabilities. I.e. instead of performing the fetching and quality controls within this
  # process. This job is spawned of to separate processors that manages the qc.
  # This will generate a number of temporary files that should be removed if self._do_remove_temporary_files=True
  # @return (a dictionary with filenames as keys and objects as values, and a string with all included nodes names)
  # 
  def _fetch_file_objects_mp(self):
    self.logger.info("tiled_compositing.py [%s]: MP Fetching (and processing) %d files for tiled compositing"%(self.compositing.info_str, len(self.compositing.filenames)))
    args = []
    ncpucores = multiprocessing.cpu_count()

    self._do_remove_temporary_files=False
    # We want to determine how many processes we are going to get prior
    # splitting the files
    #
    if self.mp_process_qc_split_evenly and self.number_of_quality_control_processes > 0:
      nobjects = len(self.compositing.filenames)
      nrfiles = len(self.compositing.filenames)
      nrprocesses = self.number_of_quality_control_processes
      if nrprocesses > ncpucores:
        nrprocesses = ncpucores
      if nrprocesses == ncpucores and ncpucores > 1:
        nrprocesses = nrprocesses - 1
      nrslices = nrfiles / nrprocesses
      for x in range(0, nrfiles, nrslices):
        args.append((self.compositing.filenames[x:x+nrslices], self.compositing.detectors, self.compositing.reprocess_quality_field, self.compositing.ignore_malfunc))
    else:
      for fname in self.compositing.filenames:
        args.append(([fname], self.compositing.detectors, self.compositing.reprocess_quality_field, self.compositing.ignore_malfunc))
      
    nobjects = len(args)
    nrprocesses = nobjects
    if nrprocesses > self.number_of_quality_control_processes:
      nrprocesses = self.number_of_quality_control_processes
    if nrprocesses > ncpucores:
      nrprocesses = ncpucores
    if nrprocesses == ncpucores and ncpucores > 1:
      nrprocesses = nrprocesses - 1 # We always want to leave at least one core for something else
    
    pool = multiprocessing.Pool(nrprocesses)
    
    results = [] # Storage for the result from the mp processes
    r = pool.map_async(execute_quality_control, args, callback=results.append)
    
    r.wait()
    pool.terminate()
    pool.join()
    
    filenames=[]
    for r in results[0]:
      filenames.extend(r[0])
      if r[1] == False:
        self.logger.info("tiled_compositing.py [%s]: MP quality control processing of %s failed."%(self.compositing.info_str, str(r[2])))
    
    self.compositing.filenames = filenames
    self._do_remove_temporary_files=True
    
    result, nodes, how_tasks, all_files_malfunc = self.compositing.fetch_objects()
    
    self.logger.info("tiled_compositing.py [%s]: MP Fetching (and processing) %d files for tiled compositing"%(self.compositing.info_str, len(self.compositing.filenames)))
    
    return (result, nodes, how_tasks, all_files_malfunc)


  ##
  # Creates the composite arguments that should be sent to one tiler.
  # @param adef: the area definition
  # @return the composite argument instance
  def _create_multi_composite_argument(self, adef=None):
    a = multi_composite_arguments()
    a.xscale = self.compositing.xscale
    a.yscale = self.compositing.yscale
    a.detectors = self.compositing.detectors
    a.ignore_malfunc = self.compositing.ignore_malfunc
    a.prodpar = self.compositing.prodpar
    a.product = self.compositing.product
    a.height = self.compositing.height
    a.elangle = self.compositing.elangle
    a.range = self.compositing.range
    a.selection_method = self.compositing.selection_method
    a.interpolation_method = self.compositing.interpolation_method
    a.qitotal_field = self.compositing.qitotal_field
    a.applygra = self.compositing.applygra
    a.zr_A = self.compositing.zr_A
    a.zr_b = self.compositing.zr_b
    a.applygapfilling = self.compositing.applygapfilling
    a.applyctfilter = self.compositing.applyctfilter
    a.quantity = self.compositing.quantity
    a.gain = self.compositing.gain
    a.offset = self.compositing.offset
    a.verbose = self.verbose
    a.dump = self.compositing.dump
    a.dumppath = self.compositing.dumppath
    a.reprocess_quality_field = self.compositing.reprocess_quality_field
    a.area_definition = adef
    
    return a
  
  ##
  # Creates an area definition used for passing to the tiler
  # @param pyarea: the python c area object
  # @return the area defintion
  def _create_tiled_area_definition(self, pyarea):
    return tiled_area_definition(pyarea.id, pyarea.projection.definition, pyarea.xscale, pyarea.yscale, pyarea.xsize, pyarea.ysize, pyarea.extent)
  
  ##
  # Creates the list of arguments to be sent to the tilers. Each item in the returned list is supposed to represent one tile
  # @param dd: date
  # @param dt: time
  # @param aid: the area id (that might or not be tiled)
  # @return a list of argument lists
  def _create_arguments(self, dd, dt, pyarea):
    tiled_areas = rave_tile_registry.get_tiled_areas(pyarea)
    
    args=[]
    
    for t in tiled_areas:
      mcomp = self._create_multi_composite_argument(self._create_tiled_area_definition(t))
      args.append([mcomp, dd, dt, t.id])
    
    # Now, make sure we have the correct files in the various areas
    self._add_files_to_argument_list(args, tiled_areas)

    # And add the radar index value to be used for each radar source so that each tile
    # have same information
    self._add_radar_index_value_to_argument_list(args)

    # We also must ensure that if any arg contains 0 files, there must be a date/time set
    if not self._ensure_date_and_time_on_args(args):
      raise Exception("tiled_compositing.py [%s]: Could not ensure existing date and time for composite"%(self.compositing.info_str))
        
    return args
  
  def _add_files_to_argument_list(self, args, tiled_areas):
    self.logger.info("tiled_compositing.py [%s]: Distributing polar objects among %d tiles"%(self.compositing.info_str, len(args)))

    # Loop through tile areas
    for i in range(len(tiled_areas)):
        p = tiled_areas[i].projection
        llx, lly, urx, ury = tiled_areas[i].extent

        # Loop through radars
        for k in self.file_objects.keys():
            v = self.file_objects[k]
            if not _polarscan.isPolarScan(v) and not _polarvolume.isPolarVolume(v):
                continue
            if _polarvolume.isPolarVolume(v):
                v = v.getScanWithMaxDistance()
            scan = v
            
            if self.compositing.quantity not in scan.getParameterNames():
                self.logger.info("tiled_compositing.py  [%s]: Quantity %s not in data from %s" % (self.compositing.info_str, self.compositing.quantity, scan.source))
                continue

            bi = scan.nbins - 1
            
            # Loop around the scan
            for ai in range(scan.nrays):
                lon, lat = scan.getLonLatFromIndex(bi, ai)
                x, y = p.fwd((lon, lat))
                
                # If this position is inside the tile, then add the radar's file string to the list and then bail
                if x >= llx and x <= urx and y >= lly and y <= ury:
                    if not k in args[i][0].filenames:
                        args[i][0].filenames.append(k)
                        break # No need to continue

    for idx in range(len(args)):
      self.logger.info("tiled_compositing.py [%s]: Tile %s contains %d files and dimensions %i x %i"%(self.compositing.info_str,
                        args[idx][0].area_definition.id, len(args[idx][0].filenames), args[idx][0].area_definition.xsize, args[idx][0].area_definition.ysize))
      
    self.logger.info("tiled_compositing.py [%s]: Finished splitting polar object"%(self.compositing.info_str))
  
  def _add_radar_index_value_to_argument_list(self, args):
    ctr = 1
    for k in self.file_objects.keys():
      v = self.file_objects[k]
      if not _polarscan.isPolarScan(v) and not _polarvolume.isPolarVolume(v):
        continue
      sourceid = v.source
      try:
        osource = odim_source.ODIM_Source(v.source)
        if osource.wmo:
          sourceid = "WMO:%s"%osource.wmo
        elif osource.rad:
          sourceid = "RAD:%s"%osource.rad
        elif osource.nod:
          sourceid = "NOD:%s"%osource.nod
      except:
        pass
            
      for arg in args:
        arg[0].radar_index_mapping[sourceid] = ctr
      ctr = ctr + 1        

  def _ensure_date_and_time_on_args(self, args):
    dtstr = None
    ddstr = None
    for k in self.file_objects.keys():
      v = self.file_objects[k]
      if not _polarscan.isPolarScan(v) and not _polarvolume.isPolarVolume(v):
        continue
      dtstr = v.time
      ddstr = v.date
      break
    
    if dtstr is None or ddstr is None:
      self.logger.warning("tiled_compositing.py [%s]: Could not determine any date and time string"%(self.compositing.info_str))
      return False
    
    for arg in args:
      if len(arg[0].filenames) == 0 and (arg[1] is None or arg[2] is None):
        arg[1] = ddstr
        arg[2] = dtstr
    
    return True
  
  def _create_lon_lat_extent(self, carg):
    pj = _projection.new("x", "y", carg.area_definition.pcsdef)
    lllon,lllat = pj.inv((carg.area_definition.extent[0], carg.area_definition.extent[1]))
    urlon,urlat = pj.inv((carg.area_definition.extent[2], carg.area_definition.extent[3]))
    return (lllon,lllat,urlon,urlat)

  ##
  # Checks if temporary (partial- / sub-) composite object is Nimbus QC.
  # These attibutes will be filled accordingly to the situation:
  #   self.nimbusQc_detected, self.whichNimbus_qc
  # @param obj: temporary composite object (CartesianVolumeCore)
  def eval_nimbus_qc_temporary_comp(self, obj, tile_area):
    # Executed only when the NIMBUS RAVE extension is enabled. Regular BALTRAD is NOT affected.
    # if os.getenv('BLT_RAVE_NIMBUS')=='1':
    #   self.eval_nimbus_qc_temporary_comp(obj, tile_area)
    # self.nimbusQc_detected = False # done in constructor
    # self.whichNimbus_qc = ""
    # NIMBUS: The QC-VOLUME objects should contain attribute ("how/task", "nimbus-qc") 
    # so that we can detect in BALTRAD that given hdf5 (radar) volume files
    # have been process with Quality-Controls.
    # Keep in mind that the applied quality controls can vary across radar-sources.
    # See config: toolbox/nimbus-rave/config/modules_by_radar.xml
    #
    # Example meta-data of COMPOSINTE from NIMBUS QC-PVOL:
    # (A) satfilter is applied 
    # group: how {
    #   :nodes = "detur,demem,deros,nldhl,nlhrw,deboo..."
    #   :task = "nimbus-qc-satfilter,ropo,beamb,qi-total,satfilter" ;
    #   :task_args = "satfilter_file:S_NWC_PC_MSG2_Europe-VISIR_20220127T100000Z.nc;
    #                 satfilter_nodes:detur,demem,deros,nldhl,...;
    #                 ropo_nodes:detur,demem,deros,nldhl,...;
    #                 beamb_nodes:detur,demem,deros,nldhl,...;..."
    # (B) satfilter NOT applied 
    # group: how {
    #   :nodes = "detur,demem,deros,nldhl,nlhrw,deboo..."
    #   :task = "nimbus-qc,ropo,beamb,qi-total" ;
    #   :task_args = "ropo_nodes:detur,demem,deros,nldhl,...;
    #                 beamb_nodes:detur,demem,deros,nldhl,...;..."

    attr_names = obj.getAttributeNames()
    # During compositing attr_names looks typically like this:
    # attr_names=['how/task', 'how/wavelength', 'what/version', 'how/scan_count', 'what/object', 'how/software', 'how/system', 'how/TXtype']
    
    # At this stage self.compositing.info_str ~= "dt=20220127T100000,area=odc_area"
    # Let's change it to a tile_area using 'tile_area'
    #     self.compositing.info_str_tile_area ~= "dt=20220127T100000,area=odc_area_se"
    try:
      self.compositing.info_str_tile_area = self.compositing.info_str.split(',area=')[0] + ',area={}'.format(tile_area)
    except:
      self.compositing.info_str_tile_area = self.compositing.info_str

    self.nimbusQc_detected = False
    if "how/task" not in attr_names:
      return

    # NOTE: obj == temporary composite object (CartesianVolumeCore)
    # The compositing(.py) object used for the composition process of this (single) tile is already gone at this stage.
    # So there is no access to the stored information in variables:
    # "nimbusQc_dict_qc_nodes", "nimbusQc_satfilter_file", or "whichNimbus_qc"
    # This information has to be reconstructed from: "how/task" and "how/task_args"

    how_task_value = obj.getAttribute("how/task")

    if not "nimbus-qc" in how_task_value:
      self.nimbusQc_short_tasks_str = how_task_value
    else:
      self.nimbusQc_detected = True
      self.whichNimbus_qc = "nimbus-qc"
      
    if not self.nimbusQc_detected:
      return

    if "nimbus-qc-satfilter" in how_task_value:
      self.whichNimbus_qc = "nimbus-qc-satfilter"

    how_task_args_value = obj.getAttribute("how/task_args")
    if self.BALTRAD_DEV_TEST_LOGGING:
      self.logger.info('tiled_compositing.py [{}]:eval_nimbus_qc(): how_task_value="{}"'.format(\
                      self.compositing.info_str_tile_area, how_task_value))
      self.logger.info('tiled_compositing.py [{}]:eval_nimbus_qc(): how_task_args_value="{}"'.format(\
                      self.compositing.info_str_tile_area, how_task_args_value))
      if self.whichNimbus_qc == "nimbus-qc-satfilter":
        self.logger.debug('tiled_compositing.py [{}]:eval_nimbus_qc(): (self.compositing) Before NIMBUS-QC tile-file: whichNimbus_qc="{}"; nimbusQc_short_tasks_str ="{}"; nimbusQc_satfilter_file="{}"; nimbusQc_dict_qc_nodes_str="{}"; nimbusQc_dict_qc_nodes="{}"; '.format(\
                          self.compositing.info_str_tile_area, \
                          self.compositing.whichNimbus_qc,\
                          self.compositing.nimbusQc_short_tasks_str, self.compositing.nimbusQc_satfilter_file,\
                          self.compositing.nimbusQc_dict_qc_nodes_str, self.compositing.nimbusQc_dict_qc_nodes ))
      else:
        self.logger.debug('tiled_compositing.py [{}]:eval_nimbus_qc(): (self.compositing) Before NIMBUS-QC tile-file: whichNimbus_qc="{}"; nimbusQc_short_tasks_str ="{}"; nimbusQc_dict_qc_nodes_str="{}"; nimbusQc_dict_qc_nodes="{}"; '.format(\
                          self.compositing.info_str_tile_area, \
                          self.compositing.whichNimbus_qc,\
                          self.compositing.nimbusQc_short_tasks_str,\
                          self.compositing.nimbusQc_dict_qc_nodes_str, self.compositing.nimbusQc_dict_qc_nodes ))
      if self.whichNimbus_qc == "nimbus-qc-satfilter":
        self.logger.debug('tiled_compositing.py [{}]:eval_nimbus_qc(): (self) Before NIMBUS-QC tile-file: whichNimbus_qc="{}"; nimbusQc_short_tasks_str ="{}"; nimbusQc_satfilter_file="{}"; nimbusQc_dict_qc_nodes_str="{}"; nimbusQc_dict_qc_nodes="{}"; '.format(\
                          self.compositing.info_str_tile_area, \
                          self.whichNimbus_qc,\
                          self.nimbusQc_short_tasks_str, self.nimbusQc_satfilter_file,\
                          self.nimbusQc_dict_qc_nodes_str, self.nimbusQc_dict_qc_nodes ))
      else:
        self.logger.debug('tiled_compositing.py [{}]:eval_nimbus_qc(): (self) Before NIMBUS-QC tile-file: whichNimbus_qc="{}"; nimbusQc_short_tasks_str ="{}"; nimbusQc_dict_qc_nodes_str="{}"; nimbusQc_dict_qc_nodes="{}"; '.format(\
                          self.compositing.info_str_tile_area, \
                          self.whichNimbus_qc,\
                          self.nimbusQc_short_tasks_str,\
                          self.nimbusQc_dict_qc_nodes_str, self.nimbusQc_dict_qc_nodes ))
    #
    # During ingestion and processing of a single radar in the tile this could look like this:
    # tiled_compositing.py [dt=20220127T100000,area=odc_area_nw]:eval_nimbus_qc(): 
    #                      how_task_value="nimbus-qc,ropo,beamb,qi-total"
    # tiled_compositing.py [dt=20220127T100000,area=odc_area_nw]:eval_nimbus_qc(): 
    #                      how_task_args_value="ropo_nodes:nlhrw;beamb_nodes:nlhrw;qi-total_nodes:nlhrw;"
    # NOTE: every (compositing)_tile has a different set of radar-sources (nodes) and different set of QC-filters applied.
    #
    # NOTE: self.compositing is the overall "compositing" object 
    #       instantiated from "rave_pgf_composite_plugin.py" as:
    #       comp = compositing(ravebdb); comp.filenames = files
    #       .. than re-instantiated as "tiled_compositing"
    #       comp = tiled_compositing(comp)

    try:
      # how_tasks_value .. is a string
      how_task_list = how_task_value.split(',')
    except:
      how_task_list = []

    # how_task_list += self.nimbusQc_tasks_main .. why to do this?
    for extra_qc_item in how_task_list:
      if 'nimbus-qc' in extra_qc_item:
        continue
      if extra_qc_item not in self.compositing.nimbusQc_short_tasks_str:
        if self.compositing.nimbusQc_short_tasks_str!="":
          self.compositing.nimbusQc_short_tasks_str +=","
        self.compositing.nimbusQc_short_tasks_str +="{}".format(extra_qc_item)

    try:
      how_task_args_list = how_task_args_value.split(';')
    except:
      how_task_args_list = []
    # There will be more radar-sources(/nodes) after <qc-filter-name>_nodes:<..>,<..>,<..>,
    # how_task_args_value = "ropo_nodes:nlhrw,deros;beamb_nodes:nlhrw,deros;qi-total_nodes:nlhrw,deros;"
    # ==>
    # how_task_args_list = ['ropo_nodes:nlhrw,deros', 'beamb_nodes:nlhrw,deros', 'qi-total_nodes:nlhrw,deros', '']
    # 
    # In case of satfilter applied there is extra "satfilter_file:<file-name>" in  'how/task_args'
    # how_task_args_value = "satfilter_file:S_NWC_PC_MSG2_Europe-VISIR_20220127T100000Z.nc;satfilter_nodes:<..>,<..>,<..>,;ropo_nodes:<..>,<..>,<..>,;..."
    # ==> 
    # how_task_args_list = ['satfilter_file:S_NWC_PC_MSG2_Europe-VISIR_20220127T100000Z.nc','satfilter_nodes:<..>,<..>,<..>','ropo_nodes:<..>,<..>,<..>', '']

    for qc_nodes_item in how_task_args_list:
      if ':' not in qc_nodes_item:
        continue
      if 'satfilter_file' in qc_nodes_item:
        try:
          self.compositing.nimbusQc_satfilter_file = qc_nodes_item.split(':')[1]
        except:
          pass
        continue
      # qc_nodes_item == '<qc-filter-name>_nodes:<..>,<..>,<..>'
      qc_nodes_list = qc_nodes_item.split(':')
      #try:
      qc_name = (qc_nodes_list[0]).split('_nodes')[0]  # example 'ropo_nodes'
      nodes_list_str = qc_nodes_list[1]  # example 'nlhrw,deros'
      nodes_list     = nodes_list_str.split(',')
      if qc_name not in self.compositing.nimbusQc_dict_qc_nodes:
        self.compositing.nimbusQc_dict_qc_nodes[qc_name] = []
      for _node in nodes_list:
        if _node not in self.compositing.nimbusQc_dict_qc_nodes[qc_name]:
          self.compositing.nimbusQc_dict_qc_nodes[qc_name].append(_node)
        
      #except:
      #  continue
      # Reconstruct every time the "string version" of "self.compositing.nimbusQc_dict_qc_nodes"
      nimbusQc_dict_qc_nodes_keys0 = list(self.compositing.nimbusQc_dict_qc_nodes.keys())
      if "satfilter" in nimbusQc_dict_qc_nodes_keys0:
        # make sure that satfilter nodes go 1st
        nimbusQc_dict_qc_nodes_keys = []
        nimbusQc_dict_qc_nodes_keys.append("satfilter")
        nimbusQc_dict_qc_nodes_keys0.remove("satfilter")
        nimbusQc_dict_qc_nodes_keys += nimbusQc_dict_qc_nodes_keys0
      else:
        nimbusQc_dict_qc_nodes_keys = nimbusQc_dict_qc_nodes_keys0
      self.compositing.nimbusQc_dict_qc_nodes_str = ""
      for qc_item in nimbusQc_dict_qc_nodes_keys:
        qc_nodes_list = self.compositing.nimbusQc_dict_qc_nodes[qc_item]
        qc_nodes_str = "{}_nodes:{};".format(qc_item, ','.join(qc_nodes_list))
        self.compositing.nimbusQc_dict_qc_nodes_str += qc_nodes_str

    self.compositing.whichNimbus_qc = self.whichNimbus_qc
    # if self.BALTRAD_DEV_TEST_LOGGING:
    #     self.logger.info('tiled_compositing.py [{}]:eval_nimbus_qc(): whichNimbus_qc="{}"; nimbusQc_short_tasks_str="{}"; nimbusQc_satfilter_file="{}"; nimbusQc_dict_qc_nodes_str="{}"; nimbusQc_dict_qc_nodes="{}"; '.format(\
    #                       self.compositing.info_str_tile_area, \
    #                       self.compositing.whichNimbus_qc,\
    #                       self.compositing.nimbusQc_short_tasks_str, self.compositing.nimbusQc_satfilter_file,\
    #                       self.compositing.nimbusQc_dict_qc_nodes_str, self.compositing.nimbusQc_dict_qc_nodes ))
    self.nimbusQc_short_tasks_str   = self.compositing.nimbusQc_short_tasks_str
    self.nimbusQc_satfilter_file    = self.compositing.nimbusQc_satfilter_file
    self.nimbusQc_dict_qc_nodes_str = self.compositing.nimbusQc_dict_qc_nodes_str
    self.nimbusQc_dict_qc_nodes     = self.compositing.nimbusQc_dict_qc_nodes

    if self.nimbusQc_detected:
      if self.whichNimbus_qc == "nimbus-qc-satfilter":
        self.logger.debug('tiled_compositing.py [{}]:eval_nimbus_qc(): Detected NIMBUS-QC tile-file: whichNimbus_qc="{}"; nimbusQc_short_tasks_str ="{}"; nimbusQc_satfilter_file="{}"; nimbusQc_dict_qc_nodes_str="{}"; nimbusQc_dict_qc_nodes="{}"; '.format(\
                          self.compositing.info_str_tile_area, \
                          self.compositing.whichNimbus_qc,\
                          self.compositing.nimbusQc_short_tasks_str, self.compositing.nimbusQc_satfilter_file,\
                          self.compositing.nimbusQc_dict_qc_nodes_str, self.compositing.nimbusQc_dict_qc_nodes ))
      else:
        self.logger.debug('tiled_compositing.py [{}]:eval_nimbus_qc(): Detected NIMBUS-QC tile-file: whichNimbus_qc="{}"; nimbusQc_short_tasks_str ="{}"; nimbusQc_dict_qc_nodes_str="{}"; nimbusQc_dict_qc_nodes="{}"; '.format(\
                          self.compositing.info_str_tile_area, \
                          self.compositing.whichNimbus_qc,\
                          self.compositing.nimbusQc_short_tasks_str,\
                          self.compositing.nimbusQc_dict_qc_nodes_str, self.compositing.nimbusQc_dict_qc_nodes ))
    else:
      self.logger.debug('tiled_compositing.py [{}]:eval_nimbus_qc(): NOT detected NIMBUS-QC tile-file'.format(self.compositing.info_str_tile_area))

  def store_dt_area_info(self,dd, dt, area):
    if area:
      try:
        self.compositing.area_id = area.id  # AreaCore object
      except:
        self.compositing.area_id = area    # string object
    else:
      self.compositing.area_id = 'None'
    try:
      if "tiled area subset " in self.compositing.area_id:
        self.compositing.area_id = (self.compositing.area_id[:]).split("tiled area subset ")[1]
    except:
      pass
    self.compositing.dt_str = "{}T{}".format(dd, dt)
    self.compositing.info_str = "dt={},area={}".format(self.compositing.dt_str,self.compositing.area_id)
  ##
  # Same as compositing generate but this is supposed to forward requests to a tiling mechanism
  # @param dd: date
  # @param dt: time
  # @param area: the area id
  def generate(self, dd, dt, area=None):
    pyarea = my_area_registry.getarea(area)
    
    self.store_dt_area_info(dd, dt, area)
    if self.preprocess_qc and self.mp_process_qc and self.number_of_quality_control_processes > 1:
      self.file_objects, self.nodes, self.how_tasks, all_files_malfunc = self._fetch_file_objects_mp()
    else:
      self.file_objects, self.nodes, self.how_tasks, all_files_malfunc = self._fetch_file_objects()
      
    if all_files_malfunc:
      self.logger.warning("tiled_compositing.py [{}]:generate(): Content of all provided files were marked as 'malfunc'. Since option 'ignore_malfunc' is set, no composite is generated!".format(self.compositing.info_str))
      return None

    args = self._create_arguments(dd, dt, pyarea)

    results = []
    
    ntiles = len(args)
    ncpucores = multiprocessing.cpu_count()

    nrprocesses = ntiles
    if not RAVE_TILE_COMPOSITING_PROCESSES is None:
       if nrprocesses > RAVE_TILE_COMPOSITING_PROCESSES:
         nrprocesses = RAVE_TILE_COMPOSITING_PROCESSES

    if nrprocesses > ncpucores:
      nrprocesses = ncpucores
    if nrprocesses == ncpucores and ncpucores > 1:
      nrprocesses = nrprocesses - 1 # We always want to leave at least one core for something else
    
    pool = multiprocessing.Pool(nrprocesses)
    
    r = pool.map_async(comp_generate, args, callback=results.append)
    
    r.wait()
    pool.terminate()
    pool.join()

    self.logger.info("tiled_compositing.py [%s]:generate(): Finished processing tiles, combining tiles."%(self.compositing.info_str))
    objects = []
    try:
      for v in results[0]:
        if v == None:
          tile_file = None # bugfix for the next line of code
          tile_area = area[:]
        else:
          tile_file = v[1]
          tile_area = v[0]
        if tile_file == None:
          self.logger.debug("tiled_compositing.py [%s]:generate(): NO partial composite for tile_area=%s was created. Therefore NOT included in complete composite."%(self.compositing.info_str, tile_area))
        else:
          self.logger.debug("tiled_compositing.py [%s]:generate(): _raveio.open(file=%s)"%(self.compositing.info_str, tile_file))
          o = _raveio.open(tile_file).object
          if os.getenv('BLT_RAVE_NIMBUS')=='1':
            # Executed only if the NIMBUS RAVE extension is enabled. Regular BALTRAD is NOT affected.
            self.eval_nimbus_qc_temporary_comp(o,tile_area)
          if _cartesianvolume.isCartesianVolume(o):
            o = o.getImage(0)
            o.objectType = _rave.Rave_ObjectType_COMP
          objects.append(o)
        
      t = _transform.new()

      self.logger.info("tiled_compositing.py [%s]:generate(): Combining %d tiles into one composite for area=%s."%(self.compositing.info_str, len(objects), area))

      result = t.combine_tiles(pyarea, objects)
      
      # Fix so that we get a valid place for /what/source and /how/nodes 
      result.source = "%s,CMT:%s"%(CENTER_ID,area)
      # NIMBUS change: self.nodes == "\'nldhl\',\'nlhrw\'"  ==> self.nodes == "nldhl,nlhrw"
      result.addAttribute('how/nodes', self.nodes)
 
      if self.how_tasks != "":
        result.addAttribute('how/task', self.how_tasks)

      self.logger.debug("tiled_compositing.py [%s]:generate(): Tiles combined"%(self.compositing.info_str))
      
      return result
    finally:
      if self._do_remove_temporary_files:
        for fname in self.compositing.filenames:
          try:
            os.unlink(fname)
          except:
            self.logger.warning("tiled_compositing.py [%s]:generate(): Failed to remove temporary file: %s"%(self.compositing.info_str, fname))
      
      if results != None:
        for v in results[0]:
          if v != None and v[1] != None and os.path.exists(v[1]):
            try:
              os.unlink(v[1])
            except Exception:
              self.logger.exception("tiled_compositing.py [%s]:generate(): Failed to unlink %s"%(self.compositing.info_str, v[1]))
    return None

##
# Function that handles the multiprocessing call for the multiprocessing
# @param args: tuple of 4 args (multi_composite_arguments, date, time, area identifier)
# @return result of multi_composite_arguments.generate
#
def comp_generate(args):
  try:
    return args[0].generate(args[1], args[2], args[3])
  except Exception:
    logger.exception("tiled_compositing.py:comp_generate(): Failed to call composite generator in tiler")
  return None

##
# Handles the multiprocessing call for the quality control section
# @param args: tuple of 4 args, ([filenames],[detectors], reprocess_quality_field, ignore_malfunc)
# @return a tuple of ([filenames], <execution status as boolean>, "filenames or source names")
def execute_quality_control(args):
  filenames,detectors,reprocess_quality_field,ignore_malfunc = args
  result = ([], False, "%s"%str(filenames))
  try:
    comp = compositing.compositing()
    comp.filenames.extend(filenames)
    comp.detectors.extend(detectors)
    comp.reprocess_quality_field = reprocess_quality_field
    comp.ignore_malfunc = ignore_malfunc
    
    objects, nodes, how_tasks, all_files_malfunc = comp.fetch_objects()
    
    objects, algorithm = comp.quality_control_objects(objects)

    status = True
    try:
      objects = store_temporary_files(objects)
    except Exception:
      status = False
    result = (objects.keys(), status, nodes)
  except Exception:
    logger.exception("tiled_compositing.py [%s]:execute_quality_control(): Failed to run quality control"%(comp.info_str))
  return result

if __name__=="__main__":
  comp = compositing.compositing()
  comp.filenames=["/projects/baltrad/rave/test/pytest/fixtures/pvol_seang_20090501T120000Z.h5",
                  "/projects/baltrad/rave/test/pytest/fixtures/pvol_searl_20090501T120000Z.h5",
                  "/projects/baltrad/rave/test/pytest/fixtures/pvol_sease_20090501T120000Z.h5",
                  "/projects/baltrad/rave/test/pytest/fixtures/pvol_sehud_20090501T120000Z.h5",
                  "/projects/baltrad/rave/test/pytest/fixtures/pvol_sekir_20090501T120000Z.h5",
                  "/projects/baltrad/rave/test/pytest/fixtures/pvol_sekkr_20090501T120000Z.h5",
                  "/projects/baltrad/rave/test/pytest/fixtures/pvol_selek_20090501T120000Z.h5",
                  "/projects/baltrad/rave/test/pytest/fixtures/pvol_selul_20090501T120000Z.h5",
                  "/projects/baltrad/rave/test/pytest/fixtures/pvol_seosu_20090501T120000Z.h5",
                  "/projects/baltrad/rave/test/pytest/fixtures/pvol_seovi_20090501T120000Z.h5",
                  "/projects/baltrad/rave/test/pytest/fixtures/pvol_sevar_20090501T120000Z.h5",
                  "/projects/baltrad/rave/test/pytest/fixtures/pvol_sevil_20090501T120000Z.h5"]
  
  comp.filenames=["/projects/baltrad/baltrad-test/fixtures4/seang_pvol_20140220T1300Z.h5",
                  "/projects/baltrad/baltrad-test/fixtures4/searl_pvol_20140220T1300Z.h5",
                  "/projects/baltrad/baltrad-test/fixtures4/sease_pvol_20140220T1300Z.h5",
                  "/projects/baltrad/baltrad-test/fixtures4/sehud_pvol_20140220T1300Z.h5",
                  "/projects/baltrad/baltrad-test/fixtures4/sekir_pvol_20140220T1300Z.h5",
                  "/projects/baltrad/baltrad-test/fixtures4/sekkr_pvol_20140220T1300Z.h5",
                  "/projects/baltrad/baltrad-test/fixtures4/selek_pvol_20140220T1300Z.h5",
                  "/projects/baltrad/baltrad-test/fixtures4/selul_pvol_20140220T1300Z.h5",
                  "/projects/baltrad/baltrad-test/fixtures4/seosu_pvol_20140220T1300Z.h5",
                  "/projects/baltrad/baltrad-test/fixtures4/seovi_pvol_20140220T1300Z.h5",
                  "/projects/baltrad/baltrad-test/fixtures4/sevar_pvol_20140220T1300Z.h5",
                  "/projects/baltrad/baltrad-test/fixtures4/sevil_pvol_20140220T1300Z.h5"]
  
  #, preprocess_qc=False, mp_process_qc=False, mp_process_qc_split_evenly=False
  tc = tiled_compositing(comp,True, True, True)
  #tc.generate("20090501","120000", "swegmaps_2000")
  tc.generate("20090501","120000", "bltgmaps_4000")
  #execute_tiled_compositing(comp, my_area_registry.getarea("swegmaps_2000"))
