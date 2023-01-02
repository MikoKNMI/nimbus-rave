'''
Copyright (C) 2010- Swedish Meteorological and Hydrological Institute (SMHI)

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
## Plugin for applying quality control on a volume that is initiated from the beast
## framework.
## Register in the RAVE PGF with: % pgf_registry -a -H http://<host>:<port>/RAVE
## --name=eu.baltrad.beast.applyqc_plugin --strings=source,date,time,anomaly-qc -m rave_pgf_apply_qc -f generate
## -d 'Apply quality controls on a polar volume'
##

## @file
## @author Mats Vernersson, SMHI
## @date 2016-04-15
## @author Michal Koutek, KNMI (Fixes needed for OPERA NIMBUS)
## @date 2022-01-06

import _raveio
import string
import os
import rave_tempfile
import rave_pgf_quality_registry
import rave_util
from rave_quality_plugin import QUALITY_CONTROL_MODE_ANALYZE_AND_APPLY
from rave_defines import RAVE_IO_DEFAULT_VERSION
import modules_selection
import odim_source
from inspect import signature

import rave_pgf_logger
logger = rave_pgf_logger.create_logger()

ravebdb = None
try:
  import rave_bdb
  ravebdb = rave_bdb.rave_bdb()
except:
  pass

## Creates a dictionary from a rave argument list
#@param arglist the argument list
#@return a dictionary
def arglist2dict(arglist):
  result={}
  for i in range(0, len(arglist), 2):
    result[arglist[i]] = arglist[i+1]
  return result

##
# Performs a quality control sequence on a volume
# @param volume: the volume to perform the quality controls on
# @param detectors: the detectors that should be run on the volume
#
def perform_quality_control(volume, detectors, qc_mode):
  for d in detectors:
    p = rave_pgf_quality_registry.get_plugin(d)
    if p != None:
      logger.debug("rave_pgf_apply_qc_plugin.py:: Processing volume with quality plugin %s. QC-mode: %s", d, qc_mode)
      volume = p.process(volume, True, qc_mode)
      if isinstance(volume,tuple):
        volume, _alg_ = volume[0],volume[1]
  return volume


##
# Performs the NIMBUS Quality Control sequence on a volume
# @param volume: the volume to perform the quality controls on
# @param detectors: the detectors that should be run on the volume
#
# Keep in mind that the applied quality controls can vary across (Opera) radar-sources.
# When detectors list (*) contains "nimbus-configured-qc" than we taken the qc-controls (qc-filters)
# from the config file: toolbox/nimbus-rave/config/modules_by_radar.xml
#
# In this case the "qc-control.groovy" script and the "QC-VOL-producer" groovy-route
# explicitely states: "nimbus-configured-qc" in the detectors list.
#
# When detectors list (*) contains "nimbus-satfilter-skip"
# then we will SKIP application of the satfilter on ALL radar-sources.
#
# When detectors list (*) contains "nimbus-satfilter-apply-all"
# then we will APPLY the satfilter on ALL radar-sources.
#
# (*) from a "qc-control.groovy" script
#
# See implementation and application:
#     config/qc-control.groovy ;     config-test/qc-control.groovy ;
#     config/NIMBUS_CFG_routes.json; config-test/NIMBUS_CFG_routes.json;
def perform_quality_control_nimbus(volume, detectors, qc_mode):
  # For testing activate by: export BLT_RAVE_NIMBUS_DEV_TEST_LOGGING=1
  BALTRAD_DEV_TEST_LOGGING = (os.getenv('BLT_RAVE_NIMBUS_DEV_TEST_LOGGING')=='1')

  node = odim_source.NODfromSource(volume)
  info_str = "dt={}T{},node={}".format(volume.date, volume.time, node)

  # Check if the list (detectors) contains "nimbus-configured-qc"
  if not "nimbus-configured-qc" in detectors:
    return perform_quality_control(volume, detectors, qc_mode), detectors

  # This part allows to choose the "modules" (alias "quality controls", alias "quality filters", alias "quality detectors") 
  # independently for each radar, as specified in from the config file: toolbox/nimbus-rave/config/modules_by_radar.xml
  modules = modules_selection.getModules(volume)
  if BALTRAD_DEV_TEST_LOGGING:
    logger.debug("rave_pgf_apply_qc_plugin.py: type(modules)={}, modules={}".format(type(modules),modules))

  if not modules:
    radar_specific_detectors_list = []
  else:
    radar_specific_detectors_list = modules.split(',')
    if "nimbus-satfilter-apply-all" in detectors:
      # In this case we will add "satfilter" to the radar_specific_detectors_list.
      logger.debug("rave_pgf_apply_qc_plugin.py [{}]: Processing volume: nimbus-satfilter-apply-all - adding satfilter to detectors: {}".format(info_str, radar_specific_detectors_list))
      if "satfilter" not in radar_specific_detectors_list:
        radar_specific_detectors_list.append("satfilter")
    elif "nimbus-satfilter-skip" in detectors:
      # In this case we will remove "satfilter" from the radar_specific_detectors_list.
      logger.debug("rave_pgf_apply_qc_plugin.py [{}]: Processing volume: nimbus-satfilter-skip - removing satfilter from detectors: {}".format(info_str, radar_specific_detectors_list))
      if "satfilter" in radar_specific_detectors_list:
        radar_specific_detectors_list.remove("satfilter")

    attr_names = volume.getAttributeNames()
    if "how/task" in attr_names:
      current_how_task_value0 = volume.getAttribute("how/task")
      # In NIMBUS the Nimbus-QC-PVOL ["how/task"] contains list of applied QC filters.
      # There is NO need to re-apply QC filters that have been already applied.
      current_how_task_list0 = current_how_task_value0.split(',')
      if len(current_how_task_list0) > 0:
        first_qc_item = current_how_task_list0[0][:]
        if "nimbus-qc" in first_qc_item:
          # remove the first_qc_item
          current_how_task_list0.remove(first_qc_item)
      # current_how_task_list will contain ONLY "how/task"-QC-items from existing registerd quality controls
      current_how_task_list = []
      for qc_item in current_how_task_list0:
        p = rave_pgf_quality_registry.get_plugin(qc_item)
        if p == None:
          # Remove (do not include) qc-items from "how/task" which are not really registered quality controls of Baltrad)
          continue
        else:
          # Add qc-items from "how/task" ONLY from existing registerd quality controls
          current_how_task_list.append(qc_item)
      for _qc in current_how_task_list:
        if _qc in radar_specific_detectors_list:
          logger.debug("rave_pgf_apply_qc_plugin.py [{}]: Processing volume: skipping {} from detectors; this QC already applied.".format(info_str, _qc))
          radar_specific_detectors_list.remove(_qc)
      current_detectors = ",".join(current_how_task_list)
      if len(current_how_task_list)>0:
        current_detectors += ","
    else:
      current_detectors = ""

  if radar_specific_detectors_list:
    logger.debug("rave_pgf_apply_qc_plugin.py [{}]: radar_specific_detectors_list={}".format(info_str, radar_specific_detectors_list))
    return_detectors = ",".join(radar_specific_detectors_list)
  else:
    logger.warning("rave_pgf_apply_qc_plugin.py [{}]: radar_specific_detectors_list=EMPTY".format(info_str))
    return_detectors = ""

  for d in radar_specific_detectors_list:
    p = rave_pgf_quality_registry.get_plugin(d)
    if p != None:
      # Before calling p.process(...) CHECK if the qc_plugin.process function has attribute "arguments" in definition.
      # arguments .. can be used to pass additional (logging) information into the process function.
      plugin_process_funct_sign = "{}".format(signature(p.process))
      # (self, obj, reprocess_quality_flag=True, quality_control_mode='analyze_and_apply', arguments=None)"
      if "arguments" not in plugin_process_funct_sign:
        logger.debug("rave_pgf_apply_qc_plugin.py [{}]: Processing volume with quality plugin {}.".format(info_str, d))
        volume = p.process(volume, True, qc_mode)
      else:
        arguments = {"info_str":info_str[:]}
        logger.debug("rave_pgf_apply_qc_plugin.py [{}]: Processing volume with quality plugin {}.".format(info_str, d))
        #logger.debug("rave_pgf_apply_qc_plugin.py [{}]: Processing volume with quality plugin {} with arguments={}.".format(info_str, d, arguments))
        volume = p.process(volume, True, qc_mode, arguments)
      if isinstance(volume,tuple):
        volume, _alg_ = volume[0],volume[1]
    else:
      logger.debug("rave_pgf_apply_qc_plugin.py [{}]: Processing volume: Could not find {} plugin.".format(info_str, d))
  if not "nimbus-satfilter-skip" in detectors:
    '''
    We are adding attribute ("how/task", "nimbus-qc")
    so that we can detect in BALTRAD routers that given hdf5 (radar) volume file
    has been process with Quality-Controls as specified for NIMBUS:

    { "type": "attr", "valueType": "STRING", "negated": false,"operator": "EQ",
      "attribute": "how/task", "value": "nimbus-qc"}
    '''
    
    if "nimbus-satfilter-apply-all" in detectors:
      #volume.addAttribute("how/task", "nimbus-qc-satfilter,*")
      return_detectors = "nimbus-qc-satfilter," + current_detectors + return_detectors[:]
    else:
      #volume.addAttribute("how/task", "nimbus-qc,*")
      return_detectors = "nimbus-qc," + current_detectors + return_detectors[:]
    volume.addAttribute("how/task", return_detectors)
    logger.info('rave_pgf_apply_qc_plugin.py [{}]: Processing volume: volume.addAttribute("how/task", "{}")'.format(info_str, return_detectors))
  else:
    '''
    If (for the testing) we skipped satfilter for all radar-sources
    then we are adding attribute ("how/task", "nimbus-qc-no-satfilter").
    In this case the related BALTRAD route has to check:

    { "type": "attr", "valueType": "STRING", "negated": false,"operator": "EQ",
      "attribute": "how/task", "value": "nimbus-qc-no-satfilter"}
    '''
    #volume.addAttribute("how/task", "nimbus-qc-no-satfilter,*")
    return_detectors = "nimbus-qc-no-satfilter," + current_detectors + return_detectors[:]
    volume.addAttribute("how/task", return_detectors)
    logger.info('rave_pgf_apply_qc_plugin.py [{}]: Processing volume: volume.addAttribute("how/task", "{}")'.format(info_str, return_detectors))
  
  return volume, return_detectors.split(',')


def generate_new_volume_with_qc(original_file, args):
  BALTRAD_DEV_TEST_LOGGING = (os.getenv('BLT_RAVE_NIMBUS_DEV_TEST_LOGGING')=='1')
  if BALTRAD_DEV_TEST_LOGGING:
    logger.debug("rave_pgf_apply_qc_plugin.py:generate_qc(): Generating new volume with quality controls from file={}; args={};".format(original_file,args))

  if ravebdb != None:
    #if BALTRAD_DEV_TEST_LOGGING:
    #  logger.debug("rave_pgf_apply_qc_plugin.py:generate_qc(): ravebdb.get_rave_object({}})".format(original_file))
    #volume = ravebdb.get_rave_object(original_file)
    # def get_rave_object(self, fname, lazy_loading=False, preloadedQuantities=None):
    # Alternative way of reading the file:
    fullpath_filename = ravebdb.path_from_uuid(original_file)
    if BALTRAD_DEV_TEST_LOGGING:
      logger.debug("rave_pgf_apply_qc_plugin.py:generate_qc(): _raveio.open({}).object".format(fullpath_filename))
    volume = _raveio.open(fullpath_filename).object
  else:
    fullpath_filename = original_file
    if BALTRAD_DEV_TEST_LOGGING:
      logger.debug("rave_pgf_apply_qc_plugin.py:generate_qc():  _raveio.open({}}).object".format(original_file))
    volume = _raveio.open(original_file).object

  if volume == None:
    logger.warning("rave_pgf_apply_qc_plugin.py:generate_qc(): Error reading file={}; volume is EMPTY; args={};".format(fullpath_filename,args))
    return None

  node = odim_source.NODfromSource(volume)
  info_str = "dt={}T{},node={}".format(volume.date, volume.time, node)

  logger.info("rave_pgf_apply_qc_plugin.py [{}]:generate_qc(): Generating new volume with quality controls from file={}; args={};".format(info_str,fullpath_filename,args))

  if "remove-malfunc" in args.keys():
    try:
      if args["remove-malfunc"].lower() in ["true", "yes", "y", "1"]:
        if BALTRAD_DEV_TEST_LOGGING:
          logger.debug("rave_pgf_apply_qc_plugin.py [{}]:generate_qc(): Checking volume for malfunc tags. Will remove scans, or complete volume, if marked malfunc.".format(info_str))
        volume = rave_util.remove_malfunc(volume)
        if volume == None:
          logger.warning("rave_pgf_apply_qc_plugin.py [{}]:generate_qc(): Malfunc volume! Since option 'remove_malfunc' is set, no new volume with QC applied will be generated! original_file={}; args={};".format(info_str,original_file,args))
          return None
    except:
      pass
  
  if "anomaly-qc" in args.keys():
    detectors = args["anomaly-qc"].split(",")
  else:
    detectors = []
    
  quality_control_mode = QUALITY_CONTROL_MODE_ANALYZE_AND_APPLY
  if "qc-mode" in args.keys():
    quality_control_mode = args["qc-mode"]

  
  if os.getenv('BLT_RAVE_NIMBUS')=='1':
    # Executed only when the NIMBUS RAVE extension is enabled. Regular BALTRAD is NOT affected.
    # NIMBUS NOTE:
    # In NIMBUS the rave-plugin execution of quality-controls is done via: 'rave_pgf_apply_qc_plugin.py'
    # BUT (!):
    # 'odc_polarQC.py' is executed ONLY from 'odc_generate.py' 
    # via command-line execution of 'odc_toolbox'.
    # For example:
    # /usr/lib/rave/bin/odc_toolbox -v -i $DIR_IN -o $DIR_OUT -q 'satfilter'
    # odc_polarQC.py is NOT used during NIMBUS processing.
    volume, detectors = perform_quality_control_nimbus(volume, detectors, quality_control_mode)
  else:
    # Regular BALTRAD; NO NIMBUS.
    volume = perform_quality_control(volume, detectors, quality_control_mode)

  
  new_time = args.get('time')
  if new_time:
    volume.time = new_time
    
  new_date = args.get('date')
  if new_date:
    volume.date = new_date

  if os.getenv('BLT_RAVE_NIMBUS')=='1':
    # Executed only when the NIMBUS RAVE extension is enabled. Regular BALTRAD is NOT affected.
    # NIMBUS NOTE: The below code will patch the date+time of the generated QC-PVOL.
    # Solution for BALTRAD issue: BALTRAD will NOT ingest 2x QC-PVOLs with same date.
    # You need this when you want to let BALTRAD produce 1st QC-PVOL without "satfilter"
    # and 2nd QC-PVOL with "satfilter".
    # We ADD 1 second to be able to store extra PVOL file (typically for QC-PVOL with satfilter).
    if "anomaly-qc" in args.keys():
      anomaly_qc = args.get('anomaly-qc')
    else:
      anomaly_qc = ""
    # NOTE: For testing the effects of satfilter there are 2 extra qc-options: nimbus-satfilter-skip, nimbus-satfilter-apply-all
    # Example args look like:
    # args={'date': '20220127', 'time': '100000','anomaly-qc': 'nimbus-configured-qc,nimbus-satfilter-skip', ...
    # args={'date': '20220127', 'time': '100000','anomaly-qc': 'nimbus-configured-qc,nimbus-satfilter-apply-all', ...
    if "nimbus-satfilter-apply-all" in anomaly_qc:
      if not new_time:
        logger.warning("rave_pgf_apply_qc_plugin.py [%s]:generate_qc(): Error changing time of PVOL (should not happen; new_time=='').",info_str)
      else:
        # volume.time .. ADD 1 second to be able to store extra PVOL file with different QC (to trick BALTRAD DB storage)
        # new_time = '123400
        seconds = int(new_time[4:])+1
        if seconds>59:
          seconds -= 2
        changed_time = new_time[:4]+str(seconds).zfill(2)
        volume.time  = changed_time
        info_str = "dt={}T{},node={}".format(volume.date, volume.time, node)
        logger.info("rave_pgf_apply_qc_plugin.py [%s]:generate_qc(): Changing time %s => %s",info_str,new_time,changed_time)

  logger.info("rave_pgf_apply_qc_plugin.py [%s]:generate_qc(): Quality controls applied on new volume: %s",info_str ,(",".join(detectors)))
  return volume
      
## Creates a new volume based on the incoming with quality controls applied to it
#@param files a list of files to apply quality controls on. currently assume only one file 
#@param arguments the arguments defining what quality controls to apply
#@return a temporary h5 file with the volume
def generate(files, arguments):
  args = arglist2dict(arguments)
  
  logger.debug("rave_pgf_apply_qc_plugin.py:generate(): rave_pgf_apply_qc_plugin called with arguments: %s", args)

  # should only be one file
  fname = files[0]


  logger.debug("rave_pgf_apply_qc_plugin.py:generate(): generate_new_volume_with_qc({},{})".format(fname, args))
  # We could extra args via args {}
  # args["extra"] = ...
  # TBD: How to pass the current  "area" and "tile-area" when pluging called during COMPOSITE processing?
  volume = generate_new_volume_with_qc(fname, args)
  
  if volume == None:
    logger.warning("rave_pgf_apply_qc_plugin.py:generate(): No volume with QC applied could be generated!")
    return None

  node = odim_source.NODfromSource(volume)
  info_str = "dt={}T{},node={}".format(volume.date, volume.time, node)

  _, outfile = rave_tempfile.mktemp(suffix='.h5', close="True")
  
  ios = _raveio.new()
  ios.object = volume
  ios.filename = outfile
  ios.version = RAVE_IO_DEFAULT_VERSION
  ios.save()
  
  logger.info("rave_pgf_apply_qc_plugin.py [{}]:generate(): Generated new volume with QC applied: {}".format(info_str, outfile))
  return outfile
  
