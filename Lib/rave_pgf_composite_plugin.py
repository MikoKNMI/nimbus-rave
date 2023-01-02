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
## Plugin for generating a composite that is initiated from the beast
## framework.
## Register in pgf with
## --name=eu.baltrad.beast.generatecomposite
## --strings=area,quantity,method,date,time,selection,anomaly-qc
## --floats=height -m rave_pgf_composite_plugin -f generate
##

## @file
## @author Anders Henja, SMHI
## @date 2010-10-15
## @author Michal Koutek, KNMI (Fixes needed for OPERA NIMBUS)
## @date 2022-01-06

import string
import rave_tempfile
import rave_pgf_logger
import _raveio
import rave_tile_registry
import os

from tiled_compositing import tiled_compositing
from compositing import compositing

from rave_defines import GAIN, OFFSET, RAVE_PGF_QUALITY_FIELD_REPROCESSING
from rave_defines import RAVE_IO_DEFAULT_VERSION

USE_AZIMUTHAL_NAVIGATION=True
try:
  from rave_defines import RAVE_PGF_AZIMUTHAL_NAVIGATION
  USE_AZIMUTHAL_NAVIGATION = RAVE_PGF_AZIMUTHAL_NAVIGATION
except:
  pass

USE_LAZY_LOADING=False
USE_LAZY_LOADING_PRELOADS=False
try:
  from rave_defines import RAVE_PGF_COMPOSITING_USE_LAZY_LOADING
  USE_LAZY_LOADING = RAVE_PGF_COMPOSITING_USE_LAZY_LOADING
except:
  pass

try:
  from rave_defines import RAVE_PGF_COMPOSITING_USE_LAZY_LOADING_PRELOADS
  USE_LAZY_LOADING_PRELOADS = RAVE_PGF_COMPOSITING_USE_LAZY_LOADING_PRELOADS
except:
  pass


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
# Converts a string into a number, either int or float
# @param sval the string to translate
# @return the translated value
# @throws ValueError if value not could be translated
#
def strToNumber(sval):
  try:
    return float(sval)
  except ValueError:
    return int(sval)


## Creates a composite
#@param files the list of files to be used for generating the composite
#@param arguments the arguments defining the composite
#@return a temporary h5 file with the composite
def generate(files, arguments):
  args = arglist2dict(arguments)
  
  comp = compositing(ravebdb)
  comp.filenames = files
  
  if "anomaly-qc" in args.keys():
    comp.detectors = args["anomaly-qc"].split(",")

  if "qc-mode" in args.keys():
    comp.set_quality_control_mode_from_string(args["qc-mode"])

  if "ignore-malfunc" in args.keys():
    try:
      if args["ignore-malfunc"].lower() in ["true", "yes", "y", "1"]:
        comp.ignore_malfunc = True
    except:
      pass

  comp.quantity = "DBZH"
  if "quantity" in args.keys():
    comp.quantity = args["quantity"]
  comp.gain = GAIN
  comp.offset = OFFSET

  comp.set_product_from_string("pcappi")
  if "method" in args.keys():
    comp.set_product_from_string(args["method"].lower())

  comp.height = 1000.0
  comp.elangle = 0.0
  comp.range = 200000.0

  if "prodpar" in args.keys():
    comp.prodpar = args["prodpar"]

  if "range" in args.keys() and comp.product == _rave.Rave_ProductType_PMAX:
    comp.range = float(args["range"])
  
  if "selection" in args.keys():
    comp.set_method_from_string(args["selection"])
    
  if "interpolation_method" in args.keys():
    interpolation_method = args["interpolation_method"].upper()
    comp.set_interpolation_method_from_string(interpolation_method)
    if interpolation_method != "NEAREST_VALUE":
      if comp.quantity == "DBZH":
        comp.minvalue = -30.0 # use a minimum value of -30.0 for DBZH
      else:
        # interpolation for other quantities than DBZH has not been tested, therefore not yet considered 
        # supported. Should in theory work, but the setting of minvalue must be adjusted for the quantity
        logger.info("rave_pgf_composite_plugin.py:generate(): Interpolation method %s is currently only supported for quantity DBZH. Provided quantity: %s. No composite generated." % (interpolation_method, comp.quantity))
        return None
  
  if "qitotal_field" in args.keys():
    comp.qitotal_field = args["qitotal_field"]
    
  if "reprocess_qfields" in args.keys():
    comp.reprocess_quality_field = args["reprocess_qfields"]
  else:
    comp.reprocess_quality_field = RAVE_PGF_QUALITY_FIELD_REPROCESSING
  
  comp.use_azimuthal_nav_information = USE_AZIMUTHAL_NAVIGATION

  comp.use_lazy_loading = USE_LAZY_LOADING
  comp.use_lazy_loading_preloads = USE_LAZY_LOADING_PRELOADS

  # We do not support best fit compositing right now
  #comp.pcsid = options.pcsid
  #comp.xscale = options.scale
  #comp.yscale = options.scale
  
  #if options.gf: Activate gap filling for rule
  #  comp.applygapfilling = True
  
  # Optional cloud-type residual non-precip filter
  if "ctfilter" in args:
    if eval(args["ctfilter"]):
      comp.applyctfilter = True
  
  if "applygra" in args:
    comp.applygra = True
  if "zrA" in args:
    comp.zr_A = float(args["zrA"])
  if "zrb" in args:
    comp.zr_b = float(args["zrb"])
  
  
  if rave_tile_registry.has_tiled_area(args["area"]):
    # NOTE: The "comp object" has been instantiated first as: comp = compositing(ravebdb)
    comp = tiled_compositing(comp)
  
  result = comp.generate(args["date"], args["time"], args["area"])

  try:
    area_id = args["area"]
  except:
    area_id = 'None'
  try:
    dt_str = "{}T{}".format(args["date"], args["time"])
  except:
    dt_str = "_dt_unset_"
  info_str = "dt={},area={}".format(dt_str,area_id)

  if result == None:
    logger.warning("rave_pgf_composite_plugin.py [{}]:generate(): No composite could be generated.".format(info_str))
    return None
  
  if os.getenv('BLT_RAVE_NIMBUS')=='1':
    # Executed only when the NIMBUS RAVE extension is enabled. Regular BALTRAD is NOT affected.
    # For testing activate by: export BLT_RAVE_NIMBUS_DEV_TEST_LOGGING=1
    BALTRAD_DEV_TEST_LOGGING = (os.getenv('BLT_RAVE_NIMBUS_DEV_TEST_LOGGING')=='1')
    #BALTRAD_DEV_TEST_LOGGING = True  # Activated only for testing
    '''
    BALTRAD DEX: Detect FULLY nimbus-qc composites:
    { "type": "attr", "valueType": "STRING", "negated": true,"operator": "LIKE",
      "attribute": "how/task", "value": "nimbus-qc-satfilter,*"},
    { "type": "attr", "valueType": "STRING", "negated": false,"operator": "LIKE",
      "attribute": "how/task", "value": "nimbus-qc,*"},
    - or -
    BALTRAD DEX: Detect nimbus-qc composites (explicitely) WITHOUT satfilter:
    { "type": "attr", "valueType": "STRING", "negated": false,"operator": "LIKE",
      "attribute": "how/task", "value": "nimbus-qc-no-satfilter*"},
    '''
    if comp.nimbusQc_detected:
      if BALTRAD_DEV_TEST_LOGGING:
        logger.debug("rave_pgf_composite_plugin.py [{}]:generate(): Detected NIMBUS_QC files used for this composite; whichNimbus_qc='{}'".format(info_str, comp.whichNimbus_qc))
        try:
          logger.debug("rave_pgf_composite_plugin.py [{}]:generate(): (tiled-)compositing.whichNimbus_qc ='{}')".format(info_str, comp.whichNimbus_qc))
          logger.debug("rave_pgf_composite_plugin.py [{}]:generate(): (tiled-)compositing.nimbusQc_short_tasks_str ='{}')".format(info_str, comp.nimbusQc_short_tasks_str))
          logger.debug("rave_pgf_composite_plugin.py [{}]:generate(): (tiled-)compositing.nimbusQc_satfilter_file ='{}')".format(info_str, comp.nimbusQc_satfilter_file))
          #logger.debug("rave_pgf_composite_plugin.py [{}]:generate(): (tiled-)compositing.nimbusQc_tasks_str ='{}')".format(info_str, comp.nimbusQc_tasks_str))
          #logger.debug("rave_pgf_composite_plugin.py [{}]:generate(): result.getAttribute('how/task') ='{}')".format(info_str, how_task_value0))
          # Note: comp.how_tasks == '' (empty) therefore not usable.
          # - or - when satfilter active
        except:
          pass
      
      try:
        nodes_attr = result.getAttribute('how/nodes')
      except:
        nodes_attr = "-empty/missing-"
      logger.info("rave_pgf_composite_plugin.py [{}]:generate(): getAttribute('how/nodes') = '{}'".format(info_str, nodes_attr))

      how_tasks = "{},{}".format(comp.whichNimbus_qc,comp.nimbusQc_short_tasks_str)
      if how_tasks != "":
        result.addAttribute('how/task', how_tasks)
        logger.info("rave_pgf_composite_plugin.py [{}]:generate(): addAttribute('how/task', '{}')".format(info_str, how_tasks))

      if comp.nimbusQc_satfilter_file !="":
        how_task_args_value = "satfilter_file:{};{}".format(comp.nimbusQc_satfilter_file, comp.nimbusQc_dict_qc_nodes_str)
      else:
        how_task_args_value = "{}".format(comp.nimbusQc_dict_qc_nodes_str)
      if how_task_args_value != "":
          result.addAttribute('how/task_args', how_task_args_value)
          logger.info("rave_pgf_composite_plugin.py [{}]:generate(): addAttribute('how/task_args', '{}')".format(info_str, how_task_args_value))

  _, outfile = rave_tempfile.mktemp(suffix='.h5', close="True")
  
  rio = _raveio.new()
  rio.object = result
  rio.filename = outfile
  rio.version = RAVE_IO_DEFAULT_VERSION
  rio.save()
  logger.info("rave_pgf_composite_plugin.py [{}]:generate(): Generated new composite: {} ".format(info_str, outfile))
  return outfile
  
  
