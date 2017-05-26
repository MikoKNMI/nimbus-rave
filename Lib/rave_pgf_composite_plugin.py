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
import string
import datetime
import rave_tempfile
import logging
import rave_pgf_logger
import rave_dom_db
import rave_util
import _raveio
import rave_tile_registry

from tiled_compositing import tiled_compositing
from compositing import compositing

from rave_defines import CENTER_ID, GAIN, OFFSET, RAVE_PGF_QUALITY_FIELD_REPROCESSING

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
  except ValueError, e:
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
    comp.detectors = string.split(args["anomaly-qc"], ",")

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
  
  if "qitotal_field" in args.keys():
    comp.qitotal_field = args["qitotal_field"]
    
  if "reprocess_qfields" in args.keys():
    comp.reprocess_quality_field = args["reprocess_qfields"]
  else:
    comp.reprocess_quality_field = RAVE_PGF_QUALITY_FIELD_REPROCESSING
  
  # We do not support best fit compositing right now
  #comp.pcsid = options.pcsid
  #comp.xscale = options.scale
  #comp.yscale = options.scale
  
  #if options.gf: Activate gap filling for rule
  #  comp.applygapfilling = True
  
  # Optional cloud-type residual non-precip filter
  if args.has_key("ctfilter"):
    if eval(args["ctfilter"]):
      comp.applyctfilter = True
  
  if args.has_key("applygra"):
    comp.applygra = True
  if args.has_key("zrA"):
    comp.zr_A = float(args["zrA"])
  if args.has_key("zrb"):
    comp.zr_b = float(args["zrb"])
  
  
  if rave_tile_registry.has_tiled_area(args["area"]):
    comp = tiled_compositing(comp)
  
  result = comp.generate(args["date"], args["time"], args["area"])
  
  if result == None:
    logger.info("No composite could be generated.")
    return None
  
  fileno, outfile = rave_tempfile.mktemp(suffix='.h5', close="True")
  
  rio = _raveio.new()
  rio.object = result
  rio.filename = outfile
  rio.save()
  
  return outfile
  
  
