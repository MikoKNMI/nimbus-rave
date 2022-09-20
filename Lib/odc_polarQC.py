#!/usr/bin/env python
'''
Copyright (C) 2012- Swedish Meteorological and Hydrological Institute (SMHI)

This file is part of RAVE.

RAVE is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RAVE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with RAVE.  If not, see <http://www.gnu.org/licenses/>.
'''

## Quality controls polar data for Odyssey

## @file
## @author Daniel Michelson, SMHI (Original Code)
## @date 2012-11-05
## @author Maud Martet, Meteo France (Add delete_GTS_header)
## @date 2014-06-17
## @author Iseline Pechin, Meteo France (Add modules_selection)
## @date 2019-06-24

import sys, os, time, traceback, shutil
import re
import logging
import multiprocessing
import datetime
import _raveio
import _polarvolume
import _polarscan
import rave_pgf_quality_registry
import rave_pgf_logger
import odc_fixIO
import modules_selection


opath = '/dev/shm/odc'  # command-line options will override these variables
algorithm_ids = None
delete = False
check = True
skip_age_seconds = 0


# -----------------------------------------------------------------------------
## Run Quality Controls
# @param pload input polar data, either a SCAN or a PVOL object
# @param new version of pload, updated with quality field(s) generated by
# the chosen QC algorithms
def QC(pload):
    # This part allows to choose the modules independently for each radar
    modules = modules_selection.getModules(pload)
    if modules != 0 : 
    	algorithm_ids=modules.split(',')
    ########
    for a in algorithm_ids:
        p = rave_pgf_quality_registry.get_plugin(a)
        if not p:
            raise AttributeError("Could not find %s plugin." % a)
        pload, _ = p.process(pload)
        if isinstance(pload, tuple):
          pload, algorithm = pload[0],pload[1]
    return pload

## Check whether we should skip (defer) processing of this scan based on the scan endtime VT 
#
#  Skip predicate is NOW() - SCAN_VT < skip_age_seconds
#
#  Only reflectivity scans (DBZH, TH) are considered.   
#
#  @return True if the scan processing be skipped for now
def should_skip_scan(scan, filename, logger):
    if not scan.hasParameter("DBZH") and not scan.hasParameter("TH"):
        return False
    
    try:
        date_str = scan.enddate + scan.endtime
    except AttributeError as e:
        rave_pgf_logger.log(logger, "warning", "%s: cannot get scan endtime: %s" % (filename, e))
        return False
        
    # correct somewhat common error where seconds==60, otherwise strptime complains
    if date_str[12:14] == '60':
        date_str = date_str[0:12] + '59'   
    try:
        date = datetime.datetime.strptime(date_str,'%Y%m%d%H%M%S')
    except Exception as e:
        rave_pgf_logger.log(logger, "warning", "%s: cannot parse filename: %s" % (filename, e))
        return False

    age = (datetime.datetime.utcnow() - date).total_seconds()
    skip = age < skip_age_seconds
    rave_pgf_logger.log(logger, "debug", "%s: age: %f skip: %s" % (filename, age, skip))
    return skip

## Check whether we should skip (defer) processing of the scan or volume 
#
#  @return True if the scan processing be skipped for now
def should_skip(rio_object, filename, logger):
    if skip_age_seconds<=0:
        return False
     
    if _polarvolume.isPolarVolume(rio_object):
        for s in range(rio_object.getNumberOfScans()):
            if should_skip_scan(rio_object.getScan(s), filename, logger):
                return True
        
    elif _polarscan.isPolarScan(rio_object):
        return should_skip_scan(rio_object, filename, logger)
        
    else:
        rave_pgf_logger.log(logger, "warning", "%s: neither PVOL nor SCAN" % (filename, e))

    return False


## Create an output file string and check for the presence of an output file.
# @param ifstr string input file string
# @return tuple containing output file string and a True/False on its presence.
# Will return False if the input file has been copied to the output directory.
def MakeCheckOfstr(ifstr):
    path, fstr = os.path.split(ifstr)    
    ofstr = fstr.split('.')[0] + '.h5'
    ofstr = os.path.join(opath, ofstr)
    newfile =  os.path.isfile(ofstr) and os.path.getsize(ofstr)
    # Check also if there's a copy of the input file in the output directory
    oldfile = os.path.join(opath, fstr)
    copied = os.path.isfile(oldfile) and os.path.getsize(oldfile)
    if newfile or copied:
        return ofstr, True
    else:
        return ofstr, False
    

## Delete GTS from files
# @param ifstr string input file string
#
# This could handle (LPMG* files from Portugal that contain head like this:
# PAZZ43 LPMG 271010
# <89>HDF  .. this is normal begin of HDF file
def delete_GTS_header(ifstr):
    fic=open(ifstr,"rb")
    # Read characters 1 to 3 of file in 'type' string
    type=""
    for i in range(1,4) :
        fic.seek(i)
        s=fic.read(1)
        type="%s%s" % (type,s)

    # If 'type' = UFR or HDF, there is no GTS in file
    # otherwise delete first line of the file
    if (type != 'UFR' and type != "HDF"):
        fic.seek(0)
        fic.readline()
        content = fic.read() 
        fic.close()
        fic=open(ifstr,"wb")
        fic.write(content)
    else :
        fic.close()
    
    
# -----------------------------------------------------------------------------
## Generator, includes gathering of timing information for benchmarking.
# @param ifstr string of the input file
# @return tuple of strings, containing input file and return status,
# either "OK", "EXISTS", "SKIPPED" or a Traceback message. In the case of "OK",
# benchmarking statistics are also included in a tuple
def generate(ifstr):
    logger = logging.getLogger("ODC")
    rave_pgf_logger.init_logger(logger)

    ofstr, done = MakeCheckOfstr(ifstr)
    path, fstr = os.path.split(ifstr)

    if check and done:
        if delete:
            os.remove(ifstr)
        rave_pgf_logger.log(logger, "debug", "%s: EXISTS" % fstr)
        return ifstr, "EXISTS"

    
    try:
        startread = time.time()
        delete_GTS_header(ifstr)
        rio = _raveio.open(ifstr)
        endread = time.time()

        odc_fixIO.Validate(rio)
        endval = time.time()

        pload = rio.object

        if should_skip(pload, fstr, logger):
            rave_pgf_logger.log(logger, "debug", "%s: SKIPPED" % fstr)
            return ifstr, "SKIPPED"

        pload = QC(pload)
        endqc = time.time()

        # Hard-wire for no compression and optimized file-creation properties.
        rio.compression_level = 0
        rio.fcp_istorek = 1
        rio.fcp_metablocksize = 0
        rio.fcp_sizes = (4,4)
        rio.fcp_symk = (1,1)
        rio.fcp_userblock = 0

        rio.object = pload
        rio.save(ofstr)
        endwrite = time.time()

        readt = endread - startread
        validt = endval- endread
        qct = endqc - endval
        writet = endwrite - endqc

        if delete:
            os.remove(ifstr)  # Clean-up input directory one file at a time.
        rave_pgf_logger.log(logger, "info", "%s: OK" % fstr)

        return ifstr, "OK", (readt, validt, qct, writet)
    except Exception:
        err_msg = traceback.format_exc()
        if delete:
            os.rename(ifstr, os.path.join(opath, fstr))
        else:
            shutil.copyfile(ifstr, os.path.join(opath, fstr))
        rave_pgf_logger.log(logger, "error", "%s: %s" % (fstr, err_msg))
        return ifstr, err_msg
  

## Distributes 'generate' among the available CPU cores on this machine.
#  @param fstrs list of input file name strings
#  @param procs int number of concurrent processes, defaults to the max allowed
#  @return list of returned tuples from \ref generate
def multi_generate(fstrs, procs=None):
    pool = multiprocessing.Pool(procs)

    results = []
    r = pool.map_async(generate, fstrs, chunksize=1, callback=results.append)
    r.wait()

    return results[0]


if __name__ == "__main__":
    pass
