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
## @author Daniel Michelson, SMHI
## @date 2012-11-05

import sys, os, time, traceback, shutil
import logging
import multiprocessing
import _raveio
import rave_pgf_quality_registry
import rave_pgf_logger
import odc_fixIO


opath = '/dev/shm/odc'  # command-line options will override these variables
algorithm_ids = None
delete = False
check = True


# -----------------------------------------------------------------------------
## Run Quality Controls
# @param pload input polar data, either a SCAN or a PVOL object
# @param new version of pload, updated with quality field(s) generated by
# the chosen QC algorithms
def QC(pload):
    for a in algorithm_ids:
        p = rave_pgf_quality_registry.get_plugin(a)
        if not p:
            raise AttributeError, "Could not find %s plugin." % a
        pload, _ = p.process(pload)
        if isinstance(pload, tuple):
          pload, algorithm = pload[0],pload[1]
    return pload


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
    
    
# -----------------------------------------------------------------------------
## Generator, includes gathering of timing information for benchmarking.
# @param ifstr string of the input file
# @return tuple of strings, containing input file and return status,
# either "OK", "EXISTS", or a Traceback message. In the case of "OK",
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
        rio = _raveio.open(ifstr)
        endread = time.time()

        odc_fixIO.Validate(rio)
        endval = time.time()

        pload = rio.object

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
