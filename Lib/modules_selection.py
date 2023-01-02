
## Fonction to select the modules needed in the pre-processing
#  for each radar according to a XML file.

## @file
## @author Iseline Pechin, Meteo-France
## @date 2019-06-24
## @author Michal Koutek, KNMI (Fixes needed for OPERA NIMBUS)
## @date 2022-01-06

import rave_defines
import xml.etree.ElementTree as ET
import odim_source

import rave_pgf_logger
logger = rave_pgf_logger.create_logger()


CONFIG_FILE_PRIMARY_LOCATION = rave_defines.RAVECONFIG + '/modules_by_radar.xml'
CONFIG_FILE_DEFAULT_LOCATION = '/etc/baltrad/rave/config/modules_by_radar.xml'

# NOTE behavior of rave_defines.py:
# --------------------------------
# RAVEROOT = os.path.split(os.path.split(rave_defines.__file__)[0])[0]
# if not RAVEROOT: RAVEROOT = '..'
#
# RAVELIB    = RAVEROOT + '/Lib'
# RAVECONFIG = RAVEROOT + '/config'
# RAVEDB     = RAVEROOT + '/db'
# RAVEBIN    = RAVEROOT + '/bin'
# RAVEETC    = RAVEROOT + '/etc'

initialized = 0
dict_modules = {}

# Initializes the dictionary with the list of the modules for each radar by 
# reading config from XML file
def init_dictionary():
    global initialized
    if initialized: return
    try:
      C = ET.parse(CONFIG_FILE_PRIMARY_LOCATION)
    except:
      C = ET.parse(CONFIG_FILE_DEFAULT_LOCATION)
      # To prevent this exeption:
      # toolbox/nimbus-rave/Lib/modules_selection.py:
      # FileNotFoundError: [Errno 2] No such file or directory: 
      #    /usr/lib/rave/config/modules_by_radar.xml
      # File located here:
      #      /etc/baltrad/rave/config/modules_by_radar.xml
      # SOLUTION:
      # ln -s /etc/baltrad/rave/config/modules_by_radar.xml /usr/lib/rave/config/modules_by_radar.xml

    OPTIONS = C.getroot()
    
    #For all the radars in the file, get the list of the modules needed
    for site in list(OPTIONS):
        for k in site.attrib.keys():
            if   k == "modules": 
                dict_modules[site.tag] = site.attrib[k]
    initialized = 1
 
# Reads in the dictionary the modules needed by the radar of the current scan or volume 
def getModules(object): 
    try:
        init_dictionary()
    except:
        print("getModules(): Failed with init_dictionary()")
        logger.error("getModules(): Failed with init_dictionary()")
        return ""
    try:
        #Get the nod_id of the radar
        nod_id = odim_source.NODfromSource(object)
    except:
        print("getModules(): Failed with  odim_source.NODfromSource(object)")
        logger.error("getModules(): Failed with  odim_source.NODfromSource(object)")
        return ""
    
    #Get the corresponding list of modules
    if nod_id in dict_modules : 
        modules = dict_modules[nod_id]
    else:
        modules = dict_modules['default']
    return modules

    logger.error("The XML file with list of modules by radars is missing. All the modules called by the toolbox will be applied to all the radars")
    return ""
