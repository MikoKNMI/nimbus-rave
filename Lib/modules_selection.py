
## Fonction to select the modules needed in the pre-processing
#  for each radar according to a XML file.

## @file
## @author Iseline Pechin, Meteo-France
## @date 2019-06-24

import rave_defines
import xml.etree.ElementTree as ET
import odim_source


CONFIG_FILE = rave_defines.RAVECONFIG + '/modules_by_radar.xml'

initialized = 0
dict_modules = {}

# Initializes the dictionary with the list of the modules for each radar by 
# reading config from XML file
def init_dictionary():
    global initialized
    if initialized: return
    
    C = ET.parse(CONFIG_FILE)
    OPTIONS = C.getroot()
    
    #For all the radars in the file, get the list of the modules needed
    for site in list(OPTIONS):       
        for k in site.attrib.keys():
            if   k == "modules": 
                dict_modules[site.tag] = site.attrib[k]    
                
    initialized = 1
 
# Reads in the dictionary the modules needed by the radar of the current scan or volume 
def getModules(object): 
    try :       
        init_dictionary()
        
        #Get the nod_id of the radar
        nod_id=odim_source.NODfromSource(object)
        
        #Get the corresponding list of modules
        if nod_id in dict_modules : 
            modules = dict_modules[nod_id]
        else : 
            modules = dict_modules['default']
        return modules
    except : 
        print "The XML file with list of modules by radars is missing. All the modules called by the toolbox will be applied to all the radars"
        return 0
    
#Tests
#init_dictionary()
#print dict_modules['default'].split(',')
#for cle in dict_modules :
#    print cle
#for valeur in ARGS.values() :
#    print valeur
#for site in list(OPTIONS) :
#    print(site)
#    for k in site.attrib.keys():
#        print k
#        print site.attrib[k]       
#if 'xyabc'in ARGS : print 1
#    
#print dict_modules['default'].split(',')
