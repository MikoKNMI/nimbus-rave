/* --------------------------------------------------------------------
Copyright (C) 2011 Swedish Meteorological and Hydrological Institute, SMHI,

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
------------------------------------------------------------------------*/
/**
 * Utilities when working with ODIM H5 files.
 * @file
 * @author Anders Henja (Swedish Meteorological and Hydrological Institute, SMHI)
 * @date 2011-09-30
 */
#include "odim_io_utilities.h"
#include "rave_debug.h"
#include "rave_alloc.h"
#include "rave_hlhdf_utilities.h"
#include <string.h>
#include <math.h>

/*@{ Private functions */

/**
 * Called when an attribute belonging to a rave field
 * is found.
 * @param[in] object - the PolarOdimArg pointing to a rave field
 * @param[in] attribute - the attribute found
 * @return 1 on success otherwise 0
 */
static int OdimIoUtilitiesInternal_loadFieldAttribute(void* object, RaveAttribute_t* attribute)
{
  RaveField_t* field = NULL;
  const char* name;
  int result = 0;
  RaveIO_ODIM_Version version = RaveIO_ODIM_Version_2_4;
  RAVE_ASSERT((object != NULL), "object == NULL");
  RAVE_ASSERT((attribute != NULL), "attribute == NULL");

  field = (RaveField_t*)((OdimIoUtilityArg*)object)->object;
  name = RaveAttribute_getName(attribute);
  version = (RaveIO_ODIM_Version)((OdimIoUtilityArg*)object)->version;
  if (name != NULL) {
    OdimIoUtilities_convertHowAttributeToInternalRave(name, version, attribute);
    result = RaveField_addAttribute(field, attribute);
  }

  return result;
}

/**
 * Called when an dataset belonging to a field parameter
 * is found.
 * @param[in] object - the PolarOdimArg pointing to a rave field
 * @param[in] xsize - the x size
 * @param[in] ysize - the y size
 * @param[in] data  - the data
 * @param[in] dtype - the type of the data.
 * @return 1 on success otherwise 0
 */
static int OdimIoUtilitiesInternal_loadFieldDataset(void* object, hsize_t xsize, hsize_t ysize, void* data, RaveDataType dtype, const char* nodeName)
{
  RaveField_t* field = NULL;
  int result = 0;

  field = (RaveField_t*)((OdimIoUtilityArg*)object)->object;
  if (data == NULL && ((OdimIoUtilityArg*)object)->lazyReader != NULL) {
    LazyDataset_t* datasetReader = RAVE_OBJECT_NEW(&LazyDataset_TYPE);
    if (datasetReader != NULL) {
      result = LazyDataset_init(datasetReader, ((OdimIoUtilityArg*)object)->lazyReader, nodeName);
    }
    if (result) {
      result = RaveField_setLazyDataset(field, datasetReader);
    }
    RAVE_OBJECT_RELEASE(datasetReader);
  } else {
    result = RaveField_setData(field, xsize, ysize, data, dtype);
  }
  return result;
}

/*@} End of Private functions */

/*@{ Interface functions */
int OdimIoUtilities_convertHowAttributeToInternalRave(const char* attrname, RaveIO_ODIM_Version origversion, RaveAttribute_t* inattr)
{
  if (origversion < RaveIO_ODIM_Version_2_4) {
    return 1;
  }

  if (RaveAttribute_getFormat(inattr) == RaveAttribute_Format_Double) {
    double v = 0.0;
    RaveAttribute_getDouble(inattr, &v);
    if (strcasecmp("how/gasattn", attrname) == 0) {
      RaveAttribute_setDouble(inattr, v * 1000.0); /* dB/m => dB/km */
    } else if (strcasecmp("how/minrange", attrname)==0 ||
        strcasecmp("how/maxrange", attrname)==0 ||
        strcasecmp("how/melting_layer_top_A", attrname)==0 ||
        strcasecmp("how/melting_layer_bottom_A", attrname)==0) {
      RaveAttribute_setDouble(inattr, v / 1000.0); /* m => km */
    } else if (strcasecmp("how/nomTXpower", attrname) == 0 ||
        strcasecmp("how/peakpwr", attrname) == 0 ||
        strcasecmp("how/avgpwr", attrname) == 0) {
      RaveAttribute_setDouble(inattr, pow(10.0, (v - 30.0)/10.0)/1000.0);
    }
  } else if (RaveAttribute_getFormat(inattr) == RaveAttribute_Format_DoubleArray && strcasecmp("how/TXpower", attrname) == 0) {
    double* darr = NULL;
    int nlen = 0, i = 0;
    RaveAttribute_getDoubleArray(inattr, &darr, &nlen);
    for (i = 0; i < nlen; i++) {
      darr[i] = pow(10.0, (darr[i] - 30.0)/10.0)/1000.0;
    }
  }

  return 1;
}

int OdimIoUtilities_convertHowAttributeFromInternalRave(const char* attrname, RaveIO_ODIM_Version outversion, RaveAttribute_t* inattr)
{
  if (outversion < RaveIO_ODIM_Version_2_4) {
    return 1;
  }
  /* Since all attribute units within rave are according to old . We have to change unit to old one */
  if (RaveAttribute_getFormat(inattr) == RaveAttribute_Format_Double) {
    double v = 0.0;
    RaveAttribute_getDouble(inattr, &v);
    if (strcasecmp("how/gasattn", attrname) == 0) {
      RaveAttribute_setDouble(inattr, v / 1000.0); /* dB/km => dB/m */
    } else if (strcasecmp("how/minrange", attrname)==0 ||
        strcasecmp("how/maxrange", attrname)==0 ||
        strcasecmp("how/melting_layer_top_A", attrname)==0 ||
        strcasecmp("how/melting_layer_bottom_A", attrname)==0) {
      RaveAttribute_setDouble(inattr, v * 1000.0); /* km => m */
    } else if (strcasecmp("how/nomTXpower", attrname) == 0 ||
        strcasecmp("how/peakpwr", attrname) == 0 ||
        strcasecmp("how/avgpwr", attrname) == 0) {
      if (v > 0) {
        RaveAttribute_setDouble(inattr, 10*log10(1000.0*v)+30); /* kW => dBm */
      }
    }
  } else if (RaveAttribute_getFormat(inattr) == RaveAttribute_Format_DoubleArray && strcasecmp("how/TXpower", attrname) == 0) {
    double* darr = NULL;
    int nlen = 0, i = 0;
    RaveAttribute_getDoubleArray(inattr, &darr, &nlen);
    for (i = 0; i < nlen; i++) {
      if (darr[i] > 0) {
        darr[i] = 10 * log10(1000.0*darr[i]) + 30; /* kW => dBm */
      }
    }
  }

  return 1;
}

int OdimIoUtilities_convertHowAttributesFromInternalRave(RaveObjectList_t* attributes, RaveIO_ODIM_Version outversion)
{
  int i = 0;
  int nlen = 0;
  if (outversion < RaveIO_ODIM_Version_2_4) {
    return 1;
  }
  nlen = RaveObjectList_size(attributes);
  for (i = 0; i < nlen; i++) {
    RaveAttribute_t* attr = (RaveAttribute_t*)RaveObjectList_get(attributes, i);
    if (attr != NULL) {
      OdimIoUtilities_convertHowAttributeFromInternalRave(RaveAttribute_getName(attr), outversion, attr);
    }
    RAVE_OBJECT_RELEASE(attr);
  }
  return 1;
}

int OdimIoUtilities_addRaveField(RaveField_t* field, HL_NodeList* nodelist, RaveIO_ODIM_Version version, const char* fmt, ...)
{
  int result = 0;
  va_list ap;
  char name[1024];
  int nName = 0;
  RaveObjectList_t* attributes = NULL;

  RAVE_ASSERT((field != NULL), "field == NULL");
  RAVE_ASSERT((nodelist != NULL), "nodelist == NULL");
  RAVE_ASSERT((fmt != NULL), "fmt == NULL");

  va_start(ap, fmt);
  nName = vsnprintf(name, 1024, fmt, ap);
  va_end(ap);
  if (nName < 0 || nName >= 1024) {
    RAVE_ERROR0("NodeName would evaluate to more than 1024 characters.");
    goto done;
  }
  if (!RaveHL_hasNodeByName(nodelist, name)) {
    if (!RaveHL_createGroup(nodelist, name)) {
      goto done;
    }
  }

  attributes = RaveField_getAttributeValues(field);

  OdimIoUtilities_convertHowAttributesFromInternalRave(attributes, version);

  if (attributes == NULL || !RaveHL_addAttributes(nodelist, attributes, name)) {
    goto done;
  }

  if (!RaveHL_addData(nodelist,
                      RaveField_getData(field),
                      RaveField_getXsize(field),
                      RaveField_getYsize(field),
                      RaveField_getDataType(field),
                      name)) {
    goto done;
  }

  result = 1;
done:
  RAVE_OBJECT_RELEASE(attributes);
  return result;
}

/**
 * Adds a list of quality fields (RaveField_t) to a nodelist.
 *
 * @param[in] fields - the list of fields
 * @param[in] nodelist - the hlhdf node list
 * @param[in] fmt - the varargs format string
 * @param[in] ... - the varargs
 * @return 1 on success otherwise 0
 */
int OdimIoUtilities_addQualityFields(RaveObjectList_t* fields, HL_NodeList* nodelist, RaveIO_ODIM_Version version, const char* fmt, ...)
{
  int result = 0;
  va_list ap;
  char name[1024];
  int nName = 0;
  int pindex = 0;
  int nrfields = 0;

  RAVE_ASSERT((fields != NULL), "fields == NULL");
  RAVE_ASSERT((nodelist != NULL), "nodelist == NULL");
  RAVE_ASSERT((fmt != NULL), "fmt == NULL");

  va_start(ap, fmt);
  nName = vsnprintf(name, 1024, fmt, ap);
  va_end(ap);
  if (nName < 0 || nName >= 1024) {
    RAVE_ERROR0("NodeName would evaluate to more than 1024 characters.");
    goto done;
  }

  result = 1;
  nrfields = RaveObjectList_size(fields);
  for (pindex = 0; result == 1 && pindex < nrfields; pindex++) {
    RaveField_t* field = (RaveField_t*)RaveObjectList_get(fields, pindex);
    if (field != NULL) {
      result = OdimIoUtilities_addRaveField(field, nodelist, version, "%s/quality%d", name, (pindex+1));
    } else {
      result = 0;
    }
    RAVE_OBJECT_RELEASE(field);
  }

done:
  return result;
}

/**
 * Loads a rave field. A rave field can be just about anything with a mapping
 * between attributes and a dataset.
 * @param[in] nodelist - the hlhdf node list
 * @param[in] nodelist - version of the file read
 * @param[in] fmt - the variable argument list string format
 * @param[in] ... - the variable argument list
 * @return a rave field on success otherwise NULL
 */
RaveField_t* OdimIoUtilities_loadField(LazyNodeListReader_t* lazyReader, RaveIO_ODIM_Version version, const char* fmt, ...)
{
  OdimIoUtilityArg arg;
  RaveField_t* field = NULL;
  RaveField_t* result = NULL;
  va_list ap;
  char name[1024];
  int nName = 0;

  RAVE_ASSERT((lazyReader != NULL), "lazyReader == NULL");
  RAVE_ASSERT((fmt != NULL), "fmt == NULL");

  va_start(ap, fmt);
  nName = vsnprintf(name, 1024, fmt, ap);
  va_end(ap);
  if (nName < 0 || nName >= 1024) {
    RAVE_ERROR0("NodeName would evaluate to more than 1024 characters.");
    goto fail;
  }

  field = RAVE_OBJECT_NEW(&RaveField_TYPE);
  if (field == NULL) {
    RAVE_CRITICAL0("Failed to allocate memory for field");
    goto fail;
  }
  arg.lazyReader = lazyReader;
  arg.nodelist = LazyNodeListReader_getHLNodeList(lazyReader);
  arg.object = (RaveCoreObject*)field;
  arg.version = version;

  if (!RaveHL_loadAttributesAndData(arg.nodelist, &arg,
                                    OdimIoUtilitiesInternal_loadFieldAttribute,
                                    OdimIoUtilitiesInternal_loadFieldDataset,
                                    name)) {
    goto fail;
  }

  result = RAVE_OBJECT_COPY(field);
fail:
  RAVE_OBJECT_RELEASE(field);
  return result;
}


int OdimIoUtilities_getIdFromSource(const char* source, const char* id, char* buf, size_t buflen)
{
  int result = 0;
  if (source != NULL && id != NULL) {
    char* p = strstr(source, id);
    if (p != NULL) {
      int len = 0;
      char* pbrk = NULL;
      p += strlen(id);
      len = strlen(p);
      pbrk = strpbrk((const char*)p, ",");

      if (pbrk != NULL) {
        len = pbrk - p;
      }
      if (len + 1 < buflen) {
        strncpy(buf, p, len);
        buf[len] = '\0';
        result = 1;
      }
    }
  }
  return result;
}

int OdimIoUtilities_getNodOrCmtFromSource(const char* source, char* buf, size_t buflen)
{
  int result = 0;
  result = OdimIoUtilities_getIdFromSource(source, "NOD:", buf, buflen);
  if (!result)
    result = OdimIoUtilities_getIdFromSource(source, "CMT:", buf, buflen);
  return result;
}
