/* --------------------------------------------------------------------
Copyright (C) 2009 Swedish Meteorological and Hydrological Institute, SMHI,

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
 * Python version of the VerticalProfile API.
 * @file
 * @author Anders Henja (Swedish Meteorological and Hydrological Institute, SMHI)
 * @date 2012-08-24
 *
 * @co-author Ulf Nordh (Swedish Meteorological and Hydrological Institute, SMHI)
 * @date 2017-10-27
 * Uppdated code with new vertical profile field definitions (HGHT, n, UWND, VWND)
 * and functions for accessing starttime, endtime, startdate, endate and product
 */
#include "pyravecompat.h"
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#define PYVERTICALPROFILE_MODULE   /**< to get correct part of pyverticalprofile,h */
#include "pyverticalprofile.h"

#include <arrayobject.h>
#include "pyrave_debug.h"
#include "pyravefield.h"
#include "rave_alloc.h"
#include "raveutil.h"
#include "rave.h"

/**
 * Debug this module
 */
PYRAVE_DEBUG_MODULE("_verticalprofile");

/**
 * Sets a python exception and goto tag
 */
#define raiseException_gotoTag(tag, type, msg) \
{PyErr_SetString(type, msg); goto tag;}

/**
 * Sets a python exception and return NULL
 */
#define raiseException_returnNULL(type, msg) \
{PyErr_SetString(type, msg); return NULL;}

/**
 * Error object for reporting errors to the python interpreeter
 */
static PyObject *ErrorObject;

/// --------------------------------------------------------------------
/// Vertical Profile
/// --------------------------------------------------------------------
/*@{ Vertical Profile */
/**
 * Returns the native VerticalProfile_t instance.
 * @param[in] pyobj - the python vertical profile instance
 * @returns the native vertical profile instance.
 */
static VerticalProfile_t*
PyVerticalProfile_GetNative(PyVerticalProfile* pyobj)
{
  RAVE_ASSERT((pyobj != NULL), "pyobj == NULL");
  return RAVE_OBJECT_COPY(pyobj->vp);
}

/**
 * Creates a python vertical profile from a native vertical profile or will create an
 * initial native VerticalProfile if p is NULL.
 * @param[in] p - the native vertical profile (or NULL)
 * @returns the python vertical profile
 */
static PyVerticalProfile* PyVerticalProfile_New(VerticalProfile_t* p)
{
  PyVerticalProfile* result = NULL;
  VerticalProfile_t* cp = NULL;

  if (p == NULL) {
    cp = RAVE_OBJECT_NEW(&VerticalProfile_TYPE);
    if (cp == NULL) {
      RAVE_CRITICAL0("Failed to allocate memory for vertical profile.");
      raiseException_returnNULL(PyExc_MemoryError, "Failed to allocate memory for vertical profile.");
    }
  } else {
    cp = RAVE_OBJECT_COPY(p);
    result = RAVE_OBJECT_GETBINDING(p); // If p already have a binding, then this should only be increfed.
    if (result != NULL) {
      Py_INCREF(result);
    }
  }

  if (result == NULL) {
    result = PyObject_NEW(PyVerticalProfile, &PyVerticalProfile_Type);
    if (result != NULL) {
      PYRAVE_DEBUG_OBJECT_CREATED;
      result->vp = RAVE_OBJECT_COPY(cp);
      RAVE_OBJECT_BIND(result->vp, result);
    } else {
      RAVE_CRITICAL0("Failed to create PyVerticalProfile instance");
      raiseException_gotoTag(done, PyExc_MemoryError, "Failed to allocate memory for vertical profile.");
    }
  }
done:
  RAVE_OBJECT_RELEASE(cp);
  return result;
}

/**
 * Deallocates the vertical profile
 * @param[in] obj the object to deallocate.
 */
static void _pyverticalprofile_dealloc(PyVerticalProfile* obj)
{
  if (obj == NULL) {
    return;
  }
  PYRAVE_DEBUG_OBJECT_DESTROYED;
  RAVE_OBJECT_UNBIND(obj->vp, obj);
  RAVE_OBJECT_RELEASE(obj->vp);
  PyObject_Del(obj);
}

/**
 * Creates a new instance of the vertical profile.
 * @param[in] self this instance.
 * @param[in] args arguments for creation (NOT USED).
 * @return the object on success, otherwise NULL
 */
static PyObject* _pyverticalprofile_new(PyObject* self, PyObject* args)
{
  PyVerticalProfile* result = PyVerticalProfile_New(NULL);
  return (PyObject*)result;
}

static PyObject* _pyverticalprofile_setLevels(PyVerticalProfile* self, PyObject* args)
{
  long l = 0;
  if (!PyArg_ParseTuple(args, "l", &l)) {
    return NULL;
  }
  if (!VerticalProfile_setLevels(self->vp, l)){
    raiseException_returnNULL(PyExc_AttributeError, "Failed to set level count");
  }
  Py_RETURN_NONE;
}

static PyObject* _pyverticalprofile_getLevels(PyVerticalProfile* self, PyObject* args)
{
  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }
  return PyLong_FromLong(VerticalProfile_getLevels(self->vp));
}


/**
 * Special version of the setFF/setFFDev/.. that takes a third argument for passing on quantities
 * and using that as switch for making correct call.
 * @param[in] self - this instance
 * @param[in] args - object, of type RaveFieldCore
 * @param[in] quantity - the quantity
 * @returns None
 */
static PyObject* _pyverticalprofile_internal_setfield(PyVerticalProfile* self, PyObject* args, const char* quantity)
{
  PyObject* inptr = NULL;
  PyRaveField* ravefield = NULL;

  if (!PyArg_ParseTuple(args, "O", &inptr)) {
    return NULL;
  }

  if (!PyRaveField_Check(inptr)) {
    raiseException_returnNULL(PyExc_TypeError,"Object must be of type RaveFieldCore");
  }

  ravefield = (PyRaveField*)inptr;

  if (strcmp("ff", quantity) == 0) {
    if (!VerticalProfile_setFF(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set ff");
    }
  } else if (strcmp("ff_dev", quantity) == 0) {
    if (!VerticalProfile_setFFDev(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set ff_dev");
    }
  } else if (strcmp("w", quantity) == 0) {
    if (!VerticalProfile_setW(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set w");
    }
  } else if (strcmp("w_dev", quantity) == 0) {
    if (!VerticalProfile_setWDev(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set w_dev");
    }
  } else if (strcmp("dd", quantity) == 0) {
    if (!VerticalProfile_setDD(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set dd");
    }
  } else if (strcmp("dd_dev", quantity) == 0) {
    if (!VerticalProfile_setDDDev(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set dd_dev");
    }
  } else if (strcmp("div", quantity) == 0) {
    if (!VerticalProfile_setDiv(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set div");
    }
  } else if (strcmp("div_dev", quantity) == 0) {
    if (!VerticalProfile_setDivDev(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set div_dev");
    }
  } else if (strcmp("def", quantity) == 0) {
    if (!VerticalProfile_setDef(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set def");
    }
  } else if (strcmp("def_dev", quantity) == 0) {
    if (!VerticalProfile_setDefDev(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set def_dev");
    }
  } else if (strcmp("ad", quantity) == 0) {
    if (!VerticalProfile_setAD(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set ad");
    }
  } else if (strcmp("ad_dev", quantity) == 0) {
    if (!VerticalProfile_setADDev(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set ad_dev");
    }
  } else if (strcmp("dbz", quantity) == 0) {
    if (!VerticalProfile_setDBZ(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set dbz");
    }
  } else if (strcmp("dbz_dev", quantity) == 0) {
    if (!VerticalProfile_setDBZDev(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set dbz_dev");
    }
  } else if (strcmp("HGHT", quantity) == 0) {
    if (!VerticalProfile_setHGHT(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set HGHT");
    }
  } else if (strcmp("n", quantity) == 0) {
    if (!VerticalProfile_setNV(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set n");
    }
  } else if (strcmp("UWND", quantity) == 0) {
    if (!VerticalProfile_setUWND(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set UWND");
    }
  } else if (strcmp("VWND", quantity) == 0) {
    if (!VerticalProfile_setVWND(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to set VWND");
    }
  } else {
    raiseException_returnNULL(PyExc_AssertionError, "Programming error");
  }

  Py_RETURN_NONE;
}

/**
 * Special version of the getFF/getFFDev/.. that takes a third argument for passing on quantities
 * and using that as switch for making correct call.
 * @param[in] self - this instance
 * @param[in] args - None
 * @param[in] quantity - the quantity
 * @returns The field
 */
static PyObject* _pyverticalprofile_internal_getfield(PyVerticalProfile* self, PyObject* args, const char* quantity)
{
  PyObject* result = NULL;
  RaveField_t* field = NULL;
  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }

  if (strcmp("ff", quantity) == 0) {
    field = VerticalProfile_getFF(self->vp);
  } else if (strcmp("ff_dev", quantity) == 0) {
    field = VerticalProfile_getFFDev(self->vp);
  } else if (strcmp("w", quantity) == 0) {
    field = VerticalProfile_getW(self->vp);
  } else if (strcmp("w_dev", quantity) == 0) {
    field = VerticalProfile_getWDev(self->vp);
  } else if (strcmp("dd", quantity) == 0) {
    field = VerticalProfile_getDD(self->vp);
  } else if (strcmp("dd_dev", quantity) == 0) {
    field = VerticalProfile_getDDDev(self->vp);
  } else if (strcmp("div", quantity) == 0) {
    field = VerticalProfile_getDiv(self->vp);
  } else if (strcmp("div_dev", quantity) == 0) {
    field = VerticalProfile_getDivDev(self->vp);
  } else if (strcmp("def", quantity) == 0) {
    field = VerticalProfile_getDef(self->vp);
  } else if (strcmp("def_dev", quantity) == 0) {
    field = VerticalProfile_getDefDev(self->vp);
  } else if (strcmp("ad", quantity) == 0) {
    field = VerticalProfile_getAD(self->vp);
  } else if (strcmp("ad_dev", quantity) == 0) {
    field = VerticalProfile_getADDev(self->vp);
  } else if (strcmp("dbz", quantity) == 0) {
    field = VerticalProfile_getDBZ(self->vp);
  } else if (strcmp("dbz_dev", quantity) == 0) {
    field = VerticalProfile_getDBZDev(self->vp);
  } else if (strcmp("HGHT", quantity) == 0) {
    field = VerticalProfile_getHGHT(self->vp);
  } else if (strcmp("n", quantity) == 0) {
    field = VerticalProfile_getNV(self->vp);
  } else if (strcmp("UWND", quantity) == 0) {
    field = VerticalProfile_getUWND(self->vp);
  } else if (strcmp("VWND", quantity) == 0) {
    field = VerticalProfile_getVWND(self->vp);
  } else {
    raiseException_gotoTag(done, PyExc_AssertionError, "Programming error");
  }
  if (field == NULL) {
    Py_RETURN_NONE;
  }
  result = (PyObject*)PyRaveField_New(field);
done:
  RAVE_OBJECT_RELEASE(field);
  return result;
}

static PyObject* _pyverticalprofile_setFF(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "ff");
}

static PyObject* _pyverticalprofile_getFF(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "ff");
}

static PyObject* _pyverticalprofile_setFFDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "ff_dev");
}

static PyObject* _pyverticalprofile_getFFDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "ff_dev");
}

static PyObject* _pyverticalprofile_setW(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "w");
}

static PyObject* _pyverticalprofile_getW(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "w");
}

static PyObject* _pyverticalprofile_setWDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "w_dev");
}

static PyObject* _pyverticalprofile_getWDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "w_dev");
}

static PyObject* _pyverticalprofile_setDD(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "dd");
}

static PyObject* _pyverticalprofile_getDD(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "dd");
}

static PyObject* _pyverticalprofile_setDDDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "dd_dev");
}

static PyObject* _pyverticalprofile_getDDDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "dd_dev");
}

static PyObject* _pyverticalprofile_setDiv(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "div");
}

static PyObject* _pyverticalprofile_getDiv(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "div");
}

static PyObject* _pyverticalprofile_setDivDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "div_dev");
}

static PyObject* _pyverticalprofile_getDivDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "div_dev");
}

static PyObject* _pyverticalprofile_setDef(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "def");
}

static PyObject* _pyverticalprofile_getDef(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "def");
}

static PyObject* _pyverticalprofile_setDefDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "def_dev");
}

static PyObject* _pyverticalprofile_getDefDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "def_dev");
}

static PyObject* _pyverticalprofile_setAD(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "ad");
}

static PyObject* _pyverticalprofile_getAD(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "ad");
}

static PyObject* _pyverticalprofile_setADDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "ad_dev");
}

static PyObject* _pyverticalprofile_getADDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "ad_dev");
}

static PyObject* _pyverticalprofile_setDBZ(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "dbz");
}

static PyObject* _pyverticalprofile_getDBZ(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "dbz");
}

static PyObject* _pyverticalprofile_setDBZDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "dbz_dev");
}

static PyObject* _pyverticalprofile_getDBZDev(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "dbz_dev");
}

static PyObject* _pyverticalprofile_setHGHT(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "HGHT");
}

static PyObject* _pyverticalprofile_getHGHT(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "HGHT");
}

static PyObject* _pyverticalprofile_setNV(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "n");
}

static PyObject* _pyverticalprofile_getNV(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "n");
}

static PyObject* _pyverticalprofile_setUWND(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "UWND");
}

static PyObject* _pyverticalprofile_getUWND(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "UWND");
}

static PyObject* _pyverticalprofile_setVWND(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_setfield(self, args, "VWND");
}

static PyObject* _pyverticalprofile_getVWND(PyVerticalProfile* self, PyObject* args)
{
  return _pyverticalprofile_internal_getfield(self, args, "VWND");
}

static PyObject* _pyverticalprofile_getFields(PyVerticalProfile* self, PyObject* args)
{
  RaveObjectList_t* fields = NULL;
  PyObject* result = NULL;

  int sz, i;

  if (!PyArg_ParseTuple(args, "")) {
    return NULL;
  }
  fields = VerticalProfile_getFields(self->vp);
  if (fields == NULL) {
    raiseException_returnNULL(PyExc_MemoryError, "Failed to get fields");
  }
  sz = RaveObjectList_size(fields);
  result = PyList_New(0);
  if (result == NULL) {
    raiseException_gotoTag(fail, PyExc_MemoryError, "Failed to create list");
  }
  for (i = 0; i < sz; i++) {
    RaveField_t* field = (RaveField_t*)RaveObjectList_get(fields, i);
    if (field != NULL) {
      PyObject* pyfield = (PyObject*)PyRaveField_New(field);
      if (pyfield == NULL || PyList_Append(result, pyfield) < 0) {
        RAVE_OBJECT_RELEASE(field);
        Py_XDECREF(pyfield);
        raiseException_gotoTag(fail, PyExc_MemoryError, "Failed to add item to list");
      }
      Py_XDECREF(pyfield);
    }
    RAVE_OBJECT_RELEASE(field);
  }

  RAVE_OBJECT_RELEASE(fields);
  return result;
fail:
  Py_XDECREF(result);
  RAVE_OBJECT_RELEASE(fields);
  return NULL;
}

static PyObject* _pyverticalprofile_addField(PyVerticalProfile* self, PyObject* args)
{
  PyObject* inptr = NULL;
  PyRaveField* ravefield = NULL;
  RaveAttribute_t* attr = NULL;

  if (!PyArg_ParseTuple(args, "O", &inptr)) {
    return NULL;
  }

  if (!PyRaveField_Check(inptr)) {
    raiseException_returnNULL(PyExc_TypeError,"Object must be of type RaveFieldCore");
  }

  ravefield = (PyRaveField*)inptr;
  attr = RaveField_getAttribute(ravefield->field, "what/quantity");
  if (attr == NULL) {
    raiseException_returnNULL(PyExc_AttributeError, "Field must contain quantity");
  }
  RAVE_OBJECT_RELEASE(attr);

  if (!VerticalProfile_addField(self->vp, ravefield->field)) {
      raiseException_returnNULL(PyExc_AttributeError, "Failed to add field");
  }

  Py_RETURN_NONE;
}

static PyObject* _pyverticalprofile_getField(PyVerticalProfile* self, PyObject* args)
{
  char* quantity = NULL;
  RaveField_t* field = NULL;
  PyObject* pyfield = NULL;

  if (!PyArg_ParseTuple(args, "s", &quantity)) {
    return NULL;
  }

  field = VerticalProfile_getField(self->vp, quantity);
  if (field != NULL) {
    pyfield = (PyObject*)PyRaveField_New(field);
  }
  RAVE_OBJECT_RELEASE(field);
  if (pyfield == NULL) {
    Py_RETURN_NONE;
  }
  return pyfield;
}
/**
 * Adds an attribute to the parameter. Name of the attribute should be in format
 * ^(how|what|where)/[A-Za-z0-9_.]$. E.g how/something, what/sthis etc.
 * Currently, the only supported values are double, long, string.
 * @param[in] self - this instance
 * @param[in] args - bin index, ray index.
 * @returns true or false depending if it works.
 */
static PyObject* _pyverticalprofile_addAttribute(PyVerticalProfile* self, PyObject* args)
{
  RaveAttribute_t* attr = NULL;
  char* name = NULL;
  PyObject* obj = NULL;
  PyObject* result = NULL;

  if (!PyArg_ParseTuple(args, "sO", &name, &obj)) {
    return NULL;
  }

  attr = RAVE_OBJECT_NEW(&RaveAttribute_TYPE);
  if (attr == NULL) {
    return NULL;
  }

  if (!RaveAttribute_setName(attr, name)) {
    raiseException_gotoTag(done, PyExc_MemoryError, "Failed to set name");
  }

  if (PyLong_Check(obj) || PyInt_Check(obj)) {
    long value = PyLong_AsLong(obj);
    RaveAttribute_setLong(attr, value);
  } else if (PyFloat_Check(obj)) {
    double value = PyFloat_AsDouble(obj);
    RaveAttribute_setDouble(attr, value);
  } else if (PyString_Check(obj)) {
    char* value = PyString_AsString(obj);
    if (!RaveAttribute_setString(attr, value)) {
      raiseException_gotoTag(done, PyExc_AttributeError, "Failed to set string value");
    }
  } else if (PyArray_Check(obj)) {
    PyArrayObject* arraydata = (PyArrayObject*)obj;
    if (PyArray_NDIM(arraydata) != 1) {
      raiseException_gotoTag(done, PyExc_AttributeError, "Only allowed attribute arrays are 1-dimensional");
    }
    if (!RaveAttribute_setArrayFromData(attr, PyArray_DATA(arraydata), PyArray_DIM(arraydata, 0), translate_pyarraytype_to_ravetype(PyArray_TYPE(arraydata)))) {
      raiseException_gotoTag(done, PyExc_AttributeError, "Failed to set array data");
    }
  } else {
    raiseException_gotoTag(done, PyExc_AttributeError, "Unsupported data type");
  }

  if (!VerticalProfile_addAttribute(self->vp, attr)) {
    raiseException_gotoTag(done, PyExc_AttributeError, "Failed to add attribute");
  }

  result = PyBool_FromLong(1);
done:
  RAVE_OBJECT_RELEASE(attr);
  return result;
}

/**
 * Returns an attribute with the specified name
 * @param[in] self - this instance
 * @param[in] args - name
 * @returns the attribute value for the name
 */
static PyObject* _pyverticalprofile_getAttribute(PyVerticalProfile* self, PyObject* args)
{
  RaveAttribute_t* attribute = NULL;
  char* name = NULL;
  PyObject* result = NULL;
  if (!PyArg_ParseTuple(args, "s", &name)) {
    return NULL;
  }

  attribute = VerticalProfile_getAttribute(self->vp, name);
  if (attribute != NULL) {
    RaveAttribute_Format format = RaveAttribute_getFormat(attribute);
    if (format == RaveAttribute_Format_Long) {
      long value = 0;
      RaveAttribute_getLong(attribute, &value);
      result = PyLong_FromLong(value);
    } else if (format == RaveAttribute_Format_Double) {
      double value = 0.0;
      RaveAttribute_getDouble(attribute, &value);
      result = PyFloat_FromDouble(value);
    } else if (format == RaveAttribute_Format_String) {
      char* value = NULL;
      RaveAttribute_getString(attribute, &value);
      result = PyString_FromString(value);
    } else if (format == RaveAttribute_Format_LongArray) {
      long* value = NULL;
      int len = 0;
      int i = 0;
      npy_intp dims[1];
      RaveAttribute_getLongArray(attribute, &value, &len);
      dims[0] = len;
      result = PyArray_SimpleNew(1, dims, PyArray_LONG);
      for (i = 0; i < len; i++) {
        *((long*) PyArray_GETPTR1(result, i)) = value[i];
      }
    } else if (format == RaveAttribute_Format_DoubleArray) {
      double* value = NULL;
      int len = 0;
      int i = 0;
      npy_intp dims[1];
      RaveAttribute_getDoubleArray(attribute, &value, &len);
      dims[0] = len;
      result = PyArray_SimpleNew(1, dims, PyArray_DOUBLE);
      for (i = 0; i < len; i++) {
        *((double*) PyArray_GETPTR1(result, i)) = value[i];
      }
    } else {
      RAVE_CRITICAL1("Undefined format on requested attribute %s", name);
      raiseException_gotoTag(done, PyExc_AttributeError, "Undefined attribute");
    }
  } else {
    raiseException_gotoTag(done, PyExc_AttributeError, "No such attribute");
  }
done:
  RAVE_OBJECT_RELEASE(attribute);
  return result;
}

/**
 * Returns if there exists an attribute with the specified name
 * @param[in] self - this instance
 * @param[in] args - name
 * @returns True if attribute exists otherwise False
 */
static PyObject* _pyverticalprofile_hasAttribute(PyVerticalProfile* self, PyObject* args)
{
  char* name = NULL;
  if (!PyArg_ParseTuple(args, "s", &name)) {
    return NULL;
  }
  return PyBool_FromLong((long)VerticalProfile_hasAttribute(self->vp, name));
}

/**
 * Returns a list of attribute names
 * @param[in] self - this instance
 * @param[in] args - N/A
 * @returns a list of attribute names
 */
static PyObject* _pyverticalprofile_getAttributeNames(PyVerticalProfile* self, PyObject* args)
{
  RaveList_t* list = NULL;
  PyObject* result = NULL;
  int n = 0;
  int i = 0;

  list = VerticalProfile_getAttributeNames(self->vp);
  if (list == NULL) {
    raiseException_returnNULL(PyExc_MemoryError, "Could not get attribute names");
  }
  n = RaveList_size(list);
  result = PyList_New(0);
  for (i = 0; result != NULL && i < n; i++) {
    char* name = RaveList_get(list, i);
    if (name != NULL) {
      PyObject* pynamestr = PyString_FromString(name);
      if (pynamestr == NULL) {
        goto fail;
      }
      if (PyList_Append(result, pynamestr) != 0) {
        Py_DECREF(pynamestr);
        goto fail;
      }
      Py_DECREF(pynamestr);
    }
  }
  RaveList_freeAndDestroy(&list);
  return result;
fail:
  RaveList_freeAndDestroy(&list);
  Py_XDECREF(result);
  return NULL;
}

/**
 * All methods a vertical profile can have
 */
static struct PyMethodDef _pyverticalprofile_methods[] =
{
  {"time", NULL},
  {"date", NULL},
  {"starttime", NULL},
  {"endtime", NULL},
  {"startdate", NULL},
  {"enddate", NULL},
  {"source", NULL},
  {"product", NULL},
  {"longitude", NULL},
  {"latitude", NULL},
  {"height", NULL},
  {"interval", NULL},
  {"maxheight", NULL},
  {"setLevels", (PyCFunction) _pyverticalprofile_setLevels, 1},
  {"getLevels", (PyCFunction) _pyverticalprofile_getLevels, 1},
  {"getFF", (PyCFunction) _pyverticalprofile_getFF, 1},
  {"setFF", (PyCFunction) _pyverticalprofile_setFF, 1},
  {"getFFDev", (PyCFunction) _pyverticalprofile_getFFDev, 1},
  {"setFFDev", (PyCFunction) _pyverticalprofile_setFFDev, 1},
  {"getW", (PyCFunction) _pyverticalprofile_getW, 1},
  {"setW", (PyCFunction) _pyverticalprofile_setW, 1},
  {"getWDev", (PyCFunction) _pyverticalprofile_getWDev, 1},
  {"setWDev", (PyCFunction) _pyverticalprofile_setWDev, 1},
  {"getDD", (PyCFunction) _pyverticalprofile_getDD, 1},
  {"setDD", (PyCFunction) _pyverticalprofile_setDD, 1},
  {"getDDDev", (PyCFunction) _pyverticalprofile_getDDDev, 1},
  {"setDDDev", (PyCFunction) _pyverticalprofile_setDDDev, 1},
  {"getDiv", (PyCFunction) _pyverticalprofile_getDiv, 1},
  {"setDiv", (PyCFunction) _pyverticalprofile_setDiv, 1},
  {"getDivDev", (PyCFunction) _pyverticalprofile_getDivDev, 1},
  {"setDivDev", (PyCFunction) _pyverticalprofile_setDivDev, 1},
  {"getDef", (PyCFunction) _pyverticalprofile_getDef, 1},
  {"setDef", (PyCFunction) _pyverticalprofile_setDef, 1},
  {"getDefDev", (PyCFunction) _pyverticalprofile_getDefDev, 1},
  {"setDefDev", (PyCFunction) _pyverticalprofile_setDefDev, 1},
  {"getAD", (PyCFunction) _pyverticalprofile_getAD, 1},
  {"setAD", (PyCFunction) _pyverticalprofile_setAD, 1},
  {"getADDev", (PyCFunction) _pyverticalprofile_getADDev, 1},
  {"setADDev", (PyCFunction) _pyverticalprofile_setADDev, 1},
  {"getDBZ", (PyCFunction) _pyverticalprofile_getDBZ, 1},
  {"setDBZ", (PyCFunction) _pyverticalprofile_setDBZ, 1},
  {"getDBZDev", (PyCFunction) _pyverticalprofile_getDBZDev, 1},
  {"setDBZDev", (PyCFunction) _pyverticalprofile_setDBZDev, 1},
  {"getHGHT", (PyCFunction) _pyverticalprofile_getHGHT, 1},
  {"setHGHT", (PyCFunction) _pyverticalprofile_setHGHT, 1},
  {"getNV", (PyCFunction) _pyverticalprofile_getNV, 1},
  {"setNV", (PyCFunction) _pyverticalprofile_setNV, 1},
  {"getUWND", (PyCFunction) _pyverticalprofile_getUWND, 1},
  {"setUWND", (PyCFunction) _pyverticalprofile_setUWND, 1},
  {"getVWND", (PyCFunction) _pyverticalprofile_getVWND, 1},
  {"setVWND", (PyCFunction) _pyverticalprofile_setVWND, 1},
  {"getFields", (PyCFunction) _pyverticalprofile_getFields, 1},
  {"addField", (PyCFunction) _pyverticalprofile_addField, 1},
  {"getField", (PyCFunction) _pyverticalprofile_getField, 1},
  {"addAttribute", (PyCFunction) _pyverticalprofile_addAttribute, 1},
  {"getAttribute", (PyCFunction) _pyverticalprofile_getAttribute, 1},
  {"hasAttribute", (PyCFunction) _pyverticalprofile_hasAttribute, 1},
  {"getAttributeNames", (PyCFunction) _pyverticalprofile_getAttributeNames, 1},
  {NULL, NULL}
};

/**
 * Returns the specified attribute in the vertical profile
 */
static PyObject* _pyverticalprofile_getattro(PyVerticalProfile* self, PyObject* name)
{
  if (PY_COMPARE_STRING_WITH_ATTRO_NAME("time", name) == 0) {
    const char* str = VerticalProfile_getTime(self->vp);
    if (str != NULL) {
      return PyString_FromString(str);
    } else {
      Py_RETURN_NONE;
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("date", name) == 0) {
    const char* str = VerticalProfile_getDate(self->vp);
    if (str != NULL) {
      return PyString_FromString(str);
    } else {
      Py_RETURN_NONE;
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("source", name) == 0) {
    const char* str = VerticalProfile_getSource(self->vp);
    if (str != NULL) {
      return PyString_FromString(str);
    } else {
      Py_RETURN_NONE;
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("product", name) == 0) {
    const char* str = VerticalProfile_getProduct(self->vp);
    if (str != NULL) {
      return PyString_FromString(str);
    } else {
      Py_RETURN_NONE;
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("starttime", name) == 0) {
    const char* str = VerticalProfile_getStartTime(self->vp);
    if (str != NULL) {
      return PyString_FromString(str);
    } else {
      Py_RETURN_NONE;
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("endtime", name) == 0) {
    const char* str = VerticalProfile_getEndTime(self->vp);
    if (str != NULL) {
      return PyString_FromString(str);
    } else {
      Py_RETURN_NONE;
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("startdate", name) == 0) {
    const char* str = VerticalProfile_getStartDate(self->vp);
    if (str != NULL) {
      return PyString_FromString(str);
    } else {
      Py_RETURN_NONE;
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("enddate", name) == 0) {
    const char* str = VerticalProfile_getEndDate(self->vp);
    if (str != NULL) {
      return PyString_FromString(str);
    } else {
      Py_RETURN_NONE;
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("longitude", name) == 0) {
    return PyFloat_FromDouble(VerticalProfile_getLongitude(self->vp));
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("latitude", name) == 0) {
    return PyFloat_FromDouble(VerticalProfile_getLatitude(self->vp));
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("height", name) == 0) {
    return PyFloat_FromDouble(VerticalProfile_getHeight(self->vp));
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("interval", name) == 0) {
    return PyFloat_FromDouble(VerticalProfile_getInterval(self->vp));
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("minheight", name) == 0) {
    return PyFloat_FromDouble(VerticalProfile_getMinheight(self->vp));
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("maxheight", name) == 0) {
    return PyFloat_FromDouble(VerticalProfile_getMaxheight(self->vp));
  }

  return PyObject_GenericGetAttr((PyObject*)self, name);
}

/**
 * Returns the specified attribute in the vertical profile
 */
static int _pyverticalprofile_setattro(PyVerticalProfile* self, PyObject* name, PyObject* val)
{
  int result = -1;
  if (name == NULL) {
    goto done;
  }
  if (PY_COMPARE_STRING_WITH_ATTRO_NAME("time", name) == 0) {
    if (PyString_Check(val)) {
      if (!VerticalProfile_setTime(self->vp, PyString_AsString(val))) {
        raiseException_gotoTag(done, PyExc_ValueError, "time must be a string (HHmmss)");
      }
    } else if (val == Py_None) {
      VerticalProfile_setTime(self->vp, NULL);
    } else {
      raiseException_gotoTag(done, PyExc_ValueError, "time must be a string (HHmmss)");
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("date", name) == 0) {
    if (PyString_Check(val)) {
      if (!VerticalProfile_setDate(self->vp, PyString_AsString(val))) {
        raiseException_gotoTag(done, PyExc_ValueError, "date must be a string (YYYYMMSS)");
      }
    } else if (val == Py_None) {
      VerticalProfile_setDate(self->vp, NULL);
    } else {
      raiseException_gotoTag(done, PyExc_ValueError, "date must be a string (YYYYMMSS)");
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("source", name) == 0) {
    if (PyString_Check(val)) {
      if (!VerticalProfile_setSource(self->vp, PyString_AsString(val))) {
        raiseException_gotoTag(done, PyExc_ValueError, "source must be a string");
      }
    } else if (val == Py_None) {
      VerticalProfile_setSource(self->vp, NULL);
    } else {
      raiseException_gotoTag(done, PyExc_ValueError, "source must be a string");
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("product", name) == 0) {
    if (PyString_Check(val)) {
      if (!VerticalProfile_setProduct(self->vp, PyString_AsString(val))) {
        raiseException_gotoTag(done, PyExc_ValueError, "product must be a string");
      }
    } else if (val == Py_None) {
      VerticalProfile_setProduct(self->vp, NULL);
    } else {
      raiseException_gotoTag(done, PyExc_ValueError, "product must be a string");
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("starttime", name) == 0) {
    if (PyString_Check(val)) {
      if (!VerticalProfile_setStartTime(self->vp, PyString_AsString(val))) {
        raiseException_gotoTag(done, PyExc_ValueError, "starttime must be a string (HHmmss)");
      }
    } else if (val == Py_None) {
      VerticalProfile_setStartTime(self->vp, NULL);
    } else {
      raiseException_gotoTag(done, PyExc_ValueError, "starttime must be a string (HHmmss)");
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("endtime", name) == 0) {
    if (PyString_Check(val)) {
      if (!VerticalProfile_setEndTime(self->vp, PyString_AsString(val))) {
        raiseException_gotoTag(done, PyExc_ValueError, "endtime must be a string (HHmmss)");
      }
    } else if (val == Py_None) {
      VerticalProfile_setEndTime(self->vp, NULL);
    } else {
      raiseException_gotoTag(done, PyExc_ValueError, "endtime must be a string (HHmmss)");
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("startdate", name) == 0) {
    if (PyString_Check(val)) {
      if (!VerticalProfile_setStartDate(self->vp, PyString_AsString(val))) {
        raiseException_gotoTag(done, PyExc_ValueError, "startdate must be a string (YYYYMMSS)");
      }
    } else if (val == Py_None) {
      VerticalProfile_setStartDate(self->vp, NULL);
    } else {
      raiseException_gotoTag(done, PyExc_ValueError, "startdate must be a string (YYYYMMSS)");
    } 
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("enddate", name) == 0) {
    if (PyString_Check(val)) {
      if (!VerticalProfile_setEndDate(self->vp, PyString_AsString(val))) {
        raiseException_gotoTag(done, PyExc_ValueError, "enddate must be a string (YYYYMMSS)");
      }
    } else if (val == Py_None) {
      VerticalProfile_setEndDate(self->vp, NULL);
    } else {
      raiseException_gotoTag(done, PyExc_ValueError, "enddate must be a string (YYYYMMSS)");
    } 
  }  else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("longitude", name) == 0) {
    if (PyFloat_Check(val)) {
      VerticalProfile_setLongitude(self->vp, PyFloat_AsDouble(val));
    } else {
      raiseException_gotoTag(done, PyExc_TypeError, "longitude must be of type float");
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("latitude", name) == 0) {
    if (PyFloat_Check(val)) {
      VerticalProfile_setLatitude(self->vp, PyFloat_AsDouble(val));
    } else {
      raiseException_gotoTag(done, PyExc_TypeError, "latitude must be of type float");
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("height", name) == 0) {
    if (PyFloat_Check(val)) {
      VerticalProfile_setHeight(self->vp, PyFloat_AsDouble(val));
    } else {
      raiseException_gotoTag(done, PyExc_TypeError, "height must be of type float");
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("interval", name) == 0) {
    if (PyFloat_Check(val)) {
      VerticalProfile_setInterval(self->vp, PyFloat_AsDouble(val));
    } else {
      raiseException_gotoTag(done, PyExc_TypeError, "interval must be of type float");
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("minheight", name) == 0) {
    if (PyFloat_Check(val)) {
      VerticalProfile_setMinheight(self->vp, PyFloat_AsDouble(val));
    } else {
      raiseException_gotoTag(done, PyExc_TypeError, "minheight must be of type float");
    }
  } else if (PY_COMPARE_STRING_WITH_ATTRO_NAME("maxheight", name) == 0) {
    if (PyFloat_Check(val)) {
      VerticalProfile_setMaxheight(self->vp, PyFloat_AsDouble(val));
    } else {
      raiseException_gotoTag(done, PyExc_TypeError, "maxheight must be of type float");
    }
  } else {
    raiseException_gotoTag(done, PyExc_AttributeError, PY_RAVE_ATTRO_NAME_TO_STRING(name));
  }

  result = 0;
done:
  return result;
}

static PyObject* _pyverticalprofile_isVerticalProfile(PyObject* self, PyObject* args)
{
  PyObject* inobj = NULL;
  if (!PyArg_ParseTuple(args,"O", &inobj)) {
    return NULL;
  }
  if (PyVerticalProfile_Check(inobj)) {
    return PyBool_FromLong(1);
  }
  return PyBool_FromLong(0);
}
/*@} End of Vertical Profile */

/// --------------------------------------------------------------------
/// Type definitions
/// --------------------------------------------------------------------
/*@{ Type definition */
PyTypeObject PyVerticalProfile_Type =
{
  PyVarObject_HEAD_INIT(NULL, 0) /*ob_size*/
  "VerticalProfileCore", /*tp_name*/
  sizeof(PyVerticalProfile), /*tp_size*/
  0, /*tp_itemsize*/
  /* methods */
  (destructor)_pyverticalprofile_dealloc, /*tp_dealloc*/
  0, /*tp_print*/
  (getattrfunc)0,               /*tp_getattr*/
  (setattrfunc)0,               /*tp_setattr*/
  0,                            /*tp_compare*/
  0,                            /*tp_repr*/
  0,                            /*tp_as_number */
  0,
  0,                            /*tp_as_mapping */
  0,                            /*tp_hash*/
  (ternaryfunc)0,               /*tp_call*/
  (reprfunc)0,                  /*tp_str*/
  (getattrofunc)_pyverticalprofile_getattro, /*tp_getattro*/
  (setattrofunc)_pyverticalprofile_setattro, /*tp_setattro*/
  0,                            /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT, /*tp_flags*/
  0,                            /*tp_doc*/
  (traverseproc)0,              /*tp_traverse*/
  (inquiry)0,                   /*tp_clear*/
  0,                            /*tp_richcompare*/
  0,                            /*tp_weaklistoffset*/
  0,                            /*tp_iter*/
  0,                            /*tp_iternext*/
  _pyverticalprofile_methods,   /*tp_methods*/
  0,                            /*tp_members*/
  0,                            /*tp_getset*/
  0,                            /*tp_base*/
  0,                            /*tp_dict*/
  0,                            /*tp_descr_get*/
  0,                            /*tp_descr_set*/
  0,                            /*tp_dictoffset*/
  0,                            /*tp_init*/
  0,                            /*tp_alloc*/
  0,                            /*tp_new*/
  0,                            /*tp_free*/
  0,                            /*tp_is_gc*/
};
/*@} End of Type definition */

/*@{ Module setup */
static PyMethodDef functions[] = {
  {"new", (PyCFunction)_pyverticalprofile_new, 1},
  {"isVerticalProfile", (PyCFunction)_pyverticalprofile_isVerticalProfile, 1},
  {NULL,NULL} /*Sentinel*/
};

MOD_INIT(_verticalprofile)
{
  PyObject *module=NULL,*dictionary=NULL;
  static void *PyVerticalProfile_API[PyVerticalProfile_API_pointers];
  PyObject *c_api_object = NULL;

  MOD_INIT_SETUP_TYPE(PyVerticalProfile_Type, &PyType_Type);

  MOD_INIT_VERIFY_TYPE_READY(&PyVerticalProfile_Type);

  MOD_INIT_DEF(module, "_verticalprofile", NULL/*doc*/, functions);
  if (module == NULL) {
    return MOD_INIT_ERROR;
  }

  PyVerticalProfile_API[PyVerticalProfile_Type_NUM] = (void*)&PyVerticalProfile_Type;
  PyVerticalProfile_API[PyVerticalProfile_GetNative_NUM] = (void *)PyVerticalProfile_GetNative;
  PyVerticalProfile_API[PyVerticalProfile_New_NUM] = (void*)PyVerticalProfile_New;

  c_api_object = PyCapsule_New(PyVerticalProfile_API, PyVerticalProfile_CAPSULE_NAME, NULL);
  dictionary = PyModule_GetDict(module);
  PyDict_SetItemString(dictionary, "_C_API", c_api_object);

  ErrorObject = PyErr_NewException("_verticalprofile.error", NULL, NULL);
  if (ErrorObject == NULL || PyDict_SetItemString(dictionary, "error", ErrorObject) != 0) {
    Py_FatalError("Can't define _verticalprofile.error");
    return MOD_INIT_ERROR;
  }

  import_array(); /*To make sure I get access to Numeric*/
  import_pyravefield();
  PYRAVE_DEBUG_INITIALIZE;
  return MOD_INIT_SUCCESS(module);
}
/*@} End of Module setup */
