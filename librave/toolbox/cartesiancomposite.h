/* --------------------------------------------------------------------
Copyright (C) 2013 Swedish Meteorological and Hydrological Institute, SMHI,

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
 * Defines the functions available when creating composites from cartesian products.
 * This object supports \ref #RAVE_OBJECT_CLONE.
 *
 * @file
 * @author Anders Henja (Swedish Meteorological and Hydrological Institute, SMHI)
 * @date 2013-10-01
 */
#ifndef CARTESIANCOMPOSITE_H
#define CARTESIANCOMPOSITE_H
#include "rave_transform.h"
#include "projection.h"
#include "area.h"
#include "rave_object.h"
#include "rave_types.h"
#include "rave_list.h"
#include "raveobject_list.h"
#include "rave_attribute.h"
#include "rave_field.h"
#include "cartesianparam.h"
#include "cartesian.h"
#include "composite.h"

/**
 * Defines a Cartesian composite generator
 */
typedef struct _CartesianComposite_t CartesianComposite_t;

/**
 * Type definition to use when creating a rave object.
 */
extern RaveCoreObjectType CartesianComposite_TYPE;

/**
 * Adds a cartesian product to the cartesian composite generator.
 * @param[in] self - self
 * @param[in] object - the item to be added to the composite
 * @returns 1 on success, otherwise 0
 */
int CartesianComposite_add(CartesianComposite_t* self, Cartesian_t* object);

/**
 * Returns the number of objects this composite will process
 * @param[in] self - self
 * @return the number of objects
 */
int CartesianComposite_getNumberOfObjects(CartesianComposite_t* self);

/**
 * Return the object at position index.
 * @param[in] self - self
 * @param[in] index - the index, should be >= 0 and < getNumberOfObjects
 * @return the object or NULL if outside range
 */
Cartesian_t* CartesianComposite_get(CartesianComposite_t* self, int index);


/**
 * Sets the nominal time.
 * @param[in] self - self
 * @param[in] value - the time in the format HHmmss
 * @returns 1 on success, otherwise 0
 */
int CartesianComposite_setTime(CartesianComposite_t* self, const char* value);

/**
 * Returns the nominal time.
 * @param[in] self - self
 * @returns the nominal time (or NULL if there is none)
 */
const char* CartesianComposite_getTime(CartesianComposite_t* self);

/**
 * Sets the nominal date.
 * @param[in] self - self
 * @param[in] value - the date in the format YYYYMMDD
 * @returns 1 on success, otherwise 0
 */
int CartesianComposite_setDate(CartesianComposite_t* self, const char* value);

/**
 * Returns the nominal date.
 * @param[in] self - self
 * @returns the nominal time (or NULL if there is none)
 */
const char* CartesianComposite_getDate(CartesianComposite_t* self);

/**
 * Sets the quantity
 * @param[in] self - self
 * @param[in] quantity - the quantity, e.g. DBZH
 * @returns 1 on success, otherwise 0
 */
int CartesianComposite_setQuantity(CartesianComposite_t* self, const char* quantity);

/**
 * Returns the quantity (default DBZH)
 * @param[in] self - self
 * @return the quantity
 */
const char* CartesianComposite_getQuantity(CartesianComposite_t* self);

/**
 * Sets the gain.
 * @param[in] self - self
 * @param[in] gain - the gain (MAY NOT BE 0.0)
 */
void CartesianComposite_setGain(CartesianComposite_t* self, double gain);

/**
 * Returns the gain
 * @param[in] self - self
 * @return the gain
 */
double CartesianComposite_getGain(CartesianComposite_t* self);

/**
 * Sets the offset
 * @param[in] self - self
 * @param[in] offset - the offset
 */
void CartesianComposite_setOffset(CartesianComposite_t* self, double offset);

/**
 * Returns the offset
 * @param[in] self - self
 * @return the offset
 */
double CartesianComposite_getOffset(CartesianComposite_t* self);

/**
 * Generates a composite according to the nearest radar principle.
 * @param[in] self - self
 * @param[in] area - the area that should be used for defining the composite.
 * @returns the generated composite.
 */
Cartesian_t* CartesianComposite_nearest(CartesianComposite_t* self, Area_t* area);

#endif /* CARTESIANCOMPOSITE_H */
