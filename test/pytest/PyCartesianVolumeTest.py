'''
Copyright (C) 2010 Swedish Meteorological and Hydrological Institute, SMHI,

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

Tests the cartesian volume module.

@file
@author Anders Henja (Swedish Meteorological and Hydrological Institute, SMHI)
@date 2010-07-03
'''
import unittest
import _cartesianvolume
import os
import _cartesian
import _projection
import _rave
import _area
import string
import numpy

class PyCartesianVolumeTest(unittest.TestCase):
  def setUp(self):
    pass

  def tearDown(self):
    pass

  def test_new(self):
    obj = _cartesianvolume.new()
    
    self.assertNotEqual(-1, str(type(obj)).find("CartesianVolumeCore"))

  def test_isCartesianVolume(self):
    obj = _cartesianvolume.new()
    self.assertEqual(True, _cartesianvolume.isCartesianVolume(obj))
    
  def test_isCartesianVolumeFalse(self):
    obj = _cartesian.new()
    self.assertEqual(False, _cartesianvolume.isCartesianVolume(obj))

  def test_attribute_visibility(self):
    attrs = ['areaextent', 'date', 'objectType', 
     'projection', 'source', 'time',
     'xscale', 'xsize', 'yscale', 'ysize']
    obj = _cartesianvolume.new()
    alist = dir(obj)
    for a in attrs:
      self.assertEqual(True, a in alist)

  def test_attributes_from_image(self):
    obj = _cartesianvolume.new()
    obj.xscale = 200.0
    obj.yscale = 200.0
    image = _cartesian.new()
    a = _area.new()
    a.xsize = 10
    a.ysize = 10
    a.xscale = 100.0
    a.yscale = 100.0
    a.extent = (1.0, 2.0, 3.0, 4.0)
    a.projection = _projection.new("x", "y", "+proj=latlong +ellps=WGS84 +datum=WGS84")
    
    image.init(a)
    image.date = "20100101"
    image.time = "100000"
    image.source = "PLC:1234"
    image.product = _rave.Rave_ProductType_CAPPI
    
    self.assertEqual(0, obj.xsize)
    self.assertEqual(0, obj.ysize)
      
    obj.addImage(image)
    self.assertEqual(10, obj.xsize)
    self.assertEqual(10, obj.ysize)
    self.assertEqual(1, obj.getNumberOfImages())

  def test_attributes_to_image(self):
    obj = _cartesianvolume.new()
    obj.xscale = 200.0
    obj.yscale = 200.0
    obj.areaextent = (1.0, 2.0, 3.0, 4.0)
    obj.projection = _projection.new("x", "y", "+proj=latlong +ellps=WGS84 +datum=WGS84")
    obj.date = "20100101"
    obj.time = "100000"
    obj.source = "PLC:1234"
    
    image = _cartesian.new()
    image.product = _rave.Rave_ProductType_CAPPI
    
    obj.addImage(image)
    
    self.assertAlmostEqual(200.0, image.xscale, 4)
    self.assertAlmostEqual(200.0, image.yscale, 4)
    self.assertEqual("20100101", image.date)
    self.assertEqual("100000", image.time)
    self.assertEqual("PLC:1234", image.source)
    self.assertAlmostEqual(1.0, image.areaextent[0], 4)
    self.assertAlmostEqual(2.0, image.areaextent[1], 4)
    self.assertAlmostEqual(3.0, image.areaextent[2], 4)
    self.assertAlmostEqual(4.0, image.areaextent[3], 4)
    
    
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()