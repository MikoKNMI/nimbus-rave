'''
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

Tests the polarvolume module.

@file
@author Anders Henja (Swedish Meteorological and Hydrological Institute, SMHI)
@date 2009-10-14
'''
import unittest
import os
import _polarvolume
import _polarscan
import _rave
import string
import _helpers
import math
import numpy

class PyPolarVolumeTest(unittest.TestCase):
  def setUp(self):
    _helpers.triggerMemoryStatus()

  def tearDown(self):
    pass
    
  def testNewVolume(self):
    obj = _polarvolume.new()
    
    result = string.find(`type(obj)`, "PolarVolumeCore")
    self.assertNotEqual(-1, result) 

  def testVolume_longitude(self):
    obj = _polarvolume.new()
    self.assertAlmostEquals(0.0, obj.longitude, 4)
    obj.longitude = 10.0
    self.assertAlmostEquals(10.0, obj.longitude, 4)

  def testVolume_longitude_typeError(self):
    obj = _polarvolume.new()
    self.assertAlmostEquals(0.0, obj.longitude, 4)
    try:
      obj.longitude = 10
      self.fail("Excepted TypeError")
    except TypeError,e:
      pass
    self.assertAlmostEquals(0.0, obj.longitude, 4)

  def testVolume_latitude(self):
    obj = _polarvolume.new()
    self.assertAlmostEquals(0.0, obj.latitude, 4)
    obj.latitude = 10.0
    self.assertAlmostEquals(10.0, obj.latitude, 4)

  def testVolume_latitude_typeError(self):
    obj = _polarvolume.new()
    self.assertAlmostEquals(0.0, obj.latitude, 4)
    try:
      obj.latitude = 10
      self.fail("Excepted TypeError")
    except TypeError,e:
      pass
    self.assertAlmostEquals(0.0, obj.latitude, 4)

  def testVolume_erroneous_member(self):
    pass

  def testVolume_height(self):
    obj = _polarvolume.new()
    self.assertAlmostEquals(0.0, obj.height, 4)
    obj.height = 10.0
    self.assertAlmostEquals(10.0, obj.height, 4)

  def testVolume_height_typeError(self):
    obj = _polarvolume.new()
    self.assertAlmostEquals(0.0, obj.height, 4)
    try:
      obj.height = 10
      self.fail("Excepted TypeError")
    except TypeError,e:
      pass
    self.assertAlmostEquals(0.0, obj.height, 4)
    
  def testVolume_addScan(self):
    obj = _polarvolume.new()
    scan = _polarscan.new()
    obj.addScan(scan)
    #self.assertEquals(1, obj.getNumberOfScans())

  def test_addScan_navigatorChanged(self):
    obj = _polarvolume.new()
    obj.longitude = 10.0
    scan1 = _polarscan.new()
    scan1.longitude = 5.0

    obj.addScan(scan1)
    self.assertAlmostEquals(10.0, scan1.longitude, 4)

    obj.longitude = 15.0
    self.assertAlmostEquals(15.0, scan1.longitude, 4)

    scan1.longitude = 20.0
    self.assertAlmostEquals(20.0, obj.longitude, 4)

  def test_addScan_refcountIncreaseOnScan(self):
    obj = _polarvolume.new()
    ids = []
    for i in range(50):
      scan = _polarscan.new()
      ids.append(`scan`)
      obj.addScan(scan)
    
    for i in range(obj.getNumberOfScans()):
      scan = obj.getScan(i)
      if `scan` != `obj.getScan(i)`:
        self.fail("Failed to verify scan consistency")


  def testVolume_getNumberOfScans(self):
    obj = _polarvolume.new()
    self.assertEquals(0, obj.getNumberOfScans())
    obj.addScan(_polarscan.new())
    self.assertEquals(1, obj.getNumberOfScans())
    obj.addScan(_polarscan.new())
    self.assertEquals(2, obj.getNumberOfScans())
    
  def testVolume_getScan(self):
    obj = _polarvolume.new()
    scan1 = _polarscan.new()
    scan2 = _polarscan.new()
    
    obj.addScan(scan1)
    obj.addScan(scan2)

    scanresult1 = obj.getScan(0)
    scanresult2 = obj.getScan(1)
    
    self.assertTrue (scan1 == scanresult1)
    self.assertTrue (scan2 == scanresult2)

  def test_getScanClosestToElevation_outside(self):
    obj = _polarvolume.new()
    scan1 = _polarscan.new()
    scan1.elangle = 0.1 * math.pi / 180.0
    obj.addScan(scan1)
    scan2 = _polarscan.new()
    scan2.elangle = 0.5 * math.pi / 180.0
    obj.addScan(scan2)
    scan3 = _polarscan.new()
    scan3.elangle = 2.0 * math.pi / 180.0
    obj.addScan(scan3)
    
    els = [(0.0, 0.1), (0.1, 0.1), (0.2, 0.1), (0.3, 0.1), (0.31, 0.5), (1.0, 0.5), (2.0, 2.0)]

    for el in els:
      elevation = el[0]*math.pi / 180.0
      result = obj.getScanClosestToElevation(elevation, 0)
      self.assertAlmostEquals(el[1], result.elangle*180.0/math.pi, 5)

  def test_getScanClosestToElevation_inside(self):
    obj = _polarvolume.new()
    scan1 = _polarscan.new()
    scan1.elangle = 0.1 * math.pi / 180.0
    obj.addScan(scan1)
    scan2 = _polarscan.new()
    scan2.elangle = 0.5 * math.pi / 180.0
    obj.addScan(scan2)
    scan3 = _polarscan.new()
    scan3.elangle = 2.0 * math.pi / 180.0
    obj.addScan(scan3)
    
    els = [(0.0, None), (0.1, 0.1), (0.2, 0.1), (0.3, 0.1), (0.31, 0.5), (1.0, 0.5), (2.0, 2.0), (2.1, None)]

    for el in els:
      elevation = el[0]*math.pi / 180.0
      result = obj.getScanClosestToElevation(elevation, 1)
      if el[1] == None:
        self.assertTrue(result == None)
      else:
        self.assertAlmostEquals(el[1], result.elangle*180.0/math.pi, 5)

  def testVolume_getNearest(self):
    obj = _polarvolume.new()
    obj.longitude = 12.0 * math.pi/180.0
    obj.latitude = 60.0 * math.pi/180.0
    obj.height = 0.0
    scan1 = _polarscan.new()
    scan1.elangle = 0.1 * math.pi / 180.0
    scan1.rstart = 0.0
    scan1.rscale = 5000.0
    scan1.nodata = 10.0
    scan1.undetect = 11.0
    data = numpy.zeros((100, 120), numpy.uint8)
    scan1.setData(data)
    
    scan2 = _polarscan.new()
    scan2.elangle = 1.0 * math.pi / 180.0
    scan2.rstart = 0.0
    scan2.rscale = 5000.0
    scan2.nodata = 10.0
    scan2.undetect = 11.0    
    data = numpy.ones((100, 120), numpy.uint8)
    scan2.setData(data)
    
    obj.addScan(scan1)
    obj.addScan(scan2)
    
    # Allow outside ranges
    t,v = obj.getNearest((12.0*math.pi/180.0, 60.45*math.pi/180.0), 1000.0, 0)
    self.assertEquals(_rave.RaveValueType_DATA, t)
    self.assertAlmostEquals(1.0, v, 4)

    t,v = obj.getNearest((12.0*math.pi/180.0, 62.00*math.pi/180.0), 1000.0, 0)
    self.assertEquals(_rave.RaveValueType_DATA, t)
    self.assertAlmostEquals(0.0, v, 4)
    
    # Only allow inside ranges
    t,v = obj.getNearest((12.0*math.pi/180.0, 60.45*math.pi/180.0), 1000.0, 1)
    self.assertEquals(_rave.RaveValueType_DATA, t)
    self.assertAlmostEquals(1.0, v, 4)

    t,v = obj.getNearest((12.0*math.pi/180.0, 62.00*math.pi/180.0), 1000.0, 1)
    self.assertEquals(_rave.RaveValueType_NODATA, t)
    
    
  def testSortByElevations_ascending(self):
    obj = _polarvolume.new()
    scan1 = _polarscan.new()
    scan1.elangle = 2.0
    scan2 = _polarscan.new()
    scan2.elangle = 3.0
    scan3 = _polarscan.new()
    scan3.elangle = 1.0
    
    obj.addScan(scan1)
    obj.addScan(scan2)
    obj.addScan(scan3)
    
    obj.sortByElevations(1)
    
    scanresult1 = obj.getScan(0)
    scanresult2 = obj.getScan(1)
    scanresult3 = obj.getScan(2)
    
    self.assertTrue (scan3 == scanresult1)
    self.assertTrue (scan1 == scanresult2)
    self.assertTrue (scan2 == scanresult3)

  def testSortByElevations_descending(self):
    obj = _polarvolume.new()
    scan1 = _polarscan.new()
    scan1.elangle = 2.0
    scan2 = _polarscan.new()
    scan2.elangle = 3.0
    scan3 = _polarscan.new()
    scan3.elangle = 1.0
    
    obj.addScan(scan1)
    obj.addScan(scan2)
    obj.addScan(scan3)
    
    obj.sortByElevations(0)
    
    scanresult1 = obj.getScan(0)
    scanresult2 = obj.getScan(1)
    scanresult3 = obj.getScan(2)
    
    self.assertTrue (scan2 == scanresult1)
    self.assertTrue (scan1 == scanresult2)
    self.assertTrue (scan3 == scanresult3)

  def testIsAscending(self):
    obj = _polarvolume.new()
    scan1 = _polarscan.new()
    scan1.elangle = 0.1
    scan2 = _polarscan.new()
    scan2.elangle = 0.3
    scan3 = _polarscan.new()
    scan3.elangle = 0.5
    obj.addScan(scan1)
    obj.addScan(scan2)
    obj.addScan(scan3)
    
    result = obj.isAscendingScans()
    self.assertEquals(True, result)
    
  def testIsAscending_false(self):
    obj = _polarvolume.new()
    scan1 = _polarscan.new()
    scan1.elangle = 0.1
    scan2 = _polarscan.new()
    scan2.elangle = 0.3
    scan3 = _polarscan.new()
    scan3.elangle = 0.5
    obj.addScan(scan1)
    obj.addScan(scan3)
    obj.addScan(scan2)
    
    result = obj.isAscendingScans()
    self.assertEquals(False, result)
    
  def testIsTransformable(self):
    obj = _polarvolume.new()
    scan1 = _polarscan.new()
    scan1.elangle = 0.1
    scan2 = _polarscan.new()
    scan2.elangle = 0.3
    scan3 = _polarscan.new()
    scan3.elangle = 0.5
    obj.addScan(scan1)
    obj.addScan(scan2)
    obj.addScan(scan3)

    result = obj.isTransformable()
    self.assertEquals(True, result)
    
  def testIsTransformable_noScans(self):
    obj = _polarvolume.new()
    result = obj.isTransformable()
    self.assertEquals(False, result)

  def testIsTransformable_oneScan(self):
    obj = _polarvolume.new()
    scan1 = _polarscan.new()
    scan1.elangle = 0.1
    obj.addScan(scan1)
    result = obj.isTransformable()
    self.assertEquals(True, result)

  def testIsTransformable_descending(self):
    obj = _polarvolume.new()
    scan1 = _polarscan.new()
    scan1.elangle = 0.1
    scan2 = _polarscan.new()
    scan2.elangle = 0.01
    obj.addScan(scan1)
    obj.addScan(scan2)
    result = obj.isTransformable()
    self.assertEquals(False, result)

if __name__ == "__main__":
  #import sys;sys.argv = ['', 'Test.testName']
  unittest.main()