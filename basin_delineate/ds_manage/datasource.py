'''
Created on 31 Oct 2020

@author: thomasgumbricht
'''

from osgeo import ogr, osr
import os
import sys

class DS():
    '''
    classdocs
    '''
    
    def __init__(self):
        '''
        Constructor
        '''
        driverName = "ESRI Shapefile" # e.g.: GeoJSON, ESRI Shapefile
        self.driver = ogr.GetDriverByName(driverName)
            
    def DelDs(self,fpn):
        ''' Delete data source file
        '''

        if os.path.exists(fpn):
            self.driver.DeleteDataSource(fpn) 
            
    def OpenDs(self, fpn): 
        ''' Open existing data source
        '''
        if self.params.verbose:
            print ('    opening data source: %s' % (fpn))
        ds = self.driver.Open(fpn, 0) # 0 = read only
            
        if ds is None:
            exitstr =  ('Could not open %s' % (fpn))
            sys.exit(exitstr)
    
        return ds
        layer = ds.GetLayer()
        

        
        '''
        self.spatialRef = layer.GetSpatialRef()
        

        featureCount = layer.GetFeatureCount()
        print ("Number of features in %s: %d" % (os.path.basename(fpn),featureCount))
        '''
        
        
        return layer
        
    def CopyDs(self, fpn):
        ''' Copy an existing data source
        '''   
                
        self.OpenDs(fpn)
                          
        # Copy all features to finalOutletL
        
        finalOutletL = [(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY(), feature.GetField("upstream"), 
                         feature.GetField("xcoord"), feature.GetField("ycoord"), feature.GetField("basin")) for feature in self.layer]
    
        return finalOutletL
      
    def WriteDs(self, fpn, spatialRef, featureD, mouth=False):
        ''' Write data source to file
        '''
        
        driverName = "ESRI Shapefile" # e.g.: GeoJSON, ESRI Shapefile

        self.driver = ogr.GetDriverByName(driverName)
        
        if spatialRef == None:
            spatialRef = osr.SpatialReference()
            spatialRef.ImportFromProj4(self.params.proj4CRS)
            
        if self.params.verbose:
            print ('        Writing distilled outlets to %s' %(fpn))
            print ('        spatialref:', spatialRef)
            
        dstds = self.driver.CreateDataSource(fpn)
        
        # Create the destination dataset layer
        dstlayer = dstds.CreateLayer('outlet', spatialRef, geom_type=ogr.wkbPoint)
        
        dstlayer.CreateField(ogr.FieldDefn('id', ogr.OFTInteger))
        dstlayer.CreateField(ogr.FieldDefn('upstream', ogr.OFTReal))
        dstlayer.CreateField(ogr.FieldDefn('xcoord', ogr.OFTReal))
        dstlayer.CreateField(ogr.FieldDefn('ycoord', ogr.OFTReal))
        if mouth:
            dstlayer.CreateField(ogr.FieldDefn('basin_id', ogr.OFTInteger))
            dstlayer.CreateField(ogr.FieldDefn('mouth_id', ogr.OFTInteger))
            dstlayer.CreateField(ogr.FieldDefn('x_orig', ogr.OFTReal))
            dstlayer.CreateField(ogr.FieldDefn('y_orig', ogr.OFTReal))
            
        defn = dstlayer.GetLayerDefn()
        i = -1
        for f in featureD:
            i +=1
            # Create a new feature (attribute and geometry)
            feat = ogr.Feature(defn)
            feat.SetField('id', i)
            feat.SetField('upstream',featureD[f][2])
            feat.SetField('xcoord',featureD[f][0])
            feat.SetField('ycoord',featureD[f][1])
            if mouth:
                feat.SetField('basin_id',featureD[f][3])
                feat.SetField('mouth_id',featureD[f][4])
                feat.SetField('x_orig',featureD[f][5])
                feat.SetField('y_orig',featureD[f][6])
            
            point = ogr.Geometry(ogr.wkbPoint)
            point.AddPoint(featureD[f][0], featureD[f][1])
            
            feat.SetGeometry(point)
            dstlayer.CreateFeature(feat)
            feat = None  # destroy feat