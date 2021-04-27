'''
Created on 21 Jan 2021
Last updated 12 Feb 2021

@author: thomasgumbricht
'''

# Standard library imports

from os import path, makedirs

# Third party imports

# Package application imports

#import geoimagine.support.karttur_dt as mj_dt

#from geoimagine.gis import kt_gis as ktgis

class LayerCommon:
    '''
    '''
    def __init__(self):
        pass
    
    def _SetDOY(self):
        self.datum.doyStr = mj_dt.YYYYDOYStr(self.datum.acqdate)
        self.datum.doy = int(self.datum.doyStr)
        
    def _SetAcqdateDOY(self):
        self.datum.acqdatedoy = mj_dt.DateToYYYYDOY(self.datum.acqdate)
        
    def _SetBounds(self, epsg, minx, miny, maxx, maxy):
        '''
        '''
        
        self.epsg = epsg
        
        self.minx = minx
        
        self.miny = miny
        
        self.maxx = maxx
        
        self.maxy = maxy
        
        self.BoundsPtL = ( (minx,maxy),(maxx,maxy),(maxx,miny), (minx,miny) )
        
    def _Update(self, upD):
        '''
        '''
        for key in upD:
            setattr(self, key, upD[key])
        
    def _Exists(self):
        """checks if the layer file exists; creates the folder path to the layer if non-existant."""
        
        if path.isfile(self.FPN):
        
            self.exists = True
            
            return True
        
        else:
        
            if not path.isdir(self.FP):
            
                makedirs(self.FP)
            
            self.exists = False
            
            return False

class Layer(LayerCommon):
    """Layer is the parentid class for all spatial layers."""
   
    def __init__(self, composition, locusD, datumD): 
        """The constructor expects an instance of the composition class."""
        
        LayerCommon.__init__(self)

        self.comp = composition
       
        self.locus = lambda: None
        
        for key, value in locusD.items():
        
            setattr(self.locus, key, value)

        self.datum = lambda: None
        
        for key, value in datumD.items():
        
            setattr(self.datum, key, value)
        
        if self.datum.acqdate:
        
            self._SetDOY()
            
            self._SetAcqdateDOY()
            
        self.path = lambda: None
        
        self.path.volume = self.comp.volume
                
        self._SetPath()
     
    def _SetPath(self):
        """Sets the complete path to region files"""
        #print ('FUCKING PATH',self.path.volume, self.comp.system, self.comp.source, self.comp.division, self.comp.folder, self.locuspath, self.datum.acqdatestr)

        self.FN = '%(prefix)s_%(prod)s_%(reg)s_%(d)s_%(suf)s%(e)s' %{'prefix':self.comp.prefix,'prod':self.comp.product,'reg':self.locus.locus, 'd':self.datum.acqdatestr, 'suf':self.comp.suffix,'e':self.comp.ext}            
        self.FP = path.join('/Volumes',self.path.volume, self.comp.system, self.comp.source, self.comp.division, self.comp.content, self.locus.path, self.datum.acqdatestr)
        self.FPN = path.join(self.FP,self.FN)
        if ' ' in self.FPN:
            exitstr = 'EXITING region FPN contains space %s' %(self.FPN)
            exit(exitstr)
                    
                  
class VectorLayer(Layer):
    def __init__(self, comp, locusD, datumD): 
        
        Layer.__init__(self, comp, locusD, datumD)
        
        pass
    
        '''
        if not 'shp' in filepath.hdr.lower():
            'Error in hdr for vector file'
            ERRORCHECK
        '''

    def CreateVectorAttributeDef(self,fieldDD): 
        '''
        '''
        
        fieldDefD = {}
        
        self.fieldDefL =[]
        
        for key in fieldDD:
        
            fieldD = fieldDD[key]
            
            if 'width' in fieldD:
            
                width = fieldD['width']
            
            else:
            
                self.width = 8
            
            if 'precision' in fieldD:
                precision = fieldD['precision']
            else:
                precision = 0
            if 'keyfield' in fieldD:
                keyfield = fieldD['keyfield']
            elif 'field' in fieldD:
                keyfield = fieldD['field']
            else:
                keyfield = False
            fieldDefD[key] = {'type':fieldD['type'].lower(), 'width':width, 'precision': precision, 'transfer': fieldD['transfer'].lower(), 'source':fieldD['source'], 'keyfield':keyfield}      
        
        for key in fieldDefD:

            self.fieldDefL.append(ktgis.FieldDef(key,fieldDefD[key]))
     
class RasterLayer(Layer):
    def __init__(self, comp, locusD, datumD): 
        Layer.__init__(self, comp, locusD, datumD)
        
    def _GetRastermetadata(self):
        ''' Get raster meta data from existing file
        '''

        self.spatialRef, self.metadata = ktgis.GetRasterMetaData(self.FPN)
        
        # Set composition cellnull
        self.comp.cellnull = self.metadata.cellnull
        
        # Set composition celltype
        self.comp.celltype = self.metadata.celltype

    def RasterOpenGetFirstLayer(self,**kwargs):
        modeD = {'mode':'read'}
        if kwargs is not None:
            for key, value in kwargs.items():
                modeD[key] = value
                #setattr(self, key, value)
        
        self.DS, self.layer = ktgis.RasterOpenGetFirstLayer(self.FPN,modeD)
        
        self.GetGeoFormatD()
        
    def EmptyLayer(self):
        
        self.layer = lambda:None
        
    def RasterCreateWithFirstLayer(self):
        '''
        '''
        dstDS = ktgis.RasterCreateWithFirstLayer(self.FPN, self)
        
        dstDS._WriteFullArray(self)
        
        # close
        dstDS = None
         
    def GetGeoFormatD(self):
        '''
        '''
        self.geoFormatD = {'lins':self.layer.lins,'cols':self.layer.cols,'projection':self.layer.projection,'geotrans':self.layer.geotrans,'cellsize':self.layer.cellsize}
        
    def SetGeoFormat(self,geoFormatD):
        """Sets the geoFormat
            Expects a dict with {['lins'],['cols'],['projection'],['geotrans'],['cellsize']}
        """ 
        for key, value in geoFormatD.items():
            setattr(self, key, value)
        
    def CreateDSWriteRasterArray(self,**kwargs):
        writeD = {'complete':True, 'of':'GTiff'}
        if kwargs is not None:
            for key, value in kwargs.items():
                writeD[key] = value
        CreateGDALraster
        ktgis.CreateDSWriteRasterArray(self, writeD)
      
    def _RetrieveLayerComp(self,session):
        '''  
        '''
        
        # Set the search to content´ and layerid
        searchItemL = ['content','layerid','product']
        
        queryD = {'content':self.comp.content, 'layerid':self.comp.layerid, 'product':self.comp.product}
        
        return session._RetrieveLayerComp(queryD, searchItemL )
    
    def CopyGeoformatFromSrcLayer(self,otherLayer):
        '''Direct copy of geoformat from srcLayer to dstLayer
        '''
        if not hasattr(self,'layer'):
            
            self.layer = lambda:None
            
        itemL = ['lins','cols','projection','geotrans','cellsize']
        
        for item in itemL:
            
            setattr(self.layer, item, getattr(otherLayer,item))
                 
class RegionLayer(Layer): 
    """layer class for arbitrary layers.""" 
    def __init__(self,comp, location, datum, movieframe = False): 

        """The constructor expects an instance of the composition class."""
        Layer.__init__(self, comp, datum)
        
        self.layertype = 'region'
        self.movieframe = movieframe
        self.location = lambda: None
        
        self.location.regionid = location
        
        #Set the filename and path
        self.SetRegionPath()
        
    def _SetLayerPath(self):
        ERRORCHECK
        self._SetRegionPath()
        
    def _SetRegionPath(self):
        """Sets the complete path to region files"""

        self.FN = '%(prefix)s_%(prod)s_%(reg)s_%(d)s%(suf)s%(e)s' %{'prefix':self.comp.prefix,'prod':self.comp.product,'reg':self.location.regionid, 'd':self.datum.acqdatestr, 'suf':self.comp.suffix,'e':self.comp.ext}            
        if self.movieframe:
            self.FP = path.join(self.comp.mainpath, self.comp.source, self.comp.division, self.comp.folder, self.location.regionid)
        else:
            self.FP = path.join(self.comp.mainpath, self.comp.source, self.comp.division, self.comp.folder, self.location.regionid, self.datum.acqdatestr)

        self.FPN = path.join(self.FP,self.FN)
        if ' ' in self.FPN:
            exitstr = 'EXITING region FPN contains space %s' %(self.FPN)
            exit(exitstr)
   
class TextLayer(Layer):
    def __init__(self, comp, locusD, datumD, filepath):
         
        Layer.__init__(self, comp, locusD, datumD, filepath)
        
        pass
