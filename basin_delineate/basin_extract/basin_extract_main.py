'''
Created on 6 Oct 2020

@author: thomasgumbricht
'''

# imports

from __future__ import division
import os
import sys
#from osgeo import ogr, osr

#from params.be_params_v03 import Params

from params.paramsjson import JsonParams

from params.bex_params import BEXparams



from ds_manage.datasource import DS

# stage 0 is the dem preparation, filling of nodata and pitfilling if required
from basin_extract import basin_extract_stage0 as s0

# stage 1 is the stream topology and hydraulic head calculations
from basin_extract import basin_extract_stage1 as s1

# stage 2
from basin_extract import basin_extract_stage2 as s2
from basin_extract import basin_extract_stage4 as s4
from basin_extract import basin_extract_stage6 as s6
#import shapely.wkt

class ListSrcDstDS(DS):
    '''
    classdocs
    '''
    
    def __init__(self, params):
        '''
        Constructor
        '''
        
        # Set all the input parameters
        self.params = params 
        
        # Set the script path
        self.scriptfp = os.path.join(self.params.srcfp)
        
        if not os.path.exists(self.scriptfp):
            
            infostr = " There are no scripts for stage %(stage)d, should have been \
                 at %(f)s" %{'stage':self.params.stage,'f': self.scriptfp}
                 
            print (infostr)
            
        # Set the path to the new basins and data to same as the source path
        self.basinfp = self.datafp = self.params.srcfp
        
        # Go to the stage for this run
        if self.params.stage == 0:
            
            print("    Files produced at stage 0 (DEM fixing)")
            
            if os.path.exists(self.hydoFillPtDEMFPN_s0):
                infostr = '        DEM with single cell filled pits DONE:\n        %s' %(self.hydoFillPtDEMFPN_s0)
            else:
                infostr = '        DEM with single cell filled pits MISSING'
            print (infostr)
            
            if os.path.exists(self.hydoFillAreaDEMFPN_s0):
                infostr = '        DEM with area filled pits DONE:\n        %s' %(self.hydoFillAreaDEMFPN_s0)
            else:
                infostr = '        DEM with area filled pits MISSING'
            print (infostr)
            
            if os.path.exists(self.fillDEMtiles_s0):
                infostr = '        Polygon file with tiling system for DEM fixing DONE:\n        %s' %(self.fillDEMtiles_s0)
            else:
                infostr = '        Polygon file with tiling system for DEM fixing MISSING'
            print (infostr)

        if self.params.stage == 2:
            
            pass
                     
                
        if self.params.stage == 4:
            
            pass
        
        if self.params.stage == 6:
            
            pass

class BasinExtract(DS):
    '''
    classdocs
    '''
    
    def __init__(self, params):
        '''
        Constructor
        '''
        
        # Initiate the data source (DS) manager
        DS.__init__(self)
        
        # Set all the input parameters
        self.params = params 
        
        # Set the script path
        self.scriptfp = os.path.join(self.params.FPNs.fp)
        
        if not os.path.exists(self.scriptfp):
            
            os.makedirs(self.scriptfp)
            
        # Set the path to the new basins and data to same as the source path
        self.basinfp = self.datafp = self.params.FPNs.fp
        
        # Go to the stage for this run
        if self.params.process.parameters.stage == 0:
            s0.ProcStage0(self.params)
            
        if self.params.process.parameters.stage == 1:
            s1.ProcStage1StreamTopology(self.params)
            
        if self.params.process.parameters.stage == 2:
            s2.ProcStage2BasinOutlets(self.params)
        
        if self.params.process.parameters.stage == 4:
            
            # The target file was deleted if <overwrite> set to True
            if not os.path.exists(self.params.FPNs.BasinOutletFpn_s4):
                
                self.ProcStage4()
                
            else:

                self.ProcStage4()


                                
        if self.params.process.parameters.stage == 6:
            
            self.ProcStage6()
                         
    def ProcStage4(self): 
        '''
        '''
        self.params.stage4datafp = os.path.join(self.datafp, 'stage4')
        
        self.params.stage4scriptfp = os.path.join(self.params.stage4datafp, 'script')
        
        if not os.path.exists(self.params.stage4scriptfp):
            
            os.makedirs(self.params.stage4scriptfp)
         
        s4.CreateGRASSprocess(self.params)
                        
    def ProcStage6(self): 
        '''
        '''

        self.stage4datafp = os.path.join(self.datafp, 'stage6')
        self.stage4scriptfp = os.path.join(self.stage4datafp, 'script')
        if not os.path.exists(self.stage4scriptfp):
            os.makedirs(self.stage4scriptfp)
        
        if not os.path.exists(self.params.FPNs.BasinOutletFpn_s4):
            exitstr ='EXITING: the output file from stage 4 %s is missing.\n' %(self.params.FPNs.BasinOutletFpn_s4)

            exit(exitstr)

        
        # Start stage 6
        s6.CleanBasinPolys(self.params)  
            
def StartUp(jsonFPN):
    '''
    '''
    ProcPar = JsonParams()
    
    processD = ProcPar._JsonObj(jsonFPN)
    
    # Inititate the extra parameters required for basin delineate
    
    ProcBEX = BEXparams()
    
    for k in range(len(processD)):
        
        print ('    ',k, processD[k])
    
        params = processD[k]['PP']
        
        # The processing has only one locus and one period
        
        params.locus = params.srcLocations.locusL[0]
        
        params.datum = params.srcPeriod.datumL[0]
        
        params.alias = ProcPar.jsonParams['userproject']['alias']

        # Get geotrans and translate cellsizes to metric units
        
        ProcBEX._SetGeoTransform(params)
        
        ProcBEX._DestinationDS(params)

        print ('params',params)
        
        #params = Params._JsonObj(jsonFPN)

        if params.process.filecheck:
            
            # list the source and destination files
  
            params.process.verbose = 2

            if params.process.verbose:
                infostr = '    Parameters from json file: %s' %(jsonFPN)
                print (infostr)
                print ('        system:',params.system)
                print ('        region:',params.region)
                
                print ('        alias:',params.alias)
                print ('        delete flag:', params.delete)
                print ('        overwrite flag:',params.overwrite)
                print ('        paramsD:',params)
                print ('        dstpath:',params.dstpath)

                print ('        compD:',params.compD)
                
            ListSrcDstDS(params)
        
        else:
            # Start main process
            
            BasinExtract(params)
    
                     
if __name__ == "__main__":
    
    jsonFPN = '/Volumes/GRASS2020/GRASSproject/COP-DEM/region/basin/nordichydro-ease2n/json/nordichydro-ease2n_drainage_outlets_stage0.json'
    
    jsonFPN = '/Volumes/GRASS2020/GRASSproject/COP-DEM/region/basin/nordichydro-ease2n/json/nordichydro-ease2n_drainage_outlets_stage1.json'

    jsonFPN = '/Volumes/GRASS2020/GRASSproject/COP-DEM/region/basin/nordichydro-ease2n/json/nordichydro-ease2n_drainage_outlets_stage2.json'

    jsonFPN = '/Volumes/GRASS2020/GRASSproject/COP-DEM/region/basin/nordichydro-ease2n/json/nordichydro-ease2n_drainage_outlets_stage4.json'
    
    jsonFPN = '/Volumes/GRASS2020/GRASSproject/COP-DEM/region/basin/nordichydro-ease2n/json/nordichydro-ease2n_drainage_outlets_stage6.json'

    StartUp(jsonFPN)
        