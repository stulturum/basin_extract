'''
Created on 6 Oct 2020

@author: thomasgumbricht
'''

# imports

from __future__ import division

import os

from ds_manage.datasource import DS

#import shapely.wkt

class ProcStage1StreamTopology(DS):
    '''
    classdocs
    '''
    
    def __init__(self, params):
        '''
        Constructor
        '''
        
        # Initiate the data source (DS) manager
        DS.__init__(self)
        
        self.params = params 

        # Set the path to the new basins and data to same as the source path
        self.basinfp = self.datafp = self.params.FPNs.fp
    
        if self.params.process.verbose:
            inforstr = '        Stage 1: Scripting up GRASS to find basin outlets'
            print (inforstr)
            
        self.stage1datafp = os.path.join(self.datafp, 'stage1')
        self.stage1scriptfp = os.path.join(self.stage1datafp, 'script')
        if not os.path.exists(self.stage1scriptfp):
            os.makedirs(self.stage1scriptfp)
    
        self._Grassscript()
        
        cmd = '\n# To run the output script you must have GRASS GIS open in a Terminal window session.\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += 'chmod 755 %(fpn)s\n' %{'fpn': self.GRASSshFPN}
        cmd += '# Then execute the command from your GRASS terminal session:\n'
        cmd += '%(fpn)s\n\n'%{'fpn': self.GRASSshFPN}
        cmd += 'The script produces stream related topology and topography:\n'
        
        print (cmd)
                         
    def _Grassscript(self):    
        '''
        '''
        GRASSshFN = '%(s)s_grass_find_basin_outlets_stage1.sh' %{'s':self.params.locus}

        self.GRASSshFPN = os.path.join(self.stage1scriptfp, GRASSshFN)
        
        cmd = '# Stream topology: GRASS watershed+stream analysis\n'
        cmd += '# Created by the Python package basin_extract\n\n'
                
        cmd += '# To run this script you must have GRASS GIS open in a Terminal window session.\n'
        cmd += '# Change the script to be executable by the command:\n'
        cmd += '# chmod 755 %(fpn)s\n' %{'fpn': self.GRASSshFPN}
        cmd += '# Then execute the command from your GRASS terminal session:\n'
        cmd += '# GRASS 7.x.3x ("region"):~ > %(fpn)s\n\n'%{'fpn': self.GRASSshFPN}

        cmd += '# Create destination path\n'
        
        cmd += 'mkdir -p %s\n\n' % (self.stage1scriptfp)
        
        cmd += '# Set region from the original DEM, if the layer is not at DEM@PERMANENT you have to manually update the layer reference.\n'
        
        cmd += 'g.region raster=%(dem)s\n\n' %{'dem': 'DEM@PERMANENT'}
        
        cmd += '# Mapcalc the DEM into  blocking layer only containing 0s.\n'
    
        cmd += 'r.mapcalc "blocking = 0" --overwrite\n\n'
        
        '''
        if self.params.process.parameters.grassDEM.lower() == 'hydro_fill_dem':
            
            if not os.path.exists(self.params.FPNs.hydroFillDEMFPN_s0):
                
                exitstr = 'EXITING - hydro_fill_dem output from stage 0 missing:\n    %(src)s' %{'src':self.params.hydroFillDEMFPN_s0}
            
                exit (exitstr)
                
            cmd += '# Import hydrologically corrected DEM from stage 0\n'
            
            cmd += 'r.in.gdal input=%(src)s output=%(dst)s\n\n' %{'src':self.params.hydroFillDEMFPN_s0, 'dst':'hydro_fill_dem'}                        
        '''  
             
        cmd += '# Multiple flow directions (MFD) watershed analysis\n'
                
        cmd += 'r.watershed -a elevation=%(dem)s max_slope_length=%(max_slope_length)s blocking=blocking accumulation=MFD_upstream drainage=MFD_flowdir stream=MFD_stream\
            tci=MFD_tci spi=MFD_spi basin=MFD_basin half_basin=MFD_half_basin length_slope=MFD_slope_length\
            slope_steepness=MFD_slope_steepness\
            threshold=%(th)d  memory=%(mem)d\
            --overwrite\n\n' %{'dem': self.params.process.parameters.grassDEM,
                             'max_slope_length':self.params.process.parameters.max_slope_length,
                             'th': self.params.process.parameters.basinCellThreshold,
                             'mem':self.params.process.parameters.MBmemory}

        cmd += '# MFD with color ramps by removing the "#" sign.\n\n'
                
        cmd += '# Convert MFD to 10 * natural log to get a Byte range\n'
        
        cmd += '# r.mapcalc "MFD_ln_upstream = 10*log(MFD_upstream)" --overwrite\n\n'
        
        cmd += '# Set color ramp\n'
        
        cmd += '# r.colors map=MFD_ln_upstream color=ryb\n\n'
        
        cmd += '# Export as geotiff \n'
        
        cmd += '# r.out.gdal -f input=MFD_ln_upstream format=GTiff type=Byte output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.FPNs.upstreamLnMFDfpn_s1}
        

        cmd += '# Run stream extract \n'
        
        cmd += 'r.stream.extract elevation=%(dem)s accumulation=MFD_upstream threshold=%(streamThreshold)d mexp=%(mexp)f stream_length=%(streamLength)s\
                 stream_rast=extract_stream stream_vector=streamvect direction=extract_flowdir memory=%(mem)d --overwrite\n\n' %{'dem':self.params.process.parameters.grassDEM,
                    'streamThreshold': self.params.process.parameters.streamThreshold, 'mexp': self.params.process.parameters.mexp, 
                    'streamLength':self.params.process.parameters.streamLength,'mem': self.params.process.parameters.MBmemory}
                 
        cmd += '# Export the stream vector as GeoJSON \n'
        
        cmd += 'v.out.ogr input=streamvect type=auto format=GeoJSON output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.FPNs.streamvectFPN_s1}

        cmd += '# Run stream order from r.watereshed streams = rivers\n'
        
        cmd += 'r.stream.order stream_rast=MFD_stream direction=MFD_flowdir elevation=%(dem)s accumulation=MFD_upstream\
                 strahler=riverorder_strahler stream_vect=riverorder_vector memory=%(mem)d --overwrite\n\n' %{'dem': self.params.process.parameters.grassDEM, 
                    'mem': self.params.process.parameters.MBmemory}
         
        cmd += '# Export riverorder_vector\n'
        
        cmd += 'v.out.ogr input=riverorder_vector type=auto format=GeoJSON output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.FPNs.riverorder_vector_s1}

        cmd += '# Run stream order from r.stream.extract streams = streams\n'
                        
        cmd += 'r.stream.order stream_rast=extract_stream direction=extract_flowdir elevation=%(dem)s accumulation=MFD_upstream\
                 strahler=streamorder_strahler stream_vect=streamorder_vector memory=%(mem)d --overwrite\n\n' %{'dem':self.params.process.parameters.grassDEM, 
                    'mem':self.params.process.parameters.MBmemory}
                 
        cmd += '# Export extract_streamorder_vector\n'
        
        cmd += 'v.out.ogr input=streamorder_vector type=auto format=GeoJSON output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.FPNs.extract_streamorder_vector_s1}
                  
        cmd += '# r.stream.distance for getting proximity and hydraulic head related to river and stream channels\n'        

        cmd += 'r.stream.distance stream_rast=extract_stream direction=extract_flowdir elevation=%(dem)s method=downstream\
                 distance=stream_proximity difference=hydraulhead --overwrite  memory=%(mem)d\n\n' %{'dem':self.params.process.parameters.grassDEM,
                 'mem':self.params.process.parameters.MBmemory}    

        cmd += '# Export proximity and hydraulic head related to river and stream channels\n'        

        cmd += 'r.out.gdal -f input=stream_proximity format=GTiff output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.FPNs.stream_proximity_s1}

        cmd += 'r.out.gdal -f input=hydraulhead format=GTiff output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.FPNs.hydraulhead_s1}
          
        cmd += '# Export tci and spi parameters\n'        

        cmd += 'r.out.gdal -f input=MFD_tci format=GTiff output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.FPNs.tci_s1}

        cmd += 'r.out.gdal -f input=MFD_spi format=GTiff output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.FPNs.spi_s1}
        
        cmd += '# Export basin and half basin\n'        

        cmd += 'r.out.gdal -f input=MFD_basin format=GTiff output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.FPNs.basin_s1}

        cmd += 'r.out.gdal -f input=MFD_half_basin format=GTiff output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.FPNs.half_basin_s1}

            
        cmd += '# Export RUSLE parameters\n'        

        cmd += 'r.out.gdal -f input=MFD_slope_length format=GTiff output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.FPNs.rusle_slope_length_s1}

        cmd += 'r.out.gdal -f input=MFD_slope_steepness format=GTiff output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.FPNs.rusle_slope_steepness_s1}


        GRASSshF = open(self.GRASSshFPN,'w')
        
        GRASSshF.write(cmd)
        
        GRASSshF.close()