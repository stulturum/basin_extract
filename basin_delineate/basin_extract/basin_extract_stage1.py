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
        self.basinfp = self.datafp = self.params.fp
    
        if self.params.verbose:
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
        GRASSshFN = '%(s)s_grass_find_basin_outlets_stage1.sh' %{'s':self.params.region}

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
        
        if self.params.grassDEM.lower() == 'hydro_fill_dem':
            
            if not os.path.exists(self.params.hydroFillDEMFPN_s0):
                
                exitstr = 'EXITING - hydro_fill_dem output from stage 0 missing:\n    %(src)s' %{'src':self.params.hydroFillDEMFPN_s0}
            
                exit (exitstr)
                
            cmd += '# Import hydrologically corrected DEM from stage 0\n'
            
            cmd += 'r.in.gdal input=%(src)s output=%(dst)s\n\n' %{'src':self.params.hydroFillDEMFPN_s0, 'dst':'hydro_fill_dem'}                        
               
        cmd += '# Multiple flow directions (MFD) watershed analysis\n'
                
        cmd += 'r.watershed -a elevation=%(dem)s accumulation=wshed_upstream drainage=wshed_flowdir stream=wshed_stream threshold=%(th)d --overwrite\n\n' %{'dem':self.params.grassDEM , 'th':self.params.basinCellThreshold}

        cmd += '# MFD with color ramps by removing the "#" sign.\n\n'
                
        cmd += '# Convert MFD to 10 * natural log to get a Byte range\n'
        
        cmd += '# r.mapcalc "wshed_ln_upstream = 10*log(wshed_upstream)" --overwrite\n\n'
        
        cmd += '# Set color ramp\n'
        
        cmd += '# r.colors map=wshed_ln_upstream color=ryb\n\n'
        
        cmd += '# Export as geotiff \n'
        
        cmd += '# r.out.gdal -f input=wshed_ln_upstream format=GTiff type=Byte output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.upstreamLnMFDfpn_s1}
        

        cmd += '# Run stream extract \n'
        
        cmd += 'r.stream.extract elevation=%(dem)s accumulation=wshed_upstream threshold=%(streamThreshold)d mexp=%(mexp)f stream_length=%(streamLength)s\
                 stream_rast=extract_stream stream_vector=streamvect direction=extract_flowdir memory=%(mem)d --overwrite\n\n' %{'dem':self.params.grassDEM,
                    'streamThreshold':self.params.streamThreshold, 'mexp':self.params.mexp, 'streamLength':self.params.streamLength, 
                    'mem':self.params.mem}
                 
        cmd += '# Export the stream vector as GeoJSON \n'
        
        cmd += 'v.out.ogr input=streamvect type=auto format=GeoJSON output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.streamvectFPN_s1}

        cmd += '# Run stream order from r.watereshed streams = rivers\n'
        
        cmd += 'r.stream.order stream_rast=wshed_stream direction=wshed_flowdir elevation=%(dem)s accumulation=wshed_upstream\
                 strahler=riverorder_strahler stream_vect=riverorder_vector memory=%(mem)d --overwrite\n\n' %{'dem':self.params.grassDEM, 
                    'mem':self.params.mem}
         
        cmd += '# Export riverorder_vector\n'
        
        cmd += 'v.out.ogr input=riverorder_vector type=auto format=GeoJSON output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.wshed_riverorder_vector_s1}

        cmd += '# Run stream order from r.stream.extract streams = streams\n'
                        
        cmd += 'r.stream.order stream_rast=extract_stream direction=extract_flowdir elevation=%(dem)s accumulation=wshed_upstream\
                 strahler=streamorder_strahler stream_vect=streamorder_vector memory=%(mem)d --overwrite\n\n' %{'dem':self.params.grassDEM, 
                    'mem':self.params.mem}
                 
        cmd += '# Export extract_streamorder_vector\n'
        
        cmd += 'v.out.ogr input=streamorder_vector type=auto format=GeoJSON output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.extract_streamorder_vector_s1}
                  
        cmd += '# r.stream.distance for getting proximity and hydraulic head related to river and stream channels\n'        

        cmd += 'r.stream.distance stream_rast=extract_stream direction=extract_flowdir elevation=%(dem)s method=downstream\
                 distance=stream_proximity difference=hydraulhead --overwrite  memory=%(mem)d\n\n' %{'dem':self.params.grassDEM,
                 'mem':self.params.mem}    

        cmd += '# Export proximity and hydraulic head related to river and stream channels\n'        

        cmd += 'r.out.gdal -f input=stream_proximity format=GTiff output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.stream_proximity_s1}

        cmd += 'r.out.gdal -f input=hydraulhead format=GTiff output=%(fpn)s --overwrite\n\n' %{'fpn':self.params.hydraulhead_s1}

        GRASSshF = open(self.GRASSshFPN,'w')
        
        GRASSshF.write(cmd)
        
        GRASSshF.close()