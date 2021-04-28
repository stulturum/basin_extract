'''
Created on 26 Apr 2021

@author: thomasgumbricht

basin extract extra parameters
'''

import os

from math import sqrt

from ds_manage.datasource import DS

class BEXparams(DS):
    ''' Extra parameters required for basin extraction
    '''
    
    def __init__(self):
        '''
        '''
        
        DS.__init__(self)
    
    
    def _SetGeoTransform(self,params):
        '''
        '''
        
        # Set the size and area of a pixel from the original DEM
        if not os.path.exists(params.srcLayerD[params.locus][params.datum]['DEM'].FPN):
            
            exitstr = 'EXITING, the input DEM does not exist: \n     %s ' %(params.srcLayerD[params.locus][params.datum]['DEM'].FPN)
            
            exit(exitstr)

        params.gt = lambda:None()
        
        # Gt GDAL geotransformation data
        gt,params.gt.cols,params.gt.rows = self.GetDsTransform(params.srcLayerD[params.locus][params.datum]['DEM'].FPN)
        
        params.gt.cellsize = (gt[1]-gt[5])/2
        params.gt.cellArea = params.gt.cellsize**2/1000000
        params.gt.nsres = -gt[5]
        params.gt.ewres = gt[1]
        params.gt.north = gt[3]
        params.gt.south = gt[3]+gt[5]*params.gt.rows
        params.gt.west = gt[0]
        params.gt.east = gt[0]+gt[1]*params.gt.cols
        
        if params.process.parameters.adjacentDist <= 0:
            
            params.process.parameters.adjacentDist = 1.01*sqrt(2)*params.gt.cellsize
            
        params.process.parameters.basinAreaThreshold = params.process.parameters.basinCellThreshold*params.gt.cellArea
        
        params.process.parameters.terminalClumpAreakm2Threshold = params.process.parameters.terminalClumpCellThreshold*params.gt.cellArea
        
        
    def _DestinationDS(self,params):
        '''
        '''
        # Set the destination paths
        self.SetDstFP(params)
        
        # Set destination files for stage 0
        self.Stage0_dstDS(params)
        
        # Set destination files for stage 1
        self.Stage1_dstDS(params)
        
        # Set destination files for stage 2
        self.Stage2_dstDS(params)
        
        # Set destination files for stage 4
        self.Stage4_dstDS(params)
        
        # Set destination files for stage 6
        self.Stage6_dstDS(params)
        
        #if self.filecheck:
        #    return
        
        params.FPNs = self
        
        
    def SetDstFP(self,params):
        '''
        '''
        
        #self.systemfp = os.path.join('Volume', params.dstPath.volume, params.procsys.dstsystem, params.compD['system']['source'], params.procsys.dstdivision, self.compD['system']['folder'], params.locus, params.datum)
              
        if 'system' in params.srcLayerD[params.locus][params.datum]:
            
            self.fp = params.srcLayerD[params.locus][params.datum]['system'].FP
            
            self.content = params.srcCompD['system'].content
            
            self.product = params.srcCompD['system'].product
            
            self.suffix = params.srcCompD['system'].suffix 
            
        else:
            
            self.fp = params.srcLayerD[params.locus][params.datum]['DEM'].FP
            
            self.content = params.srcCompD['DEM'].content
            
            self.product = params.srcCompD['DEM'].product
            
            self.suffix = params.srcCompD['DEM'].suffix 
            
    def _SetDSNameParts(self,params,dsbandid):
        '''
        '''
         
        if dsbandid in params.srcCompD:
            # Set the file name components if given

            if hasattr(params.srcCompD[dsbandid], 'source'):
                
                src = params.srcCompD[dsbandid].source
                
            else:
                
                src = self.src
    

            if hasattr(params.srcCompD[dsbandid], 'product'):
                
                product = params.srcCompD[dsbandid].product
                
            else:
                
                product = self.product
    

            if hasattr(params.srcCompD[dsbandid], 'content'):
                
                content = params.srcCompD[dsbandid].content

            else:
                
                content = self.content

            if hasattr(params.srcCompD[dsbandid], 'prefix'):
                
                prefix = params.srcCompD[dsbandid].prefix

            else:
                prefix = dsbandid
                
            if hasattr(params.srcCompD[dsbandid], 'suffix'):
                
                suffix = params.srcCompD[dsbandid].suffix
    
            else:
                
                suffix = self.suffix
                
            return (self.fp, prefix, product, suffix )
                
        else:
            
            return (self.fp, dsbandid, self.product, self.suffix )
        
    def Stage0_dstDS(self,params):
        '''
        '''
        
        # input DEM
        #fp, prefix, product, suffix = self._SetDSNameParts(params,'DEM',params)
         
               
        #inputDEM = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
        #            'reg':params.locus, 'sf':suffix}
        
        self.inputDEMFPN_s0 = params.srcLayerD[params.locus][params.datum]['DEM'].FPN
        
        # terminal-clumps
        fp, prefix, product, suffix = self._SetDSNameParts(params,'terminal-clumps')
                
        terminalClumps_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        
        self.terminalClumpsFPN_s0 = os.path.join(fp,'stage0',terminalClumps_s0)

        
        # inland-comp-DEM
        fp, prefix, product, suffix = self._SetDSNameParts(params,'inland-comp-DEM')
        inlandCompDEM_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}

        self.inlandCompDEMFPN_s0 = os.path.join(fp,inlandCompDEM_s0)
        
        # hydro-fill-pt
        fp, prefix, product, suffix = self._SetDSNameParts(params,'hydro-fill-pt')
        hydroFillPt_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        
        self.hydroFillPtFPN_s0 = os.path.join(fp,hydroFillPt_s0)
        
        # inverted-fill-pt
        fp, prefix, product, suffix = self._SetDSNameParts(params,'inverted-fill-pt')
        invertedFillPt_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        
        self.invertedFillPtFPN_s0 = os.path.join(fp,invertedFillPt_s0)
        
        '''
        # hydro-fill-pt-dem
        fp, prefix, product, suffix = self._SetDSNameParts(params,'hydro-fill-pt-dem')
        hydroFillPtDEM_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':params.locus, 'sf':suffix}
        
        self.hydroFillPtDEMFPN_s0 = os.path.join(fp,hydroFillPtDEM_s0)
        '''
        
        # hydro-fill-area
        fp, prefix, product, suffix = self._SetDSNameParts(params,'hydro-fill-area')
        hydroFillArea_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        
        self.hydroFillAreaFPN_s0 = os.path.join(fp,hydroFillArea_s0)
        
        # inverted-fill-area
        fp, prefix, product, suffix = self._SetDSNameParts(params,'inverted-fill-area')
        invertedFillArea_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        
        self.invertedFillAreaFPN_s0 = os.path.join(fp,invertedFillArea_s0)
        
        '''
        # hydro-fill-area-dem
        fp, prefix, product, suffix = self._SetDSNameParts(params,'hydro-fill-area-dem')
        hydroFillAreaDEM_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':params.locus, 'sf':suffix}
        
        self.hydroFillAreaDEMFPN_s0 = os.path.join(fp,hydroFillAreaDEM_s0)
        '''
        
        # hydro-fill-dem
        fp, prefix, product, suffix = self._SetDSNameParts(params,'hydro-fill-dem')
        hydroFillDEM_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        
        self.hydroFillDEMFPN_s0 = os.path.join(fp,hydroFillDEM_s0)
        
        # inverted-fill-dem
        fp, prefix, product, suffix = self._SetDSNameParts(params,'inverted-fill-dem')
        invertedFillDEM_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        
        self.invertedFillDEMFPN_s0 = os.path.join(fp,invertedFillDEM_s0)
        
        # fill-dem-tiles
        fp, prefix, product, suffix = self._SetDSNameParts(params,'fill-dem-tiles')
        fillDEMtiles_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        
        self.fillDEMtiles_s0 = os.path.join(fp,"stage0",fillDEMtiles_s0)
           
    def Stage1_dstDS(self,params):
        ''' Sets the names of all non internal GRASS files
        '''

        # upstream-ln-MFD upstream natural log transforrmed tif output
        fp, prefix, product, suffix = self._SetDSNameParts(params,'upstream-ln-MFD')
        upstreamLnMFD_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.upstreamLnMFDfpn_s1 = os.path.join(fp,'stage1',upstreamLnMFD_s1)
        
        # streamvectFPN_s1
        fp, prefix, product, suffix = self._SetDSNameParts(params,'streamvect')
        streamvect_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.json' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.streamvectFPN_s1 = os.path.join(fp,'stage1',streamvect_s1)
        
        # riverorder_vector_s1
        fp, prefix, product, suffix = self._SetDSNameParts(params,'river-order')
        riverorder_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.json' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.riverorder_vector_s1 = os.path.join(fp,'stage1',riverorder_s1)
        
        # extract_streamorder_vector_s1
        fp, prefix, product, suffix = self._SetDSNameParts(params,'stream-order')
        streamorder_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.json' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.extract_streamorder_vector_s1 = os.path.join(fp,'stage1',streamorder_s1)

        # stream_proximity_s1
        fp, prefix, product, suffix = self._SetDSNameParts(params,'stream-proxim')
        streamproxim_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.stream_proximity_s1 = os.path.join(fp,'stage1',streamproxim_s1)
        
        # hydraulhead_s1
        fp, prefix, product, suffix = self._SetDSNameParts(params,'hydraulichead')
        hydraulhead_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.hydraulhead_s1 = os.path.join(fp,'stage1',hydraulhead_s1)
        
        # tci_s1
        fp, prefix, product, suffix = self._SetDSNameParts(params,'tci')
        tci_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.tci_s1 = os.path.join(fp,'stage1',tci_s1)
        
        # spi_s1
        fp, prefix, product, suffix = self._SetDSNameParts(params,'spi')
        spi_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.spi_s1 = os.path.join(fp,'stage1',spi_s1)
        
        # basin_s1
        fp, prefix, product, suffix = self._SetDSNameParts(params,'basin')
        basin_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.basin_s1 = os.path.join(fp,'stage1',basin_s1)
        
        # half_basin_s1
        fp, prefix, product, suffix = self._SetDSNameParts(params,'half-basin')
        half_basin_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.half_basin_s1 = os.path.join(fp,'stage1',half_basin_s1)
        
        # rusle_slope_length_s1
        fp, prefix, product, suffix = self._SetDSNameParts(params,'rusle-slope-length')
        rusle_slope_length_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.rusle_slope_length_s1 = os.path.join(fp,'stage1',rusle_slope_length_s1)
        
        # rusle_slope_steepness_s1
        fp, prefix, product, suffix = self._SetDSNameParts(params,'rusle-slope-steepness')
        rusle_slope_steepness_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.rusle_slope_steepness_s1 = os.path.join(fp,'stage1',rusle_slope_steepness_s1)

        
    def Stage2_dstDS(self,params):
        ''' Sets the names of all non internal GRASS files
        '''
        
        # upstream-ln-SFD upstream natural log transforrmed tif output
        fp, prefix, product, suffix = self._SetDSNameParts(params,'upstream-ln-SFD')
        upstreamLnSFD_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.upstreamLnSFDfpn_s2 = os.path.join(fp,'stage2',upstreamLnSFD_s2)
        
        # upstream-ln-MFD upstream natural log transforrmed tif output
        fp, prefix, product, suffix = self._SetDSNameParts(params,'upstream-ln-MFD')
        upstreamLnMFD_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.upstreamLnMFDfpn_s2 = os.path.join(fp,'stage2',upstreamLnMFD_s2)
        
        # lowlevel-outlet-clumps 
        fp, prefix, product, suffix = self._SetDSNameParts(params,'lowlevel-outlet-clumps')
        lowlevelOutletClumps_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.lowlevelOutletClumpsfpn_s2 = os.path.join(fp,'stage2',lowlevelOutletClumps_s2)
        
        # shoreline-outlet-clumps 
        fp, prefix, product, suffix = self._SetDSNameParts(params,'shoreline-outlet-clumps')
        shorelineOutletClumps_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.shorelineOutletClumpsfpn_s2 = os.path.join(fp,'stage2',shorelineOutletClumps_s2)
        
        # all-outlets-pt 
        fp, prefix, product, suffix = self._SetDSNameParts(params,'all-outlets-pt')
        allOutlets_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.allOutletsfpn_s2 = os.path.join(fp,'stage2',allOutlets_s2)
        
        # SFD-outlets-pt 
        fp, prefix, product, suffix = self._SetDSNameParts(params,'SFD-outlets-pt')
        SFDOutlets_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.SFDOutletsfpn_s2 = os.path.join(fp,'stage2',SFDOutlets_s2)
        
        # MFD-outlets-pt 
        fp, prefix, product, suffix = self._SetDSNameParts(params,'MFD-outlets-pt')
        MFDOutlets_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.MFDOutletsfpn_s2 = os.path.join(fp,'stage2',MFDOutlets_s2)
        
        if params.process.parameters.outlet.upper() == 'MOUTH':
            params.Outletsfpn_s2 = self.allOutletsfpn_s2
            
        elif params.process.parameters.outlet.upper() == 'SFD': 
            params.Outletsfpn_s2 = self.SFDOutletsfpn_s2
            
        elif params.process.parameters.outlet.upper() == 'MFD': 
            params.Outletsfpn_s2 = self.MFDOutletsfpn_s2
        
        # thickwall 
        fp, prefix, product, suffix = self._SetDSNameParts(params,'thickwall')
        thickwall_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.thickwallfpn_s2 = os.path.join(fp,'stage2',thickwall_s2)
        
        # mouthshoreline 
        fp, prefix, product, suffix = self._SetDSNameParts(params,'mouthshoreline')
        mouthshoreline_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.mouthshorelinefpn_s2 = os.path.join(fp,'stage2',mouthshoreline_s2)
        
        # shorewall 
        fp, prefix, product, suffix = self._SetDSNameParts(params,'shorewall')
        shorewall_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.shorewallfpn_s2 = os.path.join(fp,'stage2',shorewall_s2)
        
        # shorewall-pt 
        fp, prefix, product, suffix = self._SetDSNameParts(params,'shorewall-pt')
        shorewallpt_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.shorewallptfpn_s2 = os.path.join(fp,'stage2',shorewallpt_s2)

        # shorefill-DEM 
        fp, prefix, product, suffix = self._SetDSNameParts(params,'shorefill-DEM')
        shorefillDEM_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.shorefillDEMfpn_s2 = os.path.join(fp,'stage2',shorefillDEM_s2)
        
    def Stage4_dstDS(self,params): 
        '''
        '''
    
        if params.process.parameters.outlet.upper() == 'MOUTH':
            params.process.parameters.stage4_method = 'allmouths'
            ptdsbandid = 'basin-allmouth-outlets'
            redundantpt = 'basin-allmouth-outlets-redundant'
            polydsbandid = 'basin-allmouth-areas'
            omitpt = 'basin-allmouth-outlets-omitted'
            omitpoly = 'basin-allmouth-areas-omitted'
            remainpt = 'basin-allmouth-outlets-remaining'
            duplicatept = 'basin-allmouth-outlets-duplicate'
            mouthCostgrow = 'basin-allmouth-costgrow'
            mouthRouteDEM = 'basin-allmouth-route-DEM'
            hydroDEM = 'basin-allmouth-hydro-DEM'
            watershedUpdrain = 'basin-allmouth-hydro-DEM-updrain-ln'
            
        if params.process.parameters.outlet.upper() == 'SFD':
            params.process.parameters.stage4_method = 'SFDoutlets'
            ptdsbandid = 'basin-SFDpoint-outlets'
            redundantpt = 'basin-SFDpoint-outlets-redundant'
            polydsbandid = 'basin-SFDpoint-areas'
            omitpt = 'basin-SFDpoint-outlets-omitted'
            omitpoly = 'basin-SFDpoint-areas-omitted'
            remainpt = 'basin-SFDpoint-outlets-remaining'
            duplicatept = 'basin-SFDpoint-outlets-duplicate'
            mouthCostgrow = False
            mouthRouteDEM = False
            hydroDEM = False
            watershedUpdrain = False
            
        if params.process.parameters.outlet.upper() == 'MFD':
            params.process.parameters.stage4_method = 'MFDoutlets'
            ptdsbandid = 'basin-MFDpoint-outlets'
            redundantpt = 'basin-MFDpoint-outlets-redundant'
            polydsbandid = 'basin-MFDpoint-areas'
            omitpt = 'basin-MFDpoint-outlets-omitted'
            omitpoly = 'basin-MFDpoint-areas-omitted'
            remainpt = 'basin-MFDpoint-outlets-remaining'
            duplicatept = 'basin-MFDpoint-outlets-duplicate'
            mouthCostgrow = False
            mouthRouteDEM = False
            hydroDEM = False
            watershedUpdrain = False
            
            
        # outlet points  
        fp, prefix, product, suffix = self._SetDSNameParts(params,ptdsbandid)
        
        BasinOutlet_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.BasinOutletFpn_s4 = os.path.join(fp,BasinOutlet_s4)
        
        # redundantoutlet points  
        fp, prefix, product, suffix = self._SetDSNameParts(params,redundantpt)
        
        RedundantOutlet_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.RedundantOutletFpn_s4 = os.path.join(fp,'stage4',RedundantOutlet_s4)

        # outlet areas  
        fp, prefix, product, suffix = self._SetDSNameParts(params,polydsbandid)
        
        BasinAreas_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.BasinAreasFpn_s4 = os.path.join(fp,'stage4',BasinAreas_s4)
        
        # omtted points  
        fp, prefix, product, suffix = self._SetDSNameParts(params,omitpt)
        
        BasinOutletOmitted_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.BasinOutletOmittedFpn_s4 = os.path.join(fp,'stage4',BasinOutletOmitted_s4)
        
        # omtted areas  
        fp, prefix, product, suffix = self._SetDSNameParts(params,omitpoly)
        
        BasinAreaOmitted_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.BasinAreasOmittedFpn_s4 = os.path.join(fp,'stage4',BasinAreaOmitted_s4)
        
        # Remaining outlet points  
        fp, prefix, product, suffix = self._SetDSNameParts(params,remainpt)
        
        BasinOutletRemain_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.BasinOutletRemainFpn_s4 = os.path.join(fp,'stage4',BasinOutletRemain_s4)

        # duplicate outlet points  
        fp, prefix, product, suffix = self._SetDSNameParts(params,duplicatept)
        
        BasinOutletDuplicate_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':params.locus, 'sf':suffix}
        self.BasinOutletDuplicateFpn_s4 = os.path.join(fp,'stage4',BasinOutletDuplicate_s4)
        
        if mouthCostgrow:
            
            fp, prefix, product, suffix = self._SetDSNameParts(params,mouthCostgrow)
        
            BasinMouthCostGrow_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                        'reg':params.locus, 'sf':suffix}
            self.BasinMouthCostGrowFpn_s4 = os.path.join(fp,'stage4',BasinMouthCostGrow_s4)
            
        if mouthRouteDEM:
            
            fp, prefix, product, suffix = self._SetDSNameParts(params,mouthRouteDEM)
        
            BasinMouthRouteDEM_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                        'reg':params.locus, 'sf':suffix}
            self.BasinMouthRouteDEMFpn_s4 = os.path.join(fp,'stage4',BasinMouthRouteDEM_s4)
            
        if hydroDEM:
            
            fp, prefix, product, suffix = self._SetDSNameParts(params,hydroDEM)
        
            BasinHydroDEM_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                        'reg':params.locus, 'sf':suffix}
            self.BasinHydroDEMFpn_s4 = os.path.join(fp,'stage4',BasinHydroDEM_s4)
            
        if watershedUpdrain:
            
            fp, prefix, product, suffix = self._SetDSNameParts(params,watershedUpdrain)
        
            watershedUpdrain_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                        'reg':params.locus, 'sf':suffix}
            self.watershedUpdrainFpn_s4 = os.path.join(fp,'stage4',watershedUpdrain_s4)
            
    def Stage6_dstDS(self,params): 
        '''
        '''       
        if params.process.parameters.outlet.upper() == 'MOUTH':
            params.process.parameters.stage6_method = 'allmouths'
            basinareas = 'basin-valid-mouth-areas'
            omittedareas = 'basin-omitted-mouth-areas'
            remainings1outlets = 'basin-remain-mouth-s1outlets'
            remainings2outlets = 'basin-remain-mouth-s2outlets'
            duplicates2outlets = 'basin-duplicate-mouth-s2outlets'
            
        if params.process.parameters.outlet.upper() == 'SFD':
            self.stage6_method = 'SFDoutlets'
            basinareas = 'basin-valid-SFD-areas'
            omittedareas = 'basin-omitted-SFD-areas'
            remainings1outlets = 'basin-remain-SFD-s1outlets'
            remainings2outlets = 'basin-remain-SFD-s2outlets'
            duplicates2outlets = 'basin-duplicate-SFD-s2outlets'
              
        if params.process.parameters.outlet.upper() == 'MFD':
            self.stage6_method = 'MFDoutlets'
            basinareas = 'basin-valid-MFD-areas'
            omittedareas = 'basin-omitted-MFD-areas'
            remainings1outlets = 'basin-remain-MFD-s1outets'
            remainings2outlets = 'basin-remain-MFD-s2outlets'
            duplicates2outlets = 'basin-duplicate-MFD-s2outlets'
            
        # Final map of basin areas

        fp, prefix, product, suffix = self._SetDSNameParts(params,basinareas)
        basinareasfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                        'reg':params.locus, 'sf':suffix}
        self.basin_areasFPN_s6 = os.path.join(fp,basinareasfn)
        
        fp, prefix, product, suffix = self._SetDSNameParts(params,omittedareas)
        omittedareasfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                        'reg':params.locus, 'sf':suffix}
        self.basin_omitted_areasFPN_s6 = os.path.join(fp,'stage6',omittedareasfn)

        fp, prefix, product, suffix = self._SetDSNameParts(params,remainings1outlets)
        remainings1outletsfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                        'reg':params.locus, 'sf':suffix}
        self.basin_remain_s1_ptFPN_s6 = os.path.join(fp,'stage6',remainings1outletsfn)
        
        fp, prefix, product, suffix = self._SetDSNameParts(params,remainings2outlets)
        remainings2outletsfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                        'reg':params.locus, 'sf':suffix}
        self.basin_remain_s2_ptFPN_s6 = os.path.join(fp,'stage6',remainings2outletsfn)
        
        fp, prefix, product, suffix = self._SetDSNameParts(params,duplicates2outlets)
        duplicates2outletsfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                        'reg':params.locus, 'sf':suffix}
        self.basin_duplicate_s2_ptFPN_s6 = os.path.join(fp,'stage6',duplicates2outletsfn)