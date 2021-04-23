'''
Created on 31 Oct 2020

@author: thomasgumbricht
'''

from pprint import pprint
import params.be_xml as be_xml
import os
from math import sqrt
from ds_manage.datasource import DS

class Params(DS):
    ''' Convert xml tree dict co class parameters
    '''
    def __init__(self, xmlFPN):

        DS.__init__(self)

        # Read params from XML
        d = be_xml.ReadXML(xmlFPN)

        # Loop over the dict with params
        for root in d:

            self.system = d[root]['userproj']['@system']

            self.region = d[root]['userproj']['@tractid']
            
            self.alias = d[root]['userproj']['@alias']

            delete = d[root]['process']['delete']

            if delete[0] in ['T','t','Y','y']:
                self.delete = True
            else:
                self.delete = False

            overwrite = d[root]['process']['overwrite']

            if overwrite[0] in ['T','t','Y','y']:
                self.overwrite = True
            else:
                self.overwrite = False
                
            filecheck = d[root]['process']['filecheck']

            if filecheck[0] in ['T','t','Y','y']:
                self.filecheck = True
                self.overwrite = False
                self.delete = False
            else:
                self.filecheck = False

            paramsD = d[root]['process']['parameters']

            self.stage = int(paramsD['@stage'])
            
            self.grassDEM = paramsD['@grassDEM']
            
            self.adjacentDist = float(paramsD['@adjacentDist'])
            
            self.outlet = paramsD['@outlet']
            
            self.distill = paramsD['@distill']
            
            self.clusterOutlet = paramsD['@clusterOutlet']
            
            self.basinCellThreshold = int(paramsD['@basinCellThreshold'])
            
            self.minBasinCells = int(paramsD['@minBasinCells'])
            
            
            self.streamThreshold = int(paramsD['@streamThreshold'])
            
            self.mexp = float(paramsD['@mexp'])
            
            self.streamLength = int(paramsD['@streamLength'])
            
            self.watershed = paramsD['@watershed']
            
            self.terminalClumpCellThreshold = float(paramsD['@terminalClumpCellThreshold'])
            
            self.fillDirCells = int(paramsD['@fillDirCells'])
            
            self.fillSQL = paramsD['@fillSQL']
            
            self.invertedFillDirCells = int(paramsD['@invertedFillDirCells'])
            
            self.invertedFillSQL = paramsD['@invertedFillSQL']
            
            self.invertedFillValue = paramsD['@invertedFillValue']
            
            self.tileCellSize = int(paramsD['@tileCellSize'])
            
            self.tileCellOverlap = int(paramsD['@tileCellOverlap'])

            self.proj4CRS = paramsD['@proj4CRS']
            
            self.mem = int(paramsD['@MBmemorySwapInRam'])
            
            self.verbose = int(paramsD['@verbose'])
            

            self.srcpath = '/Volumes/%s' %(d[root]['process']['srcpath']['@volume'])
            
            self.dstpath = '/Volumes/%s' %(d[root]['process']['dstpath']['@volume'])
            
            self.compD = d[root]['process']['comp']
            
            for key in self.compD:
                self.src = self.compD[key]['@source']

                self.folder = self.compD[key]['@folder']
                
                self.product = self.compD[key]['@product']
                
                self.prefix = self.compD[key]['@prefix']
                
                self.suffix = self.compD[key]['@suffix']
                
            self.systemfp = os.path.join(self.srcpath, self.system, self.compD['DEM']['@source'], 'region', self.compD['DEM']['@folder'], self.region, '0')
              
            if 'system' in self.compD:
                self.fp = os.path.join(self.dstpath, 'GRASSproject', self.src, 'region', self.folder, self.alias, '0')
            else:
                os.path.join(self.dstpath, 'GRASSproject', self.compD['DEM']['@source'], 'region', self.compD['DEM']['@folder'], self.alias, '0')
            
            if self.filecheck:
                
                self.verbose = 2

            if self.verbose:
                infostr = '    Parameters from xml file: %s' %(xmlFPN)
                print (infostr)
                print ('        system:',self.system)
                print ('        region:',self.region)
                
                print ('        alias:',self.alias)
                print ('        delete flag:', self.delete)
                print ('        overwrite flag:',self.overwrite)
                print ('        paramsD:',paramsD)
                print ('        dstpath:',self.dstpath)

                print ('        compD:',self.compD)


        if self.verbose > 1:

            pprint(d)
            
        # Set the dst files, and check if it exists
        self._DestinationDS()
        
        # Set the size and area of a pixel from the original DEM
        if not os.path.exists(self.inputDEMFPN_s0):
            exitstr = 'EXITING, the input DEM does not exist: \n     %s ' %(self.inputDEMFPN_s0)
            exit(exitstr)

        # Gt GDAL geotransformation data
        gt,self.cols,self.rows = self.GetDsTransform(self.inputDEMFPN_s0)
        self.cellsize = (gt[1]-gt[5])/2
        self.cellArea = self.cellsize**2/1000000
        self.nsres = -gt[5]
        self.ewres = gt[1]
        self.north = gt[3]
        self.south = gt[3]+gt[5]*self.rows
        self.west = gt[0]
        self.east = gt[0]+gt[1]*self.cols
        
        if self.adjacentDist <= 0:
            self.adjacentDist = 1.01*sqrt(2)*self.cellsize
        self.basinAreaThreshold = self.basinCellThreshold*self.cellArea
        
        self.terminalClumpAreakm2Threshold = self.terminalClumpCellThreshold*self.cellArea
               
    def _DestinationDS(self):
        '''
        '''
        
        # Set destination files for stage 0
        self.Stage0_dstDS()
        
        # Set destination files for stage 1
        self.Stage1_dstDS()
        
        # Set destination files for stage 2
        self.Stage2_dstDS()
        
        # Set destination files for stage 4
        self.Stage4_dstDS()
        
        # Set destination files for stage 6
        self.Stage6_dstDS()
        
        if self.filecheck:
            return
            
    def _SetDSNameParts(self,dsbandid):
        '''
        '''

        if dsbandid in self.compD:
            # Set the file name componets if given

            if '@source' in self.compD[dsbandid].keys():
                src = self.compD[dsbandid]['@source']
            else:
                src = self.src
    
            if '@product' in self.compD[dsbandid].keys():
                product = self.compD[dsbandid]['@product']
            else:
                product = self.product
    
            if '@folder' in self.compD[dsbandid].keys():
                folder = self.compD[dsbandid]['@folder']
            else:
                folder = self.folder

            if '@prefix' in self.compD[dsbandid].keys():
                prefix = self.compD[dsbandid]['@prefix']
    
            else:
                prefix = dsbandid
                
            if '@suffix' in self.compD[dsbandid].keys():
                suffix = self.compD[dsbandid]['@suffix']
    
            else:
                suffix = self.suffix
                
            fp = os.path.join(self.dstpath, 'GRASSproject', src, 'region', folder, self.alias, '0')
            
            return (fp, prefix, product, suffix )
                
        else:
            
            return (self.fp, dsbandid, self.product, self.suffix )
        
    def Stage0_dstDS(self):
        '''
        '''
        
        # input DEM
        fp, prefix, product, suffix = self._SetDSNameParts('DEM')
         
               
        inputDEM = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.inputDEMFPN_s0 = os.path.join(self.systemfp,inputDEM)
        
        # terminal-clumps
        fp, prefix, product, suffix = self._SetDSNameParts('terminal-clumps')
                
        terminalClumps_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.terminalClumpsFPN_s0 = os.path.join(fp,'stage0',terminalClumps_s0)
                
        # terminal-clumps-large
        fp, prefix, product, suffix = self._SetDSNameParts('terminal-clumps-large')
        terminalClumpsLarge_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.terminalClumpsLargeFPN_s0 = os.path.join(fp,'stage0',terminalClumpsLarge_s0)
        
        # terminal-clumps-small
        fp, prefix, product, suffix = self._SetDSNameParts('terminal-clumps-small')
        terminalClumpsSmall_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.terminalClumpsSmallFPN_s0 = os.path.join(fp,'stage0',terminalClumpsSmall_s0)
        
        # inland-comp-DEM
        fp, prefix, product, suffix = self._SetDSNameParts('inland-comp-DEM')
        inlandCompDEM_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}

        self.inlandCompDEMFPN_s0 = os.path.join(fp,inlandCompDEM_s0)
        
        # hydro-fill-pt
        fp, prefix, product, suffix = self._SetDSNameParts('hydro-fill-pt')
        hydroFillPt_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.hydroFillPtFPN_s0 = os.path.join(fp,hydroFillPt_s0)
        
        # inverted-fill-pt
        fp, prefix, product, suffix = self._SetDSNameParts('inverted-fill-pt')
        invertedFillPt_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.invertedFillPtFPN_s0 = os.path.join(fp,invertedFillPt_s0)
        
        '''
        # hydro-fill-pt-dem
        fp, prefix, product, suffix = self._SetDSNameParts('hydro-fill-pt-dem')
        hydroFillPtDEM_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.hydroFillPtDEMFPN_s0 = os.path.join(fp,hydroFillPtDEM_s0)
        '''
        
        # hydro-fill-area
        fp, prefix, product, suffix = self._SetDSNameParts('hydro-fill-area')
        hydroFillArea_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.hydroFillAreaFPN_s0 = os.path.join(fp,hydroFillArea_s0)
        
        # inverted-fill-area
        fp, prefix, product, suffix = self._SetDSNameParts('inverted-fill-area')
        invertedFillArea_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.invertedFillAreaFPN_s0 = os.path.join(fp,invertedFillArea_s0)
        
        '''
        # hydro-fill-area-dem
        fp, prefix, product, suffix = self._SetDSNameParts('hydro-fill-area-dem')
        hydroFillAreaDEM_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.hydroFillAreaDEMFPN_s0 = os.path.join(fp,hydroFillAreaDEM_s0)
        '''
        
        # hydro-fill-dem
        fp, prefix, product, suffix = self._SetDSNameParts('hydro-fill-dem')
        hydroFillDEM_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.hydroFillDEMFPN_s0 = os.path.join(fp,hydroFillDEM_s0)
        
        # inverted-fill-dem
        fp, prefix, product, suffix = self._SetDSNameParts('inverted-fill-dem')
        invertedFillDEM_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.invertedFillDEMFPN_s0 = os.path.join(fp,invertedFillDEM_s0)
        
        # fill-dem-tiles
        fp, prefix, product, suffix = self._SetDSNameParts('fill-dem-tiles')
        fillDEMtiles_s0 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.fillDEMtiles_s0 = os.path.join(fp,"stage0",fillDEMtiles_s0)
        
        
    def Stage1_dstDS(self):
        ''' Sets the names of all non internal GRASS files
        '''

        # upstream-ln-MFD upstream natural log transforrmed tif output
        fp, prefix, product, suffix = self._SetDSNameParts('upstream-ln-MFD')
        upstreamLnMFD_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.upstreamLnMFDfpn_s1 = os.path.join(fp,'stage1',upstreamLnMFD_s1)
        

        # streamvectFPN_s1
        fp, prefix, product, suffix = self._SetDSNameParts('streamvect')
        streamvect_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.json' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.streamvectFPN_s1 = os.path.join(fp,'stage1',streamvect_s1)
        
        # wshed_riverorder_vector_s1
        fp, prefix, product, suffix = self._SetDSNameParts('river-order')
        riverorder_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.json' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.wshed_riverorder_vector_s1 = os.path.join(fp,'stage1',riverorder_s1)
        
        # extract_streamorder_vector_s1
        fp, prefix, product, suffix = self._SetDSNameParts('stream-order')
        streamorder_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.json' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.extract_streamorder_vector_s1 = os.path.join(fp,'stage1',streamorder_s1)

        # stream_proximity_s1
        fp, prefix, product, suffix = self._SetDSNameParts('stream-proxim')
        streamproxim_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.stream_proximity_s1 = os.path.join(fp,'stage1',streamproxim_s1)
        
        # hydraulhead_s1
        fp, prefix, product, suffix = self._SetDSNameParts('hydraulichead')
        hydraulhead_s1 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.hydraulhead_s1 = os.path.join(fp,'stage1',hydraulhead_s1)
        
    def Stage2_dstDS(self):
        ''' Sets the names of all non internal GRASS files
        '''
        
        # upstream-ln-SFD upstream natural log transforrmed tif output
        fp, prefix, product, suffix = self._SetDSNameParts('upstream-ln-SFD')
        upstreamLnSFD_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.upstreamLnSFDfpn_s2 = os.path.join(fp,'stage2',upstreamLnSFD_s2)
        
        # upstream-ln-MFD upstream natural log transforrmed tif output
        fp, prefix, product, suffix = self._SetDSNameParts('upstream-ln-MFD')
        upstreamLnMFD_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.upstreamLnMFDfpn_s2 = os.path.join(fp,'stage2',upstreamLnMFD_s2)
        
        # lowlevel-outlet-clumps 
        fp, prefix, product, suffix = self._SetDSNameParts('lowlevel-outlet-clumps')
        lowlevelOutletClumps_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.lowlevelOutletClumpsfpn_s2 = os.path.join(fp,'stage2',lowlevelOutletClumps_s2)
        
        # shoreline-outlet-clumps 
        fp, prefix, product, suffix = self._SetDSNameParts('shoreline-outlet-clumps')
        shorelineOutletClumps_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.shorelineOutletClumpsfpn_s2 = os.path.join(fp,'stage2',shorelineOutletClumps_s2)
        
        # all-outlets-pt 
        fp, prefix, product, suffix = self._SetDSNameParts('all-outlets-pt')
        allOutlets_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.allOutletsfpn_s2 = os.path.join(fp,'stage2',allOutlets_s2)
        
        # SFD-outlets-pt 
        fp, prefix, product, suffix = self._SetDSNameParts('SFD-outlets-pt')
        SFDOutlets_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.SFDOutletsfpn_s2 = os.path.join(fp,'stage2',SFDOutlets_s2)
        
        # MFD-outlets-pt 
        fp, prefix, product, suffix = self._SetDSNameParts('MFD-outlets-pt')
        MFDOutlets_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.MFDOutletsfpn_s2 = os.path.join(fp,'stage2',MFDOutlets_s2)
        
        if self.outlet.upper() == 'MOUTH':
            self.Outletsfpn_s2 = self.allOutletsfpn_s2
            
        elif self.outlet.upper() == 'SFD': 
            self.Outletsfpn_s2 = self.SFDOutletsfpn_s2
            
        elif self.outlet.upper() == 'MFD': 
            self.Outletsfpn_s2 = self.MFDOutletsfpn_s2
        
        # thickwall 
        fp, prefix, product, suffix = self._SetDSNameParts('thickwall')
        thickwall_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.thickwallfpn_s2 = os.path.join(fp,'stage2',thickwall_s2)
        
        # mouthshoreline 
        fp, prefix, product, suffix = self._SetDSNameParts('mouthshoreline')
        mouthshoreline_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.mouthshorelinefpn_s2 = os.path.join(fp,'stage2',mouthshoreline_s2)
        
        # shorewall 
        fp, prefix, product, suffix = self._SetDSNameParts('shorewall')
        shorewall_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.shorewallfpn_s2 = os.path.join(fp,'stage2',shorewall_s2)
        
        # shorewall-pt 
        fp, prefix, product, suffix = self._SetDSNameParts('shorewall-pt')
        shorewallpt_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.shorewallptfpn_s2 = os.path.join(fp,'stage2',shorewallpt_s2)
        
        # shorefill-DEM 
        fp, prefix, product, suffix = self._SetDSNameParts('shorefill-DEM')
        shorefillDEM_s2 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.shorefillDEMfpn_s2 = os.path.join(fp,'stage2',shorefillDEM_s2)
        
    def Stage4_dstDS(self): 
        '''
        '''
    
        if self.outlet.upper() == 'MOUTH':
            self.stage4_method = 'allmouths'
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
            
        if self.outlet.upper() == 'SFD':
            self.stage4_method = 'SFDoutlets'
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
            
        if self.outlet.upper() == 'MFD':
            self.stage4_method = 'MFDoutlets'
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
        fp, prefix, product, suffix = self._SetDSNameParts(ptdsbandid)
        
        BasinOutlet_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.BasinOutletFpn_s4 = os.path.join(fp,BasinOutlet_s4)
        
        # redundantoutlet points  
        fp, prefix, product, suffix = self._SetDSNameParts(redundantpt)
        
        RedundantOutlet_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.RedundantOutletFpn_s4 = os.path.join(fp,'stage4',RedundantOutlet_s4)

        # outlet areas  
        fp, prefix, product, suffix = self._SetDSNameParts(polydsbandid)
        
        BasinAreas_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.BasinAreasFpn_s4 = os.path.join(fp,'stage4',BasinAreas_s4)
        
        # omtted points  
        fp, prefix, product, suffix = self._SetDSNameParts(omitpt)
        
        BasinOutletOmitted_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.BasinOutletOmittedFpn_s4 = os.path.join(fp,'stage4',BasinOutletOmitted_s4)
        
        # omtted areas  
        fp, prefix, product, suffix = self._SetDSNameParts(omitpoly)
        
        BasinAreaOmitted_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.BasinAreasOmittedFpn_s4 = os.path.join(fp,'stage4',BasinAreaOmitted_s4)
        
        # Remaining outlet points  
        fp, prefix, product, suffix = self._SetDSNameParts(remainpt)
        
        BasinOutletRemain_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.BasinOutletRemainFpn_s4 = os.path.join(fp,'stage4',BasinOutletRemain_s4)

        # duplicate outlet points  
        fp, prefix, product, suffix = self._SetDSNameParts(duplicatept)
        
        BasinOutletDuplicate_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.BasinOutletDuplicateFpn_s4 = os.path.join(fp,'stage4',BasinOutletDuplicate_s4)
        
        if mouthCostgrow:
            
            fp, prefix, product, suffix = self._SetDSNameParts(mouthCostgrow)
        
            BasinMouthCostGrow_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                        'reg':self.region, 'sf':suffix}
            self.BasinMouthCostGrowFpn_s4 = os.path.join(fp,'stage4',BasinMouthCostGrow_s4)
            
        if mouthRouteDEM:
            
            fp, prefix, product, suffix = self._SetDSNameParts(mouthRouteDEM)
        
            BasinMouthRouteDEM_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                        'reg':self.region, 'sf':suffix}
            self.BasinMouthRouteDEMFpn_s4 = os.path.join(fp,'stage4',BasinMouthRouteDEM_s4)
            
        if hydroDEM:
            
            fp, prefix, product, suffix = self._SetDSNameParts(hydroDEM)
        
            BasinHydroDEM_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                        'reg':self.region, 'sf':suffix}
            self.BasinHydroDEMFpn_s4 = os.path.join(fp,'stage4',BasinHydroDEM_s4)
            
        if watershedUpdrain:
            
            fp, prefix, product, suffix = self._SetDSNameParts(watershedUpdrain)
        
            watershedUpdrain_s4 = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.tif' %{'pf':prefix,'prod':product,
                        'reg':self.region, 'sf':suffix}
            self.watershedUpdrainFpn_s4 = os.path.join(fp,'stage4',watershedUpdrain_s4)
            
    def Stage6_dstDS(self): 
        '''
        '''       
        if self.outlet.upper() == 'MOUTH':
            self.stage6_method = 'allmouths'
            basinareas = 'basin-valid-mouth-areas'
            omittedareas = 'basin-omitted-mouth-areas'
            remainings1outlets = 'basin-remain-mouth-s1outlets'
            remainings2outlets = 'basin-remain-mouth-s2outlets'
            duplicates2outlets = 'basin-duplicate-mouth-s2outlets'
            
        if self.outlet.upper() == 'SFD':
            self.stage6_method = 'SFDoutlets'
            basinareas = 'basin-valid-SFD-areas'
            omittedareas = 'basin-omitted-SFD-areas'
            remainings1outlets = 'basin-remain-SFD-s1outlets'
            remainings2outlets = 'basin-remain-SFD-s2outlets'
            duplicates2outlets = 'basin-duplicate-SFD-s2outlets'
              
        if self.outlet.upper() == 'MFD':
            self.stage6_method = 'MFDoutlets'
            basinareas = 'basin-valid-MFD-areas'
            omittedareas = 'basin-omitted-MFD-areas'
            remainings1outlets = 'basin-remain-MFD-s1outets'
            remainings2outlets = 'basin-remain-MFD-s2outlets'
            duplicates2outlets = 'basin-duplicate-MFD-s2outlets'
            
        # Final map of basin areas

        fp, prefix, product, suffix = self._SetDSNameParts(basinareas)
        basinareasfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                        'reg':self.region, 'sf':suffix}
        self.basin_areasFPN_s6 = os.path.join(fp,basinareasfn)
        
        fp, prefix, product, suffix = self._SetDSNameParts(omittedareas)
        omittedareasfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                        'reg':self.region, 'sf':suffix}
        self.basin_omitted_areasFPN_s6 = os.path.join(fp,'stage6',omittedareasfn)

        fp, prefix, product, suffix = self._SetDSNameParts(remainings1outlets)
        remainings1outletsfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                        'reg':self.region, 'sf':suffix}
        self.basin_remain_s1_ptFPN_s6 = os.path.join(fp,'stage6',remainings1outletsfn)
        
        fp, prefix, product, suffix = self._SetDSNameParts(remainings2outlets)
        remainings2outletsfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                        'reg':self.region, 'sf':suffix}
        self.basin_remain_s2_ptFPN_s6 = os.path.join(fp,'stage6',remainings2outletsfn)
        
        fp, prefix, product, suffix = self._SetDSNameParts(duplicates2outlets)
        duplicates2outletsfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix,'prod':product,
                        'reg':self.region, 'sf':suffix}
        self.basin_duplicate_s2_ptFPN_s6 = os.path.join(fp,'stage6',duplicates2outletsfn)