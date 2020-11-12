'''
Created on 31 Oct 2020

@author: thomasgumbricht
'''

from pprint import pprint
import params.be_xml as be_xml
import os
import sys
from ds_manage.datasource import DS

class Params(DS):
    '''
    classdocs
    '''
    def __init__(self, xmlFPN):
        
        DS.__init__(self)
        
        # Read params from XML
        d = be_xml.ReadXML(xmlFPN)
        
        # Loop over the dict with params
        for root in d:
            
            self.system = d[root]['userproj']['@system']
            
            self.region = d[root]['userproj']['@tractid']
     
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
    
            paramsD = d[root]['process']['parameters']
            
            self.stage = int(paramsD['@stage'])
            self.adjacentdist = paramsD['@adjacentdist']
            self.outlet = paramsD['@outlet']
            self.distill = paramsD['@distill']
            self.clusteroutlet = paramsD['@clusteroutlet']
            self.verbose = int(paramsD['@verbose'])
            self.proj4CRS = paramsD['@proj4CRS']
            self.watershed = paramsD['@watershed']
            
            self.stage = int(d[root]['process']['parameters']["@stage"])
                
            self.srcpath = '/Volumes/%s' %(d[root]['process']['srcpath']['@volume'])
            
            self.dstpath = '/Volumes/%s' %(d[root]['process']['dstpath']['@volume'])
    
            self.srcCompD = d[root]['process']['srccomp']
            
            self.dstCompD = d[root]['process']['dstcomp']
            
            if self.verbose:
                infostr = '    Reading xml file: %s' %(xmlFPN)
                print (infostr)
                print ('        system:',self.system)
                print ('        region:',self.region)
                print ('        delete flag:', self.delete)
                print ('        overwrite flag:',self.overwrite)
                print ('        paramsD:',paramsD)
                print ('        srcpath:',self.srcpath)
                print ('        dstpath:',self.dstpath)
                print ('        srcCompD:',self.srcCompD)
                print ('        dstCompD:',self.dstCompD)
            
        if self.verbose > 1:
        
            pprint(d)
        
        # Set and control source Datasets
        self.SourceDS()
        
        # Arrange processing setup dependent on stage
        
        if self.stage == 0: #Stage = is the preparation of the initial GRASS processing, only generates a bare-bone script file
            # Not yet implemented
            pass
        
        elif self.stage == 1: # Stage 1 is the distillation of outlet points and the scripting of the extraction of basins
            # Not yet implemented
            pass
        
        # Set the source dataset(s) to use
        self.SourceDS()
        
        # Check for source data sets given the defined method
        self.CheckSrcMethod()

        # Set the dst files, and check if it exists
        self.DestinationDS()
        
        '''
        if stage == 1:
            return (stage, overwrite, paramsD, region, dstfpn, dstremfpn, srcfp, SFDsrcfpn, MFDsrcfpn, dsttxtfn, dstPolygonfpn, basinRoot, dsttxtfpn)
        
        if stage == 2: # Stage 2 removes overlapping basins and optionally removes incomplete basins
            return (stage, overwrite, paramsD, region, dstfpn, dstremfpn, srcfp, SFDsrcfpn, MFDsrcfpn, dsttxtfn, dstPolygonfpn, basinRoot, dsttxtfpn)
        '''
        
    def SourceDS(self):
        ''' Set and check the paths to source files
        '''
        
        # The first part of the data source name (prefix) is given
        
        if 'basin-mouth-outlet-pt' in self.srcCompD:
            
            self.Mouthsrcfpn = self.SetSrcDSFPN('basin-mouth-outlet-pt')
            
        if 'basin-SFD-outlet-pt' in self.srcCompD:
            
            self.SFDsrcfpn = self.SetSrcDSFPN('basin-SFD-outlet-pt')
                
        if 'basin-MFD-outlet-pt' in self.srcCompD:
            
            self.MFDsrcfpn = self.SetSrcDSFPN('basin-MFD-outlet-pt')
            
        if 'shorewall-pt'  in self.srcCompD:
            
            self.shorewallsrcfpn = self.SetSrcDSFPN('shorewall-pt')
    
    def SetSrcDSFPN(self, name):
        ''' Construct full path to source data set
        '''
            
        self.src = self.srcCompD[name]['@source']
        
        self.folder = self.srcCompD[name]['@folder']

        self.srcfp = os.path.join(self.srcpath, 'GRASSproject', self.src, 'region', self.folder, self.region, '0')
        
        self.prefix = self.srcCompD[name]['@prefix']
        
        self.product = self.srcCompD[name]['@product']
        
        self.suffix = self.srcCompD[name]['@suffix'] 
        
        srcfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':self.prefix, 'prod':self.product,
                    'reg':self.region, 'sf':self.suffix}
        
        srcfpn = os.path.join(self.srcfp,srcfn)
        
        if not os.path.exists(srcfpn):
            
            exitstr = 'Source file %s does not exist' %(srcfpn)
            
            sys.exit(exitstr)
        
        return srcfpn

    def CheckSrcMethod(self):
        
        if self.outlet.upper() == 'MOUTH' and not self.Mouthsrcfpn:
            exitstr = 'to run a full width mouth outlet analysis you must include a Mouth point vector file'
            sys.exit(exitstr)
               
        if self.outlet.upper() == 'SFD' and not self.SFDsrcfpn:
            exitstr = 'to run an SFD outlet analysis you must include an SFD point vector file'
            sys.exit(exitstr)
            
        if self.outlet.upper() == 'MFD' and not self.MFDsrcfpn:
            exitstr = 'to run an MFD outlet analysis you must include an MFD point vector file'
            sys.exit(exitstr)
            
        if self.outlet.upper() == 'MOUTH' and not self.Mouthsrcfpn:
            exitstr = 'to run an SFD outlet analysis you must include an SFD point vector file'
            sys.exit(exitstr)
            
        if self.outlet.upper() == 'SFD' and not self.SFDsrcfpn:
            exitstr = 'to run a full width mouth outlet analysis you must include a Mouth point vector file'
            sys.exit(exitstr)
            
        if self.outlet.upper() == 'MFD' and not self.MFDsrcfpn:
            exitstr = 'to run an MFD outlet analysis you must include an MFD point vector file'
            sys.exit(exitstr)

    def DestinationDS(self):
        '''
        '''
        # The destination DS name is by default completely defined from the source DS
        # In case there are no attributes in the tag <basin-outlet> a dummy attribute must be added
        if self.dstCompD['basin-outlet'] == '':
            self.dstCompD['basin-outlet'] = {'dummy':0}
         
        if '@source' in self.dstCompD['basin-outlet'].keys():   
            src = self.dstCompD['basin-outlet']['@source']
        else:
            src = self.src
    
        if '@product' in self.dstCompD['basin-outlet'].keys():
            product = self.dstCompD['basin-outlet']['@product'] 
        else:
            product = self.product
            
        if '@folder' in self.dstCompD['basin-outlet'].keys():
            folder = self.dstCompD['basin-outlet']['@folder']
        else:
            folder = self.folder
    
        # Note that for all intermediate files, the path is GRASSProject, i.e. the 'srcpath'
        self.dstfp = os.path.join(self.srcpath, 'GRASSproject', src, 'region', folder, self.region, '0')
        
    
        if '@prefix' in self.dstCompD['basin-outlet'].keys():
            prefix = self.dstCompD['basin-outlet']['@prefix']
    
        else:
            #Force prefix
            prefix = 'basin-%s-outlet' %(self.outlet.lower()) 
            if self.outlet.lower() == 'sfd' and self.distill.lower() == 'mfd':
                prefix = 'basin-sfdxmfd-outlet' 
            if self.outlet.lower() == 'mfd' and self.distill.lower() == 'mfd':
                prefix = 'basin-mfd2-outlet'
            if self.outlet.lower() == 'mouth' and self.distill.lower() == 'mouth':
                prefix = 'basin-mouth2-outlet'
                
        pfL = prefix.split('-') 
        pfpoly = pfL[0]
        for p in range(1,len(pfL)-1):
            pfpoly = '%s_%s' %(pfpoly, pfL[p])
        self.basinRoot = pfpoly
         
        pfpoly = '%s_%s' %(pfpoly, 'area')
    
        if '@suffix' in self.dstCompD['basin-outlet'].keys():
            suffix = self.dstCompD['basin-outlet']['@suffix']
        else:
            suffix = self.suffix 
        
        dstfn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix, 'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.dstfpn = os.path.join(self.dstfp,dstfn)
        
        self.grassbasefn = '%(pf)s_%(prod)s_%(reg)s_0_%(sf)s' %{'pf':'replaceme', 'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.grassdstfn = prefix.replace('-','_')
        
        dstPolygonfn = '%(pf)s-area_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':prefix, 'prod':product,
                    'reg':self.region, 'sf':suffix} 
        self.dstPolygonfpn = os.path.join(self.dstfp,dstPolygonfn)
        
        # depressionfill shape file
        dstDeprissionfn = '%(pf)s-dem_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'pf':'flowdir', 'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.dstDepressionfpn = os.path.join(self.dstfp,dstDeprissionfn)
        
        dstflowdirptfn = 'flowdir-pt-dem_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.dstflowdirptfpn = os.path.join(self.dstfp,dstflowdirptfn)
        
        dstflowdirareafn = 'flowdir-pt-dem_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        inlandCompDEM = 'inland-comp-demm_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.dstflowdirareafpn = os.path.join(self.dstfp,dstflowdirareafn)
        
        dstflowdirareaptfn = 'flowdir-pt-dem_%(prod)s_%(reg)s_0_%(sf)s.shp' %{'prod':product,
                    'reg':self.region, 'sf':suffix}
        
        self.dstflowdirareaptfpn = os.path.join(self.dstfp,dstflowdirareaptfn)
        
        self.inlandCompDEMFPN = os.path.join(self.dstfp,inlandCompDEM)
        
        if self.overwrite and os.path.exists(self.dstPolygonfpn):
            self.DelDs(self.dstPolygonfpn)
    
        dsttxtfn = '%(pf)s-pt_%(prod)s_%(reg)s_0_%(sf)s.txt' %{'pf':prefix, 'prod':product,
                    'reg':self.region, 'sf':suffix}
        self.dsttxtfpn = os.path.join(self.dstfp,dsttxtfn)
        
        self.dstremfpn = False
        if self.outlet.lower() == 'sfd' and self.distill.lower() == 'mfd':     
            dstremfn = 'rem-%s' %(dstfn)
            self.dstremfpn = os.path.join(self.dstfp,dstremfn)
            
        if self.outlet.lower() == 'mfd' and self.distill.lower() == 'mfd':
            dstremfn = 'rem-%s' %(dstfn)
            self.dstremfpn = os.path.join(self.dstfp,dstremfn)
    
        if self.stage == 1 and os.path.exists(self.dstfpn):
            if self.overwrite or self.delete:
                if self.verbose:
                    infostr = '    Deleting existing shape file %s' %(self.dstfpn)
                    print (infostr)
                self.DelDs(self.dstfpn)
                if self.delete:
                    sys.exit()
            
            elif not os.path.exists(self.dstfp):
                os.makedirs(self.dstfp)