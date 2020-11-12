'''
Created on 6 Oct 2020

@author: thomasgumbricht
'''

#imports

from __future__ import division
import os
import sys
#from osgeo import ogr, osr
import numpy as np
from scipy.spatial.distance import cdist
from operator import itemgetter
from params.be_params import Params
from ds_manage.datasource import DS
#import shapely.wkt

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
        self.params = params 
        
        if self.params.stage == 1:
        
            if not os.path.exists(params.dstfpn):
                self.ProcStage1()
            else:
                self.GetDistilled()
            
            self.scriptfp = os.path.join(self.params.srcfp)
            if not os.path.exists(self.scriptfp):
                os.makedirs(self.scriptfp)
            
            # Set the path to the new basns to same as the source path
            self.basinfp = self.params.srcfp
            #stats = False
            #discharge = False
            #iteration = 'qstn'
            #iteration = '02'
            #version = 'v993'
            
        
            #set the filenames
            #coordLFN = '%(s)s_%(typ)s_outlet_test.txt' %{'s':region, 'typ':basinRoot}
            
            #txtFN = '%(s)s_%(typ)s_outlet_final.txt' %{'s':region, 'typ':basinRoot}
            
            self.GRASSshFN = '%(s)s_grass_%(typ)s.sh' %{'s':self.params.region, 'typ':self.params.basinRoot}
            
            self. OGRshFN = '%(s)s_ogr2ogr_%(typ)s.sh' %{'s':self.params.region, 'typ':self.params.basinRoot}
            
            #WriteGrassCmd(finalOutlets, scriptfp, GRASSshFN, OGRshFN, dstPolygonfpn, basinRoot, dsttxtfpn, srcfp, overwrite, paramsD)

            self.WriteGrassCmd()
            
    def ProcStage1(self): 
        '''
        '''

        if self.params.outlet.upper() == 'MOUTH':
            if self.params.distill.upper() == 'MOUTH':
                if self.params.verbose:
                    infostr = '    Identifying outlets from full width mouth data distilled using mouth and basin clusters'
                    print (infostr)

                self.MouthMouthOutlets()
                
            else:
                # copies all input mouth cells to become basin outlet points
                if self.params.verbose:
                    infostr = '    Outlets set equal to full width mouth data'
                    print (infostr)
                FIX
                self.CopyDs()

        elif self.params.outlet == 'SFD':
            if self.params.distill == 'MFD':
                if self.params.verbose:
                    infostr = '    Identifying outlets from SFD data distilled using MFD clusters'
                    print (infostr)
                    
                #finalOutlets, removedOutlets, spatialRef = SFDMFDoutlets(SFDsrcfpn, MFDsrcfpn, verbose, paramsD)
                
            else:
                # This is alternative is simply using the existing SFD
                if self.params.verbose:
                    infostr = '    Outlets set equal to SFD input data'
                    print (infostr)
                #finalOutlets, spatialRef = CopyDs(SFDsrcfpn, paramsD['@proj4CRS'])
                
        elif self.params.outlet == 'MFD':
            if self.params.distill == 'MFD':
                if self.params.verbose:
                    infostr = '    Identifying outlets by distilling MFD clusters to unique outlets'
                    print (infostr)
                #finalOutlets, removedOutlets, spatialRef = MFDoutlets(MFDsrcfpn, verbose, paramsD)
            else:
                if self.params.verbose:
                    infostr = '    Outlets set equal to MFD input data'
                    print (infostr)
                # This is alternative is simply using the existing MFD
                #finalOutlets, spatialRef = CopyDs(MFDsrcfpn, paramsD['@proj4CRS']) 
        else:
            exitstr = 'EXITING, the parameter outlet must be set to either mouth, MFD or SFD' 
            exit(exitstr)
                
        if self.params.verbose:
            print ('    Total nr of final outlets', len(self.finalOutletsD))
                
        self.WriteDs(self.params.dstfpn, self.spatialRef, self.finalOutletsD, True)
   
        if self.dstremfpn:
            self.WriteDs(self.dstremfpn, self.spatialRef, self.removedOutlets)

    def MouthMouthOutlets(self):
        ''' Distill one outlet point per full width mouth
        '''
        
        # Open the dataset with outlet points
        
        ds = self.OpenDs(self.params.Mouthsrcfpn)
        
        srcLayer = ds.GetLayer()
        
        # Open the shore wall points
        
        swds = self.OpenDs(self.params.shorewallsrcfpn)
        
        swLayer = swds.GetLayer()
        
        self.GetSpatialRef(srcLayer)
        
        self.featureCount = srcLayer.GetFeatureCount()
        if self.params.verbose:
            print ("Number of candidate outlets: %d" % (self.featureCount))

        mouthtab = 'mouth_id'
        basintab = 'basin_id'  
        mouthtab = 'sepmouth'
        basintab = 'linkmouth'   
        self.UniqueDbPt(srcLayer, swLayer, mouthtab, basintab)

    def SFDMFDoutlets(self, SFDsrcfpn, MFDsrcfpn):
        ''' Distill SFD outlets from MFD clusters
        '''
        
        if self.params.verbose: 
            infostr = '    Cleaning candidate SFD outlets from clustered MFD outlets'
            print (infostr)
           
        # Open the SFD DS and get the point layer
        SFDsrclayer = self.OpenDs(SFDsrcfpn)
        
        # Open the MFD DS and get the point layer
        MFDsrclayer = self.OpenDs(MFDsrcfpn)
        
        # Get all outlet points that are clustered in the MFD layer
        clusteredD = self.ClusteredOutlets()
        
        # get the complete feature data for the SFD outlets
        self.featsL = [(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY(), feature.GetField("upstream"),
                   feature.GetField("xcoord"), feature.GetField("ycoord"), feature.GetField("basin")) for feature in SFDsrclayer]
        
        if self.params.verbose: 
            infostr = '    Initial number of SFD outlets: %s' %(len(self.featsL))
            print (infostr)
    
        # For each cluster, get all the SFD within the cluster, and only retain the one with the largest upstream
        SFDremoveL = []
        for cluster in clusteredD:
            SFDclusterL = []
            for cell in clusteredD[cluster]:
                for x, slot in enumerate(self.featsL):
                    if cell[0] == slot[0] and cell[1] == slot[1]:
                        SFDclusterL.append( ( x, slot[0], slot[1], slot[2], slot[3] ) )
            if len(SFDclusterL)> 1:
                clusterL = sorted(SFDclusterL, key=itemgetter(3), reverse=True)
                
                for i, item in enumerate(clusterL):
    
                    if i > 0:
                        SFDremoveL.append(item[0])
    
        SFDremoveL.sort(reverse=True)
            
        self.removeL = [self.featsL.pop(x) for x in SFDremoveL]
    
        if self.params.verbose:
            infostr = '    Removing %s SFD outlets identified as clustered from MFD outlets' %(len(SFDremoveL))
            print (infostr)
                                               
    def MFDoutlets(self):
        
                    
        srcLayer = self.OpenDs(self.MFDsrcfpn) 
        
        # Get the closest point of each candidate basin outlet as a list of np.array [x y]    
        self.finalOutletL, self.removedOutlets = self.UniqueOutlets(srcLayer)
    
        # Save and close everything
        ds = layer = feat = geom = None
     
    def UniqueDbPt(self, srcLayer, swLayer, mouthtab, basintab):
        '''
        '''
        
        # get the complete feature data for the input layer
        self.featsL = [(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY(), feature.GetField("upstream"),
                    feature.GetField(mouthtab), feature.GetField(basintab)) for feature in srcLayer]

        #shorewallL = [(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY()) for feature in swLayer]
 
        # Create a numpy array of the shorewall points
        
        [(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY()) for feature in srcLayer]
        
        #li = [ [feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY() ] for feature in swLayer] 

        swCoordA = np.asarray(  [ [feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY() ] for feature in swLayer] )
        
        mouthD = {}
        outletD = {}
        self.widestMouth = 0
        for f in self.featsL:
            if not f[3] in mouthD:
                mouthD[f[3]] = []
            mouthD[f[3]].append(f)
        
        if self.params.verbose:
                
            print ('        identified nr of mouths', len(mouthD))
        
        for m in mouthD:
            if len(mouthD[m]) == 1:
                outletD[m] = mouthD[m][0]
            else:
                if len(mouthD[m]) > self.widestMouth:
                    self.widestMouth = len(mouthD[m])
                if self.params.clusteroutlet[0].lower() == 'c':
                    # Construct a numpy array of coords
                    coordA = np.asarray( [ [p[0],p[1]] for p in mouthD[m] ] )
                    
                    # Calculate the average outlet pt
                    avgXY = np.average(coordA, axis=0)
                    
                    # Find the mouth point closest to the average
                    c = cdist([avgXY], coordA).argmin()
                    
                    # Set the outlet to the point closest to the average
                    outletD[m] = mouthD[m][c]
                else:
                    u = 0
                    for p in mouthD[m]: 
                        if p[2] > u:
                            u = p[2]
                            outletD[m] = p
        
        if self.params.verbose:
            print ('    Identified %s outlet points' %(len(outletD)))
            print ('    Widest mouth = %s pixels' %(self.widestMouth))
            
        # Find the point in the wall that is closest to the original outlet, and move the point
        swD = {}
        for o in outletD:
            oa = np.asarray( [ [outletD[o][0],outletD[o][1]] ] )
            n = cdist(oa, swCoordA).argmin()
            swD[o] = [swCoordA[n][0],swCoordA[n][1],outletD[o][2],outletD[o][3],outletD[o][4],outletD[o][0],outletD[o][1] ]

        self.finalOutletsD = swD
        self.dstremfpn = self.removedOutletsD = 0
             
    def ClusterAdjacent(self):
        ''' Cluster adjacent cells 
        '''
        
        # Convert the list of coordinates to an array
        adjA = np.array(self.adjacentL)

        # Create a Dict to hold the clusters of candidate basin outlet cells
        self.clusterD = {}
        
        # Create a counter for the Dict clusterD
        nrKeys = 0
        
        # Get the total number of candidate basin outlets
        totalAdjCells = adjA.shape[0]
        
        if self.params.verbose:
            print ('    Total number of cells to cluster', totalAdjCells)
            
        # Loop all candidate nodes
        for i,node in enumerate(adjA):
            if self.params.verbose > 1:
                print ('    Processing cell nr', i)
                
            # Get all nodes, except the node itself
            nodes = np.delete(adjA, [i], 0)
                    
            # Retrieve the cells that fall within dt as a Boolean array
            within = cdist([node], nodes) < self.params.thresholddist 
            
            # Mask the nodes (original coordinates) with the Boolean array
            adjacents = nodes[within[0],:]
    
            # Append the original node to the array of nodes adjacent to the source node of this loop
            aA = np.append([node], adjacents, axis=0)
    
            # Loop cluster Dicts to check if this cluster is already registered
            cellin = False
            for cell in aA:
                for key in self.clusterD:  
                    for node in self.clusterD[key]:
                        # check if the cell and the node are identical
                        if cell[0] == node[0] and cell[1] == node[1]:
                            # if a single cell of this new cluster is already registerd, the whole cluster belongs to the existing cluster
                            cellin = True
                            theKey = key
                            break
            if cellin: # This cluster is already in the Dict
                # Loop over the cells in the new cluster to check that they are actually in the Dict
                for cell in aA:
                    cellin = False
                    for item in self.clusterD[theKey]:
                        if cell[0] == item[0] and cell[1] == item[1]:
                            cellin = True
                    if not cellin:
                        # if the individual cell of this cluster is not registered, add it
                        self.clusterD[theKey] = np.append(self.clusterD[theKey],[cell], axis=0)
                        
            else: 
                # A new cluster is needed, add 1 to nrKeys
                nrKeys += 1
                theKey = nrKeys
                
                # Set Dict cluster to all of the cells in the cluster
                self.clusterD[nrKeys] = aA
                
        # End of for i,node in enumerate(adjA)
        
        # Calculate the number of clusters and the total number of nodes in those clusters
        if not self.CalcClusterContent(totalAdjCells):
       
            if self.params.verbose:
                print ('    Some nodes are registered as duplicates, cleaning clusters')
            
            # If there are duplicates within the clusters, these clusters need to be joined
            
            # Declare variables for holding clusters identied to have duplicates
            keyDelL = []
            keyLinkD ={}
            keyDoneL = []
            duplFound = 0
        
            # Loop over all the clusters
            for key in self.clusterD: 
                
                # Inside the first loop, loop all clusters again
                for keyTwo in self.clusterD:
                    
                    # If the outer and inner loop is for the same, continue
                    if key == keyTwo:
                        continue
                    
                    # if the other loop key is already processed in the inner loop, continue
                    if key in keyDoneL:
                        continue 
                    
                    # Compare all the cells of the other loop with all the cells of the inner loop
                    for cell in self.clusterD[key]:
                        for item in self.clusterD[keyTwo]:
                            #If any cell in the outher and inner llop clusters is the same, note that
                            if cell[0] == item[0] and cell[1] == item[1]:
                                keyDoneL.append(keyTwo)
                                duplFound += 1
                                keyDelL.append(keyTwo)
                                if not key in keyLinkD:
                                    keyLinkD[key] = []
                                if not keyTwo in keyLinkD[key]:
                                    keyLinkD[key].append(keyTwo)
                                
            if self.params.verbose:
                printstr = '    Identified %(d)d duplicate nodes, transferring and deleting' %{'d':duplFound}
                print (printstr)
            
            # Loop over the Dict listing clusters having duplicates
            for key in keyLinkD: 
                # Loop over the duplicated clusters of the outer loop cluster
                for keyTwo in keyLinkD[key]:
                    #loop over the cells in the clusterD[keyTwo], 
                    # if that cell is not in clusterD[key] - move it
                    for cell in self.clusterD[keyTwo]:
                        cellin = False
                        for item in self.clusterD[key]:
                            if cell[0] == item[0] and cell[1] == item[1]:
                                cellin = True
                        if not cellin:
                            self.clusterD[key] = np.append(self.clusterD[key],[cell], axis=0)
                            
                # Delete the redundant cluster
                del self.clusterD[keyTwo]
            
            if not self.CalcClusterContent(totalAdjCells):  
                sys.exit('    Error in the number of clustered cells') 
                
        #return clusterD
        
    def CentralClusterOutlet(self):
        ''' Identify the most geometrically central cell of each cluster
        '''
        
        # Declare Dict to hold the cluster cell closest to the cluster center
        clusterOutletD ={}
        
        for key in self.clusterD:
            avgXY = np.average(self.clusterD[key], axis=0)
            clusterOutletD[key] = self.clusterD[key][cdist([avgXY], self.clusterD[key]).argmin()]
                
        return clusterOutletD
    
    def LargestClusterOutlet(self,featsL):
        ''' Identify the cell with the highest accumulation
        '''
        
        # Declare Dict to hold the cluster cell with he largest upstream area
        clusterOutletD ={}
        
        for key in self.clusterD:
            maxupstream = 0
            for cell in self.clusterD[key]:
                for f in featsL:
                    if f[0] == cell[0] and f[1] == cell[1]:
                        if f[2] > maxupstream:
                            maxupstream = f[2]
                            clusterOutlet = cell
                
     
            clusterOutletD[key] = clusterOutlet
                
        return clusterOutletD
        
    def UniqueOutlets(self, layer, dt, clusterOutlet):
        ''' Separate basin outlets defined as single cells compared to multiple adjacent vells
        '''
        
        # get coordinates of points as 2d array
        coords = np.array([(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY()) for feature in layer])
        
        layer.ResetReading()
        
        # get the feature data as well
        #featsL = [(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY(), feature.GetField("upstream"), feature.GetField("basin")) for feature in layer]
        featsL = [(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY(), feature.GetField("upstream"),
                   feature.GetField("xcoord"), feature.GetField("ycoord"), feature.GetField("basin")) for feature in layer]
        layer.ResetReading()
        
        # Create two lists, holding uniquely separated outlets and those that are adjacent
        adjacentL = [] # list of basin outlets that are at close proximity
        
        uniqueL = [] # list of single cell basin outlets 
        
        # Calculate distance matrix between points
        for i,pt in enumerate(coords):
            
            self.ClosestNode(pt, np.delete(coords, [i], 0))
            
        self.ClusterAdjacent()
            
        if clusterOutlet.lower()[0] == 'c':
            clusterOutletD = self.CentralClusterOutlet()
            if self.params.verbose:
                print ('    Looking for central outlet in clusters')
        else:
            clusterOutletD = self.LargestClusterOutlet(featsL)
            if self.params.verbose:
                print ('    Looking for maximum outlet in clusters')
            
        # Join the clusterMouths to the uniqueL
        for key in clusterOutletD:
            uniqueL.append(clusterOutletD[key])
            
        # Restore the original order of outlets, and add the upstream area
        self.finaloutletL = []
        self.removedOutlets = []
    
        for candpt in featsL:
            candin = False
            for outlet in uniqueL:
                if outlet[0] == candpt[0] and outlet[1] == candpt[1]:
                    self.finaloutletL.append(candpt)
                    candin = True
            if not candin:
                self.removedOutlets.append(candpt)  
                               
    def ClusteredOutlets(self, layer):
        ''' Separate basin outlets defined as single cells compared to multiple adjacent vells
        '''
        
        # get coordinates of points as 2d array
        coords = np.array([(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY()) for feature in layer])
        
        layer.ResetReading()
        
        # Create two lists, holding uniquely separated outlets and those that are adjacent
        self.adjacentL = [] # list of basin outlets that are at close proximity
        
        self.uniqueL = [] # list of single cell basin outlets 
        
        # Calculate distance matrix between points
        for i,pt in enumerate(coords):
            
            self.ClosestNode(pt, np.delete(coords, [i], 0))
            
        clusterD = self.ClusterAdjacent()
        
        return clusterD

    def GetDistilled(self, dstfpn, remdstfpn): 
    
        dstlayer = self.OpenDs(dstfpn)
                          
        # Copy all features to finalOutletL 
        self.finalOutletL = [(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY(), feature.GetField("upstream")) for feature in dstlayer]
        
        self.removedOutletL = False
        
        if remdstfpn:
            remdstlayer = self.OpenDs(remdstfpn)
            self.removedOutletL = [(feature.GetGeometryRef().GetX(), feature.GetGeometryRef().GetY(), feature.GetField("upstream")) for feature in remdstlayer]
        
    def ClosestNode(self, node, nodes):
        ''' Separate nodes based on the closest distance to any other node
        '''
        
        argmindist = cdist([node], nodes).argmin()
        mindist = cdist([node], nodes)[0,argmindist]
        
        if mindist > self.params['threholddist']:
            # The distance to the nearest other point is longer than the given threshold
            self.uniqueL.append(node)
            
        else:
            # The distance to the nearest other point is shorter than or equal to the given threshold
            self.adjacentL.append(node)
    
    def CalcClusterContent(self, totalAdjCells, clusterD, verbose):
        ''' Compare the nr of clusters and then nr of basin starting points
        '''
        
        totalCellsInClusters = 0  
        nClusters= 0    
        
        for key in clusterD:
            nClusters += 1
            totalCellsInClusters += clusterD[key].shape[0]
        if verbose:
            print ('    Total nr of cells in cluster',totalCellsInClusters)
            print ('    Total number of cells to cluster',totalAdjCells)
            
        return totalCellsInClusters == totalAdjCells
           
    def GetSpatialRef(self, layer):
        '''
        '''
        
        self.spatialRef = layer.GetSpatialRef()
        
    def GRASSDEM(self):
        ''' Creates script for hydrological corrected river mouth DEM 
        '''
        
        # Calculate the max distance for the cost grow
        if self.params.clusteroutlet[0] == 'c':
            d = self.widestMouth*0.71
        else:
            d = self.widestMouth
            
        # Set a baseprefix 
        bpL =  self.params.grassdstfn.split('_') 
        
        grassbaseprefex = '%s_%s' %(bpL[1], 'distilled')
            
        print (self.params.grassdstfn)
        
        vin = '%s_%s' %(grassbaseprefex, 'shorewall_pt')
        # import the final outlet points
        cmd = 'v.in.ogr -o input=%(in)s output=%(out)s --overwrite\n' %{'in':self.params.dstfpn,'out':vin}
        print (cmd)
        self.GRASSshF.write(cmd)
        
        # Rasterize the points
        
        vtorast = '%s_%s' %(grassbaseprefex, 'shorewall_pt') 
        cmd = 'v.to.rast input=%(in)s output=%(out)s type=point use=val value=%(d)d --overwrite\n' %{'in':vin, 'out':vtorast, 'd':0 }
        print (cmd)
        
        
        # 
        cmd = 'r.mapcalc "outletlowlevel_DEM_p1 = if(isnull(%(vtorast)s),lowlevel_outlet_costgrow+1,%(vtorast)s)" --overwrite\n' %{'vtorast':vtorast}
        print (cmd)
        
        '''
        # Overlay this on lowleveDEM fixing it at the same time
        DOES NOT WORK AS THE WALL DISAPPEARS AFTER SOSTGROW
                
        cmd = 'r.mapcalc "walledlowlevel_DEM_p1 = if(isnull(%(punch)s),lowlevel_outlet_costgrow+1,%(punch)s)" --overwrite\n' %{'punch':vtorast}
        print (cmd)
        
        vtorast2 = vtorast.replace('_','-')
        
        outfn = '%s.tif' %(self.params.grassbasefn.replace('replaceme',vtorast2))
        outfpn = os.path.join(self.params.dstfp, outfn)
        cmd = 'r.out.gdal -f input=%(in)s format=GTiff type=Int16 output=%(out)s --overwrite' %{'in':vtorast, 'out':outfpn}
        
        # Create a DEM = 1 for the outlet basins
        cmd = 'r.mapcalc "lowlevel_DEM_p1 = lowlevel_outlet_costgrow+1" --overwrite\n'
        self.GRASSshF.write(cmd)
        '''
        
        rcost = '%s_%s' %(grassbaseprefex, 'mouth_dist')     
        # Costgrow from the outlets over the lowlevel_DEM_p1    
        cmd = 'r.cost -n input=outletlowlevel_DEM_p1 output=%(out)s start_points=%(pt)s max_cost=%(d)d --overwrite' %{'pt': self.params.grassdstfn,'d':d, 'out':rcost}
        print (cmd)
        self.GRASSshF.write(cmd)
        
        rcost2 = rcost.replace('_','-')
        outfn = '%s.tif' %(self.params.grassbasefn.replace('replaceme',rcost2))
        outfpn = os.path.join(self.params.dstfp, outfn)
        cmd = 'r.out.gdal -f input=%(in)s format=GTiff type=Int16 output=%(out)s --overwrite' %{'in':rcost, 'out':outfpn}
        print (cmd)

        
        routedem =  '%s_%s' %(grassbaseprefex, 'routing_dem') 
        # Invert lowlevel_DEM_p1 to flow towards the outlet 
        cmd = 'r.mapcalc "%(out)s = int(%(in)s-%(d)d)" --overwrite' %{'out':routedem,'in':rcost,'d':d+2}
        print (cmd)
        self.GRASSshF.write(cmd)
        
        
        routedem2 = routedem.replace('_','-')
        outfn = '%s.tif' %(self.params.grassbasefn.replace('replaceme',routedem2))
        outfpn = os.path.join(self.params.dstfp, outfn)
        cmd = 'r.out.gdal -f input=%(in)s format=GTiff type=Int16 output=%(out)s --overwrite' %{'in':routedem, 'out':outfpn}
        print (cmd)
  
        
        hydrodem = '%s_%s' %(grassbaseprefex, 'basin_dem')
        # Combine the original DEM and the mouth routing DEM
        cmd = 'r.mapcalc "%(hydrodem)s = if(isnull(%(routedem)s),DEM@PERMANENT,%(routedem)s)" --overwrite\n' %{'hydrodem':hydrodem,'routedem':routedem}
        print (cmd)
        self.GRASSshF.write(cmd)
        
        
        
        hydrodem2 = hydrodem.replace('_','-')
        outfn = '%s.tif' %(self.params.grassbasefn.replace('replaceme',hydrodem2))
        outfpn = os.path.join(self.params.dstfp, outfn)
        cmd = 'r.out.gdal -f input=%(hydrodem)s format=GTiff type=Int16 output=%(out)s --overwrite' %{'hydrodem':hydrodem, 'out':outfpn}
        print (cmd)
        self.GRASSshF.write(cmd)
        

        
        if self.params.watershed == 'SFD':
            cmd = 'r.watershed -as elevation=%(hydrodem)s accumulation=%(b)s_acc drainage=%(b)s_draindir threshold=2000'

            
        elif self.params.watershed == 'MFD':
            cmd = 'r.watershed -a elevation=%(hydrodem)s accumulation=%(b)s_acc drainage=%(b)s_draindir threshold=2000'

            
        else:
            convergence = int(self.params.watershed[3])
            
            #r.watershed -a convergence=1 elevation=DEM@PERMANENT accumulation=MFD01_upstream drainage=MFD01_drainage threshold=2000 --overwrite
            
            outprefix = 'no_wall'
            cmd = 'r.watershed -a convergence=%(cnv)d elevation=%(hydrodem)s accumulation=%(b)s_acc drainage=%(b)s_draindir threshold=2000' %{'cnv':convergence, 'hydrodem':hydrodem, 'b':outprefix}
        print (cmd)
        self.GRASSshF.write(cmd)
        
        
        # Mapcalc to crate a byte version for visualisation
        
        cmd = 'r.mapcalc "%(b)s_acc_ln = 10*log(%(b)s_acc)" --overwrite' %{'b':outprefix}
        print (cmd)
        
        cmd = 'r.colors map=%(b)s_acc_ln color=ryb' %{'b':outprefix}
        
        grassname = '%(b)s_acc_ln' %{'b':outprefix}
        print (cmd)
        
        acctif = grassname.replace('_','-')
        outfn = '%s.tif' %(self.params.grassbasefn.replace('replaceme',acctif))
        outfpn = os.path.join(self.params.dstfp, outfn)
        
        cmd = 'r.out.gdal -f input=%(acc)s format=GTiff type=Byte output=%(out)s --overwrite' %{'acc':grassname, 'out':outfpn}
        print (cmd)
        self.GRASSshF.write(cmd)
        
        print (cmd)
        
        BALLE

        
        
        


        SNULLE
        print (self.params.basinRoot)
        SNULE    
        
    def WriteGrassCmd(self):
        '''
        '''
        
        # Set the names of the script and text files
        GRASSshFPN = os.path.join(self.scriptfp, self.GRASSshFN)
        self.GRASSshF = open(GRASSshFPN,'w')
        
        #vFN = os.path.join(datafp, txtFN)
        #vF = open(vFN,'w')
        
        OGRshFPN = os.path.join(self.scriptfp, self.OGRshFN)
        OGRshF = open(OGRshFPN,'w')
        
        '''
        if discharge:
            #import the vector file for the point data
            pvFN = 'discharge_hydro_%(reg)s__%(v)s.shp' %{'reg':region,'v':hydroversion}
            dsn = os.path.join(datafp,pvFN)
        
            printStr = 'v.in.ogr dsn=%(dsn)s output=runoff_stations type=boundary --overwrite\n' %{'dsn':dsn}
            GRASSshF.write(printStr)
            #write the db connections to connect point data to area data
            printStr = 'db.tables driver=dbf database=/Users/thomasg/Documents/grassdata/%(region)s/DEM/dbf\n' %{'region': region}
            GRASSshF.write(printStr)
            printStr = 'db.connect -p\n'
            GRASSshF.write(printStr)
        '''
         
        self.GRASSDEM()
         
        x = 0
        
        print ('    Run these script from your GRASS project command line')
        print ('    chmod u+x %(s)s' %{'s':GRASSshFPN})
        print ('    %(s)s' %{'s':GRASSshFPN})
        
        for outlet in finalOutlets:
            x += 1
            if x < 10:
                xStr= '00000%s' %(x)
            elif x < 100:
                xStr= '0000%s' %(x)
            elif x < 1000:
                xStr= '000%s' %(x)
            elif x < 10000:
                xStr= '00%s' %(x)
            elif x < 100000:
                xStr= '0%s' %(x)
            else:
                xStr= '%s' %(x)
    
            layername = '%(typ)s_%(d)s' %{'typ': basinRoot,'d':xStr}
            
            # Set ESRI shape filename
            shpfn = '%sc.shp' %(layername) 
            shpfpn = os.path.join(datafp, shpfn)
            
            if os.path.exists(shpfpn) and not overwrite:
                continue
            
            if discharge:
                #copy the point data db table to a table to be attached to this vector 
                printStr = 'db.copy from_table=runoff_stations to_table=%(ln)s_table where=id=%(d)d\n' %{'d':row[2], 'ln':layername}
                GRASSshF.write(printStr)
                printStr = 'echo UPDATE %(ln)s_table SET cat=1 WHERE id=%(d)d | db.execute\n' %{'d':row[2], 'ln':layername}
                GRASSshF.write(printStr)
    
            cmdStr = 'r.water.outlet input=%(sm)s_drainage output=%(ln)s coordinates=%(x)f,%(y)f --overwrite\n' %{'sm': paramsD['@outlet'], 'x': outlet[0], 'y':outlet[1], 'ln':layername}
            GRASSshF.write(cmdStr)
    
            cmdStr = '# in 7.8 NULL is set in r.water.outlet # r.null %(ln)s setnull=0\n' %{'ln': layername}
            GRASSshF.write(cmdStr)
            
            cmdStr = 'r.to.vect input=%(ln)s output=%(ln)s type=area --overwrite\n' %{'ln': layername}
            GRASSshF.write(cmdStr)
            
            cmdStr = 'g.remove -f type=raster name=%(ln)s --quiet\n' %{'ln': layername}
            GRASSshF.write(cmdStr)
    
            #clean the basin vector
            cmdStr = 'v.clean input=%(ln)s output=%(ln)sc type=area tool=prune,rmdupl,rmbridge,rmline,rmdangle thresh=0,0,0,0,-1 --overwrite\n' %{'ln': layername}
            GRASSshF.write(cmdStr)
                        
            #Export the cleaned vector to get rid of all non-boundary features 
            #shpfn = '%sc.shp' %(layername) 
            #shpfpn = os.path.join(datafp, shpfn) 
                       
            #cmdStr = 'v.out.ogr input=%(ln)s type=area format=ESRI_Shapefile output=%(shpfpn)s\n\n' %{'ln': layername,'shpfpn':shpfpn}
            #shF.write(cmdStr)
                
            if discharge:
                #connect the fixed table to the vector - the cat should be the same
                printStr = 'v.db.connect -o map=%(ln)sf table=%(ln)s_table\n' %{'ln': layername}
                GRASSshF.write(printStr)
                #add a column to the db
                #printStr = 'v.db.addcol map=%(ln)sf columns="id INT"\n' %{'ln': layername}
                #shF.write(printStr)
                #add the id
                printStr = 'v.to.db map=%(ln)sf layer=1 column=id value=%(d)d\n' %{'ln': layername, 'd':row[2]}
                GRASSshF.write(printStr)        
            
            # Add column for xcoord to vector db
            cmdStr = 'v.db.addcolumn map=%(ln)sc columns="xcoord DOUBLE PRECISION"\n' %{'ln': layername}
            GRASSshF.write(cmdStr)
            
            # add the x coordinate
            cmdStr = 'v.db.update map=%(ln)sc layer=1 column=xcoord value=%(x)f\n' %{'ln': layername, 'x':outlet[0]}
            GRASSshF.write(cmdStr)
            
            # Add column for ycoord to vector db
            cmdStr = 'v.db.addcolumn map=%(ln)sc columns="ycoord DOUBLE PRECISION"\n' %{'ln': layername}
            GRASSshF.write(cmdStr)
            
            # add the y coordinate
            cmdStr = 'v.db.update map=%(ln)sc layer=1 column=ycoord value=%(y)f\n' %{'ln': layername, 'y':outlet[1]}
            GRASSshF.write(cmdStr)        
                            
            # Add the basin square kilometers
            cmdStr = 'v.to.db map=%(ln)sc  type=centroid option=area columns=area_km2 units=kilometers\n\n' %{'ln': layername}
            GRASSshF.write(cmdStr)
            
            #export the vector as an area
            cmdStr = 'v.out.ogr input=%(ln)sc type=area format=ESRI_Shapefile output=%(dsn)s --overwrite\n\n' %{'ln': layername,'dsn':shpfpn}
            GRASSshF.write(cmdStr)
        
            if stats:
                for l in range(7):
                    if l == 0: layer = 'cell_precip_250'; prefix = 'p0'
                    if l == 1: layer = 'cell_runoff_250'; prefix = 'q0' 
                    if l == 2: layer = 'cell_restevap_250'; prefix = 'et0'              
                    if l == 3: layer = 'SFD_upstream_log'; prefix = 'sfd'
                    if l == 4: layer = 'DEM_SRTM_250'; prefix = 'dem'
                    if l == 5: layer = 'DEM_SRTM_slope'; prefix = 'slope'
                    if l == 6: layer = 'DEM_SRTM_pcurve5'; prefix = 'profc'
                    '''                
                    if l == 6: layer = 'cell_runoff_adjusted'; prefix = 'q1'
                    if l == 7: layer = 'cell_restevap_net_log'; prefix = 'et1'
                    if l == 8: layer = 'cell_restevap_adjusted'; prefix = 'et2'
                    if l == 9: layer = 'cell_runoff_final'; prefix = 'q2'
                    if l == 10: layer = 'cell_restevap_remain'; prefix = 'etr'
                    '''
                    printStr = 'v.rast.stats vector=%(ln)sf raster=%(l)s colprefix=%(pf)s\n' %{'ln': layername, 'l':layer, 'pf':prefix} 
                    GRASSshF.write(printStr)
            
                #Export the basin vector with all attributes added       
                dsn = '%(ln)s_attr' %{'ln': layername}
                dsn = os.path.join(datafp,dsn)        
                printStr = 'v.out.ogr input=%(ln)sf type=area dsn=%(dsn)sf\n\n' %{'ln': layername,'dsn':dsn}
                GRASSshF.write(printStr)
                
            '''   
            dstFN = '%(typ)s_%(s)sf.shp' %{'typ': basinRoot,'s': region} 
            dstFPN = os.path.join(datafp,dstFN)
            srcFN = '%(ln)sf.shp' %{'ln': layername}
            srcFPN = os.path.join(dsn,srcFN)
            
            
            if x == 0 and iteration == '01':
                linestr = '/Library/Frameworks/GDAL.framework/Versions/1.9/Programs/ogr2ogr ' 
                linestr = '%(s1)s -skipfailures %(s2)s %(s3)s\n' %{'s1':linestr, 's2':tarFPN, 's3':srcFPN}
            else:
                linestr = '/Library/Frameworks/GDAL.framework/Versions/1.9/Programs/ogr2ogr ' 
                linestr = '%(s1)s -append -skipfailures %(s2)s %(s3)s\n' %{'s1':linestr, 's2':tarFPN, 's3':srcFPN}
            ogrF.write(linestr)
            '''
            #cmdStr = '%(x)f:%(y)f:%(d)d\n' %{'x': row[0], 'y':row[1], 'd':row[2]}
            #vF.write(cmdStr)
    
            if x == 1:
                cmdStr = 'ogr2ogr -skipfailures -nlt MULTIPOLYGON %(dstFPN)s %(srcFPN)s\n' %{'dstFPN': dstPolygonfpn, 'srcFPN':shpfpn}
                OGRshF.write(cmdStr)
            else:
                cmdStr = 'ogr2ogr  -append -skipfailures -nlt MULTIPOLYGON %(dstFPN)s %(srcFPN)s\n' %{'dstFPN': dstPolygonfpn, 'srcFPN':shpfpn}
                OGRshF.write(cmdStr)
                
        GRASSshF.close()
        #vF.close()
        OGRshF.close()
        
def InsidePoly(polyShp,ptL):
    '''
    '''
    pass
    '''
    driver = ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(polyShp,0)
    # Get the only layer within it
    if dataSource is None:
        exitstr = ('Could not open %s' % (polyShp))
        sys.exit(exitstr)

    layer = dataSource.GetLayer(0)
    x = 0
    for feature in layer:
        x += 1
        print (x)
        popL = []
        geom = feature.GetGeometryRef()
        shape = shapely.wkt.loads(geom.ExportToWkt())
        for x,pt in enumerate(ptL):
            if pt.within(shape):
                popL.append(x)
            #pop all points inside polygon
            popL.sort()
            popL.reverse()
            for p in popL:
        
                if len(ptL) == 0:
                    sys.exit('No more mouth points')
                ptL.pop(p)

    '''
    

      
if __name__ == "__main__":
    
    xmlFPN = '/Users/thomasgumbricht/GitHub/geoimagine/script/grass_drainage_outlets.xml'
    
    # Get the parameters from the xml file
    
    params = Params(xmlFPN)
    
    BasinExtract(params)
    
    # dt is the distance threshpöd for signling out a unique bsin outlet
    # dstfpn, dstremfpn, srcfp, SFDsrcfpn, MFDsrcfpn, region, overwrite, paramsD = ReadXML(xmlFPN)
    # stage, overwrite, paramsD, region, dstfpn, dstremfpn, srcfp, SFDsrcfpn, MFDsrcfpn, dsttxtfn, dstPolygonfpn, basinRoot, dsttxtfpn = ReadXML(xmlFPN)
    
    
    
            