'''
Created on 31 Dec 2020

@author: thomasgumbricht
'''

import json

from sys import exit

from os import path

#from geoimagine.support import ConvertHVstring, ConvertHVinteger, ConvertXYstring, ConvertXYinteger

#from pprint import pprint

from osgeo.gdal import ColorTable

#from geoimagine.ktpandas import PandasTS

from .layers import VectorLayer, RasterLayer

#import geoimagine.support.karttur_dt as mj_dt
#from palettable import palette


def UpdateDict(mainD, defaultD, jsonFPN = False):
    '''
    '''

    if len(mainD) == 0:

        return defaultD

    elif mainD == 'processid' and jsonFPN:

        exitstr = 'It seems that the process clause of the json file is not a list\n    %s' %(jsonFPN)

        exit(exitstr)

    else:

        d = {key: defaultD.get(key, mainD[key]) for key in mainD}

        for key in defaultD:

            if key not in d:

                mainD[key] = defaultD[key]

        return mainD

class RasterPalette():
    def __init__(self):
        """The constructor is just an empty container.""" 
        
        pass
            
    def SetTuplePalette(self,paletteT):
        self.paletteT = paletteT
        self.FixGDALPalette()
        ct = ColorTable() 

        for c in range(1,len(self.PcR)):
            ct.CreateColorRamp(self.PcR[c-1][0],self.PcR[c-1][1],self.PcR[c][0],self.PcR[c][1])
        self.colortable = ct
       
    def FixGDALPalette(self):
        '''
        '''
        
        if len(self.paletteT) == 0:
            
            return False
        PcR = []
        AT = []
        ATvalue = []
        ATattr = []
        for c in range(len(self.paletteT)):
            cr = (self.paletteT[c][0],(self.paletteT[c][1],self.paletteT[c][2],self.paletteT[c][3]))
            PcR.append(cr)
            ATvalue.append(self.paletteT[c][0])
            ATattr.append(self.paletteT[c][5])
        for a in range(self.paletteT[c][0]): #Up to the highest attibute
            if a in ATvalue:
                x = ATvalue.index(a)
                AT.append(ATattr[x])
            else:
                AT.append('')
        self.PcR = PcR
        self.AT = AT
        self.maxAT = self.paletteT[c][0]
        
class Composition:
    '''
    classdocs
    '''
    def __init__(self, compD, parameters, system, division, defPath):
        '''
        '''
        checkL  = ['source','product','content','layerid','prefix','suffix']

        for key in compD:

            if key in checkL:

                if '_' in compD[key]:

                    exitstr = 'the "%s" parameter can not contain underscore (_): %s ' %(key, compD[key])

                    exit(exitstr)

                if compD[key][0:10] == 'parameter:':

                    param = compD[key].split(':')[1]

                    if not hasattr(parameters,param):

                        exitstr = 'EXITING - error in Composition when identifying default parameter'

                    compD[key] =  getattr(parameters, param)

            setattr(self, key, compD[key])

        if not hasattr(self, 'layerid'):

            exitstr = 'All compositions must contain a layerid'

            exit(exitstr)

        if not hasattr(self, 'content'):

            exitstr = 'All compositions must have a content'

            exit(exitstr)

        if not hasattr(self, 'suffix'):

            self.suffix = '0'

        if self.suffix == '':

            self.suffix = '0'

        self._SetVolume(defPath)

        self._SetExtention(defPath)

        self._SetCompid()

        self._SetSystem(system)

        self._SetDivision(division)

    def _SetVolume(self, defPath):
        ''' Set the volume for each composition
        '''

        if not hasattr(self, 'volume'):

            self.volume = defPath.volume
            
        if not path.exists(path.join('/Volumes',self.volume)):
            
            exitstr = 'Exiting: the volume %s is not attached' %(self.volume)
            
            exit(exitstr)

    def _SetExtention(self, defPath):
        ''' Set the volume for each composition
        '''

        self.ext = defPath.hdr

        if not self.ext[0] == '.':

            self.ext ='.%s' %(self.ext)


        # Redo tbe same procedure for dat file


        self.dat = defPath.hdr

        if len(self.dat) >= 2:

            if self.dat[0] == '.':

                self.dat = self.dat

            else:

                self.dat = '.%s' %(self.dat)

    def _SetSystem(self,system):
        '''
        '''
        self.system = system

    def _SetDivision(self,division):
        '''
        '''
        self.division = division

    def _SetCompid(self):
        '''
        '''
        self.compid = '%(f)s_%(b)s' %{'f':self.content, 'b':self.layerid}

    def _Update(self, compD):
        '''
        '''
        for key in compD:
            if key in self.checkL:
                if '_' in compD[key]:
                    exitstr = 'the "%s" parameter can npot contain underscore (_): %s ' %(key, compD[key])
                    exit(exitstr)
            setattr(self, key, compD[key])
            
    def _SetPalette(self, palettename, session):
        '''
        '''
        
        # The layer "shade" is by default a brightness layer
        if self.layerid in ['shade','hillshade']:
            
            self._SelectPaletteColors('shade', session)
            
        else:
            
            self._SelectPaletteColors(palettename, session)  
 
    def _SelectPaletteColors(self, palettename,session):

        self.palette = RasterPalette()
        
        if palettename == 'default':
            #Look for a default palette for this composition
            query = {'compid':self.compid}
            self.palettename = self.session._SelectCompDefaultPalette(query)
            if self.palette.palettename == None:
                exitstr = 'No default palette for compid %(c)s' %{'c':self.compid}
                exit(exitstr)
        else:
            self.palette.palettename = palettename
            
        # Re create the palette variable
        
        
        
        self.palette.svgPalette = {}
        
        self.palette.svgColor = {}
        
        #Get the palette
        queryD = {'palette':self.palette.palettename}
        
        paramL = ['value','red','green','blue','alpha','label','hint']
               
        recs = session._MultiSearch(queryD, paramL, 'layout', 'rasterpalcolors')
        
        if len(recs) < 2:
            
            exitstr = 'EXITING - palette "%s" not found' %(self.palette.palettename)
            
            print (exitstr)
            
            PALETTEERROR
            
            exit(exitstr)
        
        recs.sort() 
        
      
        for rec in recs:
            
            if 0 <= rec[0] <= 250:
                
                self.palette.svgPalette[rec[0]] = 'rgb(%(r)d,%(g)d,%(b)d)' %{ 'r':rec[1], 'g':rec[2], 'b':rec[3] }
            
            else:
            
                self.palette.svgColor[rec[0]] = 'rgb(%(r)d,%(g)d,%(b)d)' %{ 'r':rec[1], 'g':rec[2], 'b':rec[3] }

        #self.palette.items = recs
        
        self.palette.paletteL = []
        
        print ('recs',recs)
        
        self.palette.SetTuplePalette(recs)
        
        for i in range(256):  
                       
            self.palette.paletteL.append(self.palette.colortable.GetColorEntry(i))
            
class Location:
    ''' Set location
    '''

    def __init__(self, parameters, processid, defregid, procsys, session, srcLocation=True):
        '''
        '''
        
        if srcLocation:
            
            system = procsys.srcsystem; division = procsys.srcdivision; epsg = procsys.srcepsg
            
        else:
            
            system = procsys.dstsystem; division = procsys.dstdivision; epsg = procsys.dstepsg
        
        self.defregid = defregid

        self.system = system

        self.division = division

        self.epsg = epsg

        self.locusD = {}

        self.locusL = []

        if division in ['NA','none','None','na']:

            #No spatial data involved
            pass

        elif division == 'region':
            
            # If defreg id exists as a parameter it is swapped here
            if hasattr(parameters, 'defregid'):
                
                self.defregid = parameters.defregid
                
            self.locusL.append(self.defregid)

            self.locusD[self.defregid] = {'locus':self.defregid, 'path':self.defregid}

        elif division == 'tiles' and system.lower() == 'modis':

            print (session.name)

            paramL = ['htile','vtile']

            #tiles = session._SelectRegionTiles('regionid' ,'modis', self.defregid , paramL)
     
            tiles = session._SelectRegionTiles({'schema':'modis', 'regionid':self.defregid}, paramL)

            for tile in tiles:
                
                hvD = ConvertHVinteger(tile[0],tile[1])

                self.locusL.append(hvD['prstr'])
    
                self.locusD[hvD['prstr']] = {'locus':hvD['prstr'], 'path':hvD}

        elif division == 'region' and system.lower() == 'modis':

            SNULLE
            tiles = session._SelectModisRegionTiles({'regionid':self.defregid})

            for tile in tiles:

                hvD = convTile(tile)

                self.locusL.append(hvD['prstr'])

                self.locusD[hvD['prstr']] = {'locus':hvD['prstr'], 'path':hvD['prstr']}
                
        elif division == 'tiles' and system.lower() == 'ease2n':
            
            paramL = ['xtile','ytile']
            
            tiles = session._SelectRegionTiles({'schema':system, 'regionid':self.defregid}, paramL)            
            
            for tile in tiles:

                xyD = ConvertXYinteger(tile[0],tile[1])

                self.locusL.append(xyD['prstr'])

                self.locusD[xyD['prstr']] = {'locus':xyD['prstr'], 'path':xyD['prstr']}
                
        elif division == 'tiles' and system.lower() == 'export':
            
            if procsys.srcsystem.lower()[0:4] == 'ease':
            
                paramL = ['xtile','ytile']
                
                tiles = session._SelectRegionTiles({'schema':procsys.srcsystem, 'regionid':self.defregid}, paramL)            
                
                for tile in tiles:
    
                    xyD = ConvertXYinteger(tile[0],tile[1])
    
                    self.locusL.append(xyD['prstr'])
    
                    self.locusD[xyD['prstr']] = {'locus':xyD['prstr'], 'path':xyD['prstr']}
            
            else:
                
                print ('export system', division, system, processid)
            
                ERRORCHECK
                    
                                
        elif division == 'tiles' and system.lower() == 'sentinel' and processid[0:7] in ['downloa', 'explode','extract','geochec','findgra','reorgan']:
            
            self.locusL.append('unknown')
        
            self.locusD['unknown'] = {'locus':'unknwon', 'path':'unknown'}
        
        elif division == 'scenes' and system.lower() == 'landsat' and processid[0:7] in ['downloa', 'explode','extract','geochec','findgra','reorgan']:
            
            self.locusL.append('unknown')
            
            self.locusD['unknown'] = {'locus':'unknown', 'path':'unknown'}

        else:
            
            print ('add division, system', division, system, processid)
            
            ERRORCHECK

class TimeSteps:
    """Sets the time span, seasonality and timestep to process data for.
    """

    def __init__(self, period, verbose):
        """The constructor expects a string code for periodicty
        Additionally dates for strat and end, seasonality and addons can be given
        """

        self.period = period

        self.verbose = verbose

        self.datumL = []

        self.datumD = {}

        if not period or period == None:

            # Set to static
            self.SetStaticTimeStep()

        elif self.period.timestep == 'static':

            # Set to static
            self.SetStaticTimeStep()

        elif self.period.timestep == 'singledate':

            self.SingleDateTimeStep()

        elif self.period.timestep == 'singleyear':

            self.SingleYearTimeStep()

        elif self.period.timestep == 'staticmonthly':

            self.SingleStaticMonthlyStep()

        elif self.period.timestep == 'fiveyears':

            self.FiveYearStep()

        # All timestep data containing periods
        else:

            # Get the startdata and enddate of the period
            self.SetStartEndDates()

            self.SetSeasonStartEndDates()

            if self.period.timestep in ['M','MS','monthly','monthlyday']:

                '''
                self.MonthlyTimeStep()

                self.startdatestr = self.startdatestr[0:6]

                self.enddatestr = self.enddatestr[0:6]

                self.pandasCode = 'MS'
                '''

                self.SetMstep()

            elif self.period.timestep == 'varying':

                self.Varying()

            elif self.period.timestep == 'allscenes':

                self.AllScenes()

            elif self.period.timestep == 'inperiod':

                self.InPeriod()

            elif self.period.timestep == 'ignore':

                self.Ignore()

            elif self.period.timestep[len(self.period.timestep)-1] in ['D', '1D']:

                self.SetDstep()

            elif self.period.timestep == '8D':

                self.SetDstep()

            elif self.period.timestep == '16D':

                self.SetDstep()

            else:
                exitstr = 'Unrecognized timestep in class TimeSteps %s' %(self.period.timestep)
                exit(exitstr)

            if self.verbose > 1:

                print ('    Timestep set: %s' %(self.period.timestep))

    def SetStartEndDates(self):
        ''' Set stardata and end date for period
        '''

        self.startdate = mj_dt.IntYYYYMMDDDate(self.period.startyear,
                                self.period.startmonth, self.period.startday)

        self.enddate = mj_dt.IntYYYYMMDDDate(self.period.endyear,
                                self.period.endmonth, self.period.endday)

        self.startdatestr = mj_dt.DateToStrDate(self.startdate)

        self.enddatestr = mj_dt.DateToStrDate(self.enddate)

        if self.enddate < self.startdate:

            exitstr = 'period starts after ending'

            exit(exitstr)

    def SetSeasonStartEndDates(self):
        '''
        '''
        self.startdoy = self.enddoy = 0

        if hasattr(self.period, 'seasonstartmonth') and self.period.seasonstartmonth > 0:

            if hasattr(self.period, 'seasonstartday') and self.period.seasonstartday > 0:

                seasonstart = mj_dt.IntYYYYMMDDDate(2001, self.period.seasonstartmonth, self.period.seasonstartday )

                self.startdoy = int(mj_dt.YYYYDOYStr(seasonstart))

        if hasattr(self.period, 'seasonendmonth') and self.period.seasonendmonth > 0:

            if hasattr(self.period, 'seasonendday') and self.period.seasonendday > 0:

                seasonend = mj_dt.IntYYYYMMDDDate(2001, self.period.seasonendmonth, self.period.seasonendday)

                self.enddoy = int(mj_dt.YYYYDOYStr(seasonend))

        if self.enddoy < self.startdoy:

            errorigen

    def SetStaticTimeStep(self):
        ''' Set to static
        '''
        self.datumL.append('0')

        self.datumD['0'] = {'acqdate':False, 'acqdatestr':'0'}
        
    def SingleDateTimeStep(self):
        
        acqdate = mj_dt.IntYYYYMMDDDate(self.period.startyear,self.period.startmonth,self.period.startday)
        
        acqdatestr = mj_dt.DateToStrDate(acqdate)
        
        self.datumL.append(acqdatestr)

        self.datumD[acqdatestr] = {'acqdate':acqdate, 'acqdatestr':acqdatestr}

    def SingleYearTimeStep(self):
        '''
        '''
        if not self.period.startyear == self.period.endyear:

            exitstr = 'error in period: year'

            exit(exitstr)

        acqdatestr = '%(y)d' %{'y':self.period.startyear}

        if not len(acqdatestr) == 4 or not acqdatestr.isdigit:

            exitstr = 'len(acqdatestr) != 4'

            exit(exitstr)

        self.datumL.append(acqdatestr)

        acqdate = mj_dt.SetYYYY1Jan(int(acqdatestr))

        self.datumD[acqdatestr] = {'acqdate':acqdate, 'acqdatestr':acqdatestr}

    def FiveYearStep(self,periodD):
        if not periodD['startyear'] < periodD['endyear'] or periodD['startyear'] < 1000 or periodD['endyear'] > 9999:
            exitstr = "periodD['startyear'] < periodD['endyear'] or periodD['startyear'] < 1000 or periodD['endyear'] > 9999"
            exit(exitstr)
        for y in range(periodD['startyear'],periodD['endyear']+1,5):
            acqdatestr = '%(y)d' %{'y':y}
            if not len(acqdatestr) == 4:
                exitstr = 'len(acqdatestr) != 4'
                exit(exitstr)



    def SingleStaticMonthlyStep(self,periodD):
        if periodD['endmonth'] < periodD['startmonth'] or periodD['startmonth'] > 12 or periodD['endmonth'] > 12:
            exitstr = "periodD['endmonth'] < periodD['startmonth'] or periodD['startmonth'] > 12 or periodD['endmonth'] > 12"
            exit(exitstr)
        for m in range(periodD['startmonth'],periodD['endmonth']+1):
            if m < 10:
                mstr = '0%(m)d' %{'m':m}
            else:
                mstr = '%(m)d' %{'m':m}
            ERRORCHECK



    def MonthlyTimeStep(self):
        '''
        '''

        mstr = self.MonthToStr( self.period.startmonth )

        yyyymmdd = '%(yyyy)s%(mm)s01' %{'yyyy':self.period.startyear, 'mm':mstr }

        startdate = mj_dt.yyyymmddDate(yyyymmdd)

        #get end date
        mstr = self.MonthToStr(self.period.endmonth)

        yyyymm = '%(yyyy)s%(mm)s' %{'yyyy':self.period.endyear,'mm':mstr }

        enddate = mj_dt.YYYYMMtoYYYYMMDD(yyyymm,32)

        acqdatestr = mj_dt.DateToStrDate(startdate)

        acqdate = startdate

        self.datumL.append(acqdatestr[0:6])

        if len(self.datumL) > 200:

            BALLE

        self.datumD[acqdatestr[0:6]] = {'acqdate':acqdate, 'acqdatestr':acqdatestr[0:6], 'timestep':'M'}

        while True:

            acqdate = mj_dt.AddMonth(acqdate,1)

            if acqdate > enddate:

                break

            acqdatestr = mj_dt.DateToStrDate(acqdate)

            self.datumL.append(acqdatestr[0:6])

            self.datumD[acqdatestr[0:6]] = {'acqdate':acqdate, 'acqdatestr':acqdatestr[0:6], 'timestep':'M'}


    def SetDstep(self):
        '''
        '''

        pdTS = PandasTS(self.period.timestep)

        npTS = pdTS.SetDatesFromPeriod(self.period, self.startdate, self.enddate, pdTS.centralday)

        for d in range(npTS.shape[0]):

            acqdate = npTS[d].date()

            #acqdatestr = mj_dt.DateToStrDate(acqdate)

            acqyyydoy = mj_dt.DateToYYYYDOY(acqdate)

            #acqdatestr = '%(d)d' %{'d':acqyyydoy}

            self.datumL.append(acqyyydoy)

            self.datumD[acqyyydoy] = {'acqdate':acqdate, 'acqdatestr':acqyyydoy}

    def SetMstep(self):
        '''
        '''

        pdTS = PandasTS(self.period.timestep)

        npTS = pdTS.SetMonthsFromPeriod(self)

        for d in range(npTS.shape[0]):

            acqdate = npTS[d].date()

            acqdatestr = mj_dt.DateToStrDate(npTS[d])[0:6]
            
            self.datumL.append(acqdatestr)

            self.datumD[acqdatestr] = {'acqdate':acqdate, 'acqdatestr':acqdatestr}

    def Varying(self):
        self.datumL.append({'acqdatestr':'varying', 'timestep':'varying'})
        ERRORCHECK

    def AllScenes(self):
        '''
        '''
        #self.SetStartEndDates( )

        #self.SetSeasonStartEndDates( )

        #self.datumL.append({'acqdatestr':'allscenes', 'timestep':'allscenes'})
        self.datumL.append('all')
        self.datumD['all'] = {'acqdate':'all', 'acqdatestr':'all', 'startdate':self.startdate, 'enddate':self.enddate, 'startdoy':self.startdoy, 'enddoy':self.enddoy}

        print (self.datumD)

    def Ignore(self):
        self.datumL.append({'acqdatestr':'ignore', 'timestep':'ignore'})
        ERRORCHECK

    def InPeriod(self):
        self.datumL.append({'acqdatestr':'inperiod', 'timestep':'inperiod','startdate':self.startdate, 'enddate':self.enddate})

    def FindVaryingTimestep(self,path):
        ERRORCHECK
        if os.path.exists(path):
            folders = [f for f in os.listdir(path) if os.path.isdir(os.path.join(path, f))]
            self.datumL = []
            for f in folders:
                try:
                    int(f)
                    self.datumL.append({'acqdatestr':f, 'timestep':'varying'})
                except:
                    pass

    def MonthToStr(self,m):
        if m < 10:
            mstr = '0%(m)d' %{'m':m}
        else:
            mstr = '%(m)d' %{'m':m}
        return mstr

    def SetAcqDateDOY(self):
        ERRORCHECK
        for d in self.datumL:
            acqdate = mj_dt.yyyymmddDate(d['acqdatestr'])
            #d['acqdatedaystr'] = mj_dt.DateToYYYYDOY( acqdate)

    def SetAcqDate(self):
        NALLE
        for d in self.datumL:
            pass
            #d['acqdate'] = mj_dt.yyyymmddDate(d['acqdatestr'])

class ProcessParams():
    ''' Class for setting all process parameters
    '''

    def __init__(self, process, pnr, jsonFPN):
        '''
        '''
        
        # Initiate the Basin extract Xtra parameters
        
        #BEXparams.__init__(self)
        
        # Define the process

        self.process = process

        self.pnr = pnr

        self.jsonFPN = jsonFPN

        if self.process.delete and self.process.overwrite:

            strD = {'processid': self.process.processid, 'pnr':self.pnr, 'jsonFPN':self.jsonFPN}

            exitstr = 'EXITING: the subprocess %(processid)s has both delete and overwrite set to true, only one can he true\n' %strD

            exitstr += '    json file with error: %(jsonFPN)s\n' %strD

            exitstr += '    process sequence number in file: %(pnr)d' %strD

            exit(exitstr)

    def _AssembleParameters(self, session):
        ''' Update missing parameters from default setting
        '''

        if not hasattr(self.process,'processid'):
            
            exitstr = 'Process lacking processid - most likely errors in the first dict cluase defining the json object'
            
            exit(exitstr)
            
        if self.process.parameters == None or not hasattr(self.process, 'parameters'):
            
            return

        jsonParams = dict( list( self.process.parameters.__dict__.items() ) )
        
        if not session:
            
            # No db connection, just set parameters and return
            
            self.process.parameters = Struct(jsonParams)
            
            return

        queryD = {'subprocid':self.process.processid, 'parent':'process',
                  'element': 'parameters'}

        paramL =['paramid', 'defaultvalue', 'required', 'paramtyp']
        
        paramRecs = session._MultiSearch(queryD, paramL, 'process','processparams')

        # Create a dict with non-required parameters
        #defaultD  = dict( [ (i[0],i[1]) for i in paramRecs if not i[2] ] )

        defaultL  = [ (i[0],int( i[1] )) for i in paramRecs if not i[2] and i[3].lower()[0:3] == 'int' ]

        defaultL.extend([ (i[0], float( i[1] )) for i in paramRecs if not i[2] and i[3].lower()[0:3] in ['flo','rea'] ] )

        defaultL.extend([ (i[0], i[1]) for i in paramRecs if not i[2] and i[3].lower()[0:3] not in ['int','flo','rea'] ] )

        defaultD = dict (defaultL)
        
        typeD = dict ( [ ( i[0],i[3] ) for i in paramRecs ] )

        # Create a dict with compulsory parameters
        compulsD = dict( [ (i[0],i[1]) for i in paramRecs if i[2] ] )

        # Check that all compulsory parameters are included
        for key in compulsD:

            if not hasattr(self.process.parameters, key):

                exitstr = 'EXITING parameter %s missing for process %s' %(key,self.process.processid)

                exit(exitstr)
                        

        # Update the parameters
        paramD = UpdateDict(jsonParams, defaultD)

        
        # Set the type of all params
        for p in paramD:
            
            try:
                
                if typeD[p].lower()[0:4] == 'bool' and isinstance(paramD[p], str):
                    
                    pass
            except:
                               
                if not p in typeD:
                    
                    exitstr = 'The parameter %s does not exists for the process %s' %(p,self.process.processid)
                    
                    exit(exitstr)
                    
                print ('Errors in typeD', typeD)
                exitstr = 'EXITING - error in identifying boolean variabel for parameter: %s (%s)' %(p, typeD[p])
            
                print (exitstr)
                
            if typeD[p].lower()[0:4] == 'bool' and isinstance(paramD[p], str):
                
                # This is a string that must be converted to boolean

                if paramD[p].lower()[0] == 'f':
                    
                    paramD[p] = False
                    
                else:
                    
                    paramD[p] = True

        self.process.parameters = Struct(paramD)

    def _Verbose(self):
        ''' Set the level of text response
        '''

        if not hasattr(self.process,'verbose'):

            self.process.verbose = 0

    def _SetDefRegion(self, defregion):
        ''' Set default region
        '''

        self.defregion = defregion

    def _SetDb(self, postgresdb):
        ''' Set postgres dataabse as a Struct
        '''

        self.postgresdb = Struct(postgresdb)

    def _SetUserProject(self, userproject):
        ''' Set user project  as a Struct
        '''

        self.userproject = Struct(userproject)

    def _SetLayers(self):
        '''
        '''
        self.dstLayerD = {}
        self.srcLayerD = {}
        self._SetDstLayers()
        self._SetSrcLayers()

    def _SetSrcLayers(self):
        '''
        '''

        self.srcCompL = []
        self.srcLayerExistD = {}; self.srcLayerNonExistD = {};
        self.srcLayerDateExistD = {}; self.srcLayerDateNonExistD = {}
        self.srcLayerCreateD = {}; self.srcLayerDateCreateD = {};

        for comp in self.srcCompD:

            self.srcCompL.append(comp)

        for locus in self.srcLocations.locusL:

            self.srcLayerExistD[locus] = {}; self.srcLayerNonExistD[locus] = {};
            self.srcLayerDateExistD[locus] = {}; self.srcLayerDateNonExistD[locus] = {}
            self.srcLayerCreateD[locus] = {}; self.srcLayerDateCreateD[locus] = {};

            for comp in self.srcCompD:

                self.srcLayerExistD[locus][comp] = []; self.srcLayerNonExistD[locus][comp] = [];
                self.srcLayerDateExistD[locus][comp] = []; self.srcLayerDateNonExistD[locus][comp] = []
                self.srcLayerCreateD[locus][comp] = []; self.srcLayerDateCreateD[locus][comp] = [];

        for locus in self.srcLocations.locusL:

            self.srcLayerD[locus] = {}

            for datum in self.srcPeriod.datumL:

                self.srcLayerD[locus][datum] = {}

                for comp in self.srcCompD:

                    if self.srcCompD[comp].ext.lower() in ['.txt','.csv','.none']:

                        self.srcLayerD[locus][datum][comp] = TextLayer(self.srcCompD[comp], self.srcLocations.locusD[locus], self.srcPeriod.datumD[datum])

                    elif self.srcCompD[comp].ext.lower() in [ '.shp']:

                        self.srcLayerD[locus][datum][comp] = VectorLayer(self.srcCompD[comp], self.srcLocations.locusD[locus], self.srcPeriod.datumD[datum])

                    else:

                        self.srcLayerD[locus][datum][comp] = RasterLayer(self.srcCompD[comp], self.srcLocations.locusD[locus], self.srcPeriod.datumD[datum])

                    if path.exists(self.srcLayerD[locus][datum][comp].FPN):

                        self.srcLayerExistD[locus][comp].append(self.srcLayerD[locus][datum][comp].FPN)
                        self.srcLayerDateExistD[locus][comp].append( datum )
                        #self.session._InsertLayer(self.pp.dstLayerD[self.locus][datum][comp],self.process.overwrite,self.process.delete)

                    else:

                        self.srcLayerNonExistD[locus][comp].append(self.srcLayerD[locus][datum][comp].FPN)
                        self.srcLayerDateNonExistD[locus][comp].append( datum )

                        warnstr = 'WARNING, source file (%(c)s) does not exist: %(p)s' %{'c': comp, 'p':self.srcLayerD[locus][datum][comp].FPN}

                        print (warnstr)


    def _SetDstLayers(self):
        '''
        '''

        self.dstCompL = []
        self.dstLayerExistD = {}; self.dstLayerNonExistD = {};
        self.dstLayerDateExistD = {}; self.dstLayerDateNonExistD = {}
        self.dstLayerCreateD = {}; self.dstLayerDateCreateD = {};

        for comp in self.dstCompD:

            self.dstCompL.append(comp)

        for locus in self.dstLocations.locusL:

            self.dstLayerExistD[locus] = {}; self.dstLayerNonExistD[locus] = {};
            self.dstLayerDateExistD[locus] = {}; self.dstLayerDateNonExistD[locus] = {}
            self.dstLayerCreateD[locus] = {}; self.dstLayerDateCreateD[locus] = {};

            for comp in self.dstCompD:

                self.dstLayerExistD[locus][comp] = []; self.dstLayerNonExistD[locus][comp] = [];
                self.dstLayerDateExistD[locus][comp] = []; self.dstLayerDateNonExistD[locus][comp] = []
                self.dstLayerCreateD[locus][comp] = []; self.dstLayerDateCreateD[locus][comp] = [];

        for locus in self.dstLocations.locusL:

            self.dstLayerD[locus] = {}

            for datum in self.dstPeriod.datumL:

                self.dstLayerD[locus][datum] = {}

                for comp in self.dstCompD:

                    if self.dstCompD[comp].ext == '.shp':

                        self.dstLayerD[locus][datum][comp] = VectorLayer(self.dstCompD[comp], self.dstLocations.locusD[locus], self.dstPeriod.datumD[datum])

                    else:

                        self.dstLayerD[locus][datum][comp] = RasterLayer(self.dstCompD[comp], self.dstLocations.locusD[locus], self.dstPeriod.datumD[datum])
                    
                    if path.exists(self.dstLayerD[locus][datum][comp].FPN):

                        self.dstLayerExistD[locus][comp].append(self.dstLayerD[locus][datum][comp].FPN)
                        self.dstLayerDateExistD[locus][comp].append( datum )

                        if self.process.overwrite:

                            self.dstLayerCreateD[locus][comp].append(self.dstLayerD[locus][datum][comp].FPN)
                            self.dstLayerDateCreateD[locus][comp].append( datum )

                        else:
                            pass
                            #self.session._InsertLayer(self.pp.dstLayerD[self.locus][datum][comp],self.process.overwrite,self.process.delete)


                    else:

                        self.dstLayerNonExistD[locus][comp].append(self.dstLayerD[locus][datum][comp].FPN)
                        self.dstLayerDateNonExistD[locus][comp].append( datum )

                        self.dstLayerCreateD[locus][comp].append(self.dstLayerD[locus][datum][comp].FPN)
                        self.dstLayerDateCreateD[locus][comp].append( datum )

                    '''
                    #self.srcLayerD[locus][datum][comp] = self.srcCompD[comp]
                    #self.dstLayerD[locus][datum][comp]['comp'] = self.dstCompD[comp]
                    if (self.dstLayerD[locus][datum][comp]['comp'].celltype == 'vector'):
                        self.dstLayerD[locus][datum][comp]['layer'] = VectorLayer(self.dstCompD[comp], self.dstLocations.locusD[locus], self.srcPeriod.datumD[datum], self.dstpath)
                    else:
                        self.dstLayerD[locus][datum][comp]['layer'] = RasterLayer(self.dstCompD[comp], self.dstLocations.locusD[locus], self.srcPeriod.datumD[datum], self.dstpath)
                    '''

    def _SetPaths(self, session):
        ''' Set source path and destination path
        '''

        def AssemblePath(srcdstpath, jsonPathD):
            ''' Sub processes for assembling paths by combining db entries and json objects
            '''

            if session:
                                
                # Set query for retrieving path parameters for this process from the db
                queryD = {'subprocid':self.process.processid, 'element': srcdstpath}
    
                # Set the params to retrieve from the db
                paramL = ['paramid', 'paramtyp', 'required', 'defaultvalue']
    
                compRecs = session._MultiSearch(queryD, paramL, 'process','processparams')
    
                # Convert the nested list of paths to lists of parameters and default values
                # only include parameters that are not required by the user
                params = [ i[0] for i in compRecs if not i[2] ]
    
                values = [ i[3] for i in compRecs if not i[2] ]
    
                # Convert the db retrieved parameters and values to a dict
                defaultD = dict( zip( params, values) )
    
                # If this destination path is in the jsonObj, prioritize the given parameters
    
                if jsonPathD:
    
                    jsonPathD = UpdateDict(jsonPathD, defaultD)
    
                else:
    
                    jsonPathD = defaultD

            if srcdstpath == 'srcpath':

                self.srcPath = Struct(jsonPathD)

            else:

                self.dstPath = Struct(jsonPathD)

        ''' Source path'''

        # Get the json object list of source path
        if hasattr(self.process, 'srcpath'):

            jsonPathD = dict( list( self.process.srcpath.__dict__.items() ) )

        else:

            jsonPathD = False

        AssemblePath('srcpath', jsonPathD)


        ''' Destination path'''

        # Get the json object list of destination compositions
        if hasattr(self.process, 'dstpath'):

            jsonPathD = dict( list( self.process.dstpath.__dict__.items() ) )

        else:
            
            jsonPathD = False

        AssemblePath('dstpath', jsonPathD)

    def _SetCompositions(self, session, jsonProcess=False):
        '''
        '''

        # Set the dictionaries to hold the source /src) and destination (dst) compositions
        self.srcCompD = {}; self.dstCompD = {}
        
        if not session:
                
            if 'srccomp' in jsonProcess:
                
                for layerid in jsonProcess['srccomp'][0]:
                    
                    print ('layerid',layerid)
                    
                    self.srcCompD[layerid] = Composition(jsonProcess['srccomp'][0][layerid], self.process.parameters, self.procsys.srcsystem, self.procsys.srcdivision, self.srcPath)

                
            if 'dstcomp' in jsonProcess:
                
                for layerid in jsonProcess['dstcomp'][0]:
                    
                    self.srcCompD[layerid] = Composition(jsonProcess['dstcomp'][0][layerid], self.process.parameters, self.procsys.srcsystem, self.procsys.srcdivision, self.srcPath)

            return

        def GetComps(srcdstcomp):
            '''
            '''
            
            
                
            jsonCompD = {}

            # Select the composition(s) from the database
            queryD = {'subprocid':self.process.processid, 'parent': 'process', 'element': srcdstcomp, 'paramtyp': 'element'}

            paramL = ['paramid', 'defaultvalue', 'required']

            processComps = session._MultiSearch(queryD, paramL, 'process','processparams')

            if len(processComps) > 0:

                if hasattr(self.process, srcdstcomp):

                    if srcdstcomp == 'srccomp':

                        compsL = self.process.srccomp

                    else:

                        compsL = self.process.dstcomp

                    if not isinstance(compsL, list):

                        exitstr = 'Either scrcomp or dstcomp is not a list'

                        exit(exitstr)

                    for jsonComp in compsL:

                        # Convert the jsonComposition from Struct to dict_items
                        dct = jsonComp.__dict__.items()

                        # Loop over the dict_items (only contains 1 item)
                        for item in dct:

                            # convert the dict_item to an ordinary dict
                            jsonCompD[item[0]] = dict ( list(item[1].__dict__.items() ) )

            return processComps, jsonCompD

        def AssembleComp(srcdstcomp, jsonCompD, defaultvalue, layerId):
            ''' Sub processes for assembling compostions by combining db entries and json objects
            '''

            # Set query for retrieving compositions for this process from the db
            queryD = {'subprocid':self.process.processid, 'parent': srcdstcomp, 'element': defaultvalue}

            # Set the params to retrieve from the db
            paramL = ['paramid', 'paramtyp', 'required', 'defaultvalue']

            # Retrieve all compositions for this process from the db
            compParams = session._MultiSearch(queryD, paramL, 'process','processparams')

            # Convert the nested list of compositions to lists of parameters and default values
            # only include parameters that are not required by the user
            params = [ i[0] for i in compParams if not i[2] ]

            values = [ i[3] for i in compParams if not i[2] ]

            # Convert the db retrieved parameters and values to a dict
            defaultD = dict( zip( params, values) )

            # Convert the nested list of compositions to lists of parameters and default values
            # only include parameters that are not required by the user
            params = [ i[0] for i in compParams if i[2] ]

            values = [ i[3] for i in compParams if i[2] ]

            # Convert the user required parameters and values to a dict
            requiredD = dict( zip( params, values) )

            # If this composition is in the jsonObj, prioritize the given json parameters
            if defaultvalue == '*' or defaultvalue == layerId:

                mainD = UpdateDict( jsonCompD, defaultD)

            else:

                mainD = defaultD

            # Check that all required parameters are given
            for key in requiredD:

                if not key in mainD:

                    exitstr = 'EXITING, the required parameter %s in the process %s missing in %s for layer' %(key, self.process.processid, layerId)

                    exit ( exitstr )

            if srcdstcomp == 'srccomp':

                self.srcCompD[layerId] = Composition(mainD, self.process.parameters, self.procsys.srcsystem, self.procsys.srcdivision, self.srcPath)

            else:

                self.dstCompD[layerId] = Composition(mainD, self.process.parameters, self.procsys.dstsystem, self.procsys.dstdivision, self.dstPath)
                
                if hasattr(self.process.parameters,'palette') and self.process.parameters.palette:
                    
                    self.dstCompD[layerId]._SetPalette(self.process.parameters.palette, session)
                
        ''' Source compositions'''

        processComps, jsonCompD = GetComps('srccomp')

        # Start processing all the required compositions

        if len(processComps) > 0:

            # if there is only one comp in the db and with default value == '*',
            if len(processComps) == 1 and processComps[0][1] == '*':

                if len(jsonCompD) > 0:

                    for compkey in jsonCompD:

                        AssembleComp('srccomp', jsonCompD[compkey], '*',compkey)

                else:

                    exitstr = 'Exiting, the process %s need at least one source composition' %(self.pp.processid)

                    exit(exitstr)
            else:

                for rc in processComps:
                    
                    if not rc[1] in jsonCompD:
                        
                        exitstr = 'EXITING - the default compositon %s missing' %(rc[1])
                    
                    #Assemble compositon
                    # Not sure why there are 3 params sent, 2 should be fine 20210329
                    AssembleComp('srccomp', jsonCompD[rc[1]], rc[0] , rc[1])

        ''' Destination compositions'''

        processComps, jsonCompD = GetComps('dstcomp')

        # Loop over all compositions
        if len(processComps) > 0:

            # if there is only one comp in the db and with default value == '*',
            if len(processComps) == 1 and processComps[0][1] == '*':

                if len(jsonCompD) > 0:

                    for compkey in jsonCompD:

                        AssembleComp('dstcomp', jsonCompD[compkey], '*',compkey)

                else:
                    
                    AssembleComp('dstcomp', {}, '*', '*')

            else:

                for rc in processComps:

                    SNULLE

    def _CopyCompositions(self, session, jsonProcess=False):
        ''' For processes where the source comp is copied to the destination comd
        '''
        if not session:
                   
            if 'dstcopy' in jsonProcess:
                
                NOTYETDONE
                for layerid in jsonProcess['srccomp'][0]:
                    
                    self.srcCompD[layerid] = Composition(jsonProcess['dstcomp'][0][layerid], self.process.parameters, self.procsys.srcsystem, self.procsys.srcdivision, self.srcPath)

            return
        
        def AssembleDstComp(jsonCompD, compId, srcD):
            ''' Sub processes for assembling compostions by combining db entries and json objects
            '''

            # Set query for retrieving compositions for this process from the db
            queryD = {'subprocid':self.process.processid, 'parent': 'dstcopy'}

            # Set the params to retrieve from the db
            paramL = ['paramid', 'paramtyp', 'required', 'defaultvalue']

            # Retrieve all compositions for this process from the db
            compParams = session._MultiSearch(queryD, paramL, 'process','processparams')

            # Convert the nested list of compositions to lists of parameters and default values
            # only include parameters that are not required by the user
            params = [ i[0] for i in compParams if not i[2] ]

            values = [ i[3] for i in compParams if not i[2] ]

            # Convert the db retrieved parameters and values to a dict
            defaultD = dict( zip( params, values) )

            # Combine any data given in the json parameter file and the db

            mainD = UpdateDict( jsonCompD, defaultD )
            
            for item in mainD:
                
                if mainD[item] == 'copy':
                    
                    mainD[item] = srcD[item]
                
                if mainD[item] == 'auto':
                    
                    # Reset to original in case there is no other change definition
                    if item in srcD:
                    
                        mainD[item] = srcD[item]
                    
                    paramL = ['prefix','postfix','find','replace']
                    
                    queryD = {'subprocid':self.process.processid, 'compitem':item}
                    
                    if hasattr(self.process.parameters,'mode'):
                        
                        queryD['mode'] = self.process.parameters.mode
                        
                    compLabel = session._SingleSearch(queryD, paramL, 'process','complabeldefaults')
  
                    if compLabel == None:
                        
                        exitstr = 'EXITING, default composition labels missing in table process.complabeldefaults for\n    process: %(subprocid)s, comp: %(compitem)s' %queryD
                    
                        exit (exitstr)
                        
                    if compLabel[0]:
                        
                        mainD[item] = '%s%s' %(compLabel[0], srcD[item])
                        
                    if compLabel[1]:
                        
                        mainD[item] = '%s%s' %(srcD[item],compLabel[1] )
                        
                    if compLabel[2] and compLabel[3]:
                        
                        mainD[item] =  srcD[item].replace(compLabel[2], compLabel[3])
                        
                    elif compLabel[3]:
                        
                        mainD[item] =  compLabel[3]
                                        
            self.dstCompD[compId] = Composition(mainD, self.process.parameters, self.procsys.dstsystem, self.procsys.dstdivision, self.dstPath)
            
            if hasattr(self.process.parameters,'palette') and self.process.parameters.palette:
                    
                self.dstCompD[compId]._SetPalette(self.process.parameters.palette, session)
                    
        # Select the composition(s) from the database
        queryD = {'subprocid':self.process.processid,
                  'parent': 'process', 'element': 'dstcopy', 'paramid': 'srccomp'}

        paramL =['paramid', 'defaultvalue', 'required']

        copyRec = session._MultiSearch(queryD, paramL, 'process','processparams')

        jsonCompD = {}
        
        # srcCompD is a shallow copy of self.src.CompD
        srcCompD = {}

        for c in copyRec:
                
            if not hasattr(self.process,'dstcopy'):
                
                exitstr = 'EXITING - parameters for dstcopy missing for process %s' %(self.process.processid)
                
                exit(exitstr)
                
            # Loop over the destination compositions 
            for item in self.process.dstcopy:
                                            
                compsL = list( item.__dict__.items() )
                    
                    
                print (self.srcCompD)
                
                if len(self.srcCompD) == 1 and '*' in self.srcCompD:
                    
                    # the SrcComp is copied to any number of dstComp
                    # create a shallow copy of 
                    srcCompD[compsL[0][0]] = self.srcCompD['*']
                    
                
                elif not compsL[0][0] in self.srcCompD:
                    
                    exitstr = 'EXITING, can not find source comp for copycomp id = %s' %(compsL[0][0])
                    
                    exit(exitstr) 
                    
                else: 

                    srcCompD[compsL[0][0]] = self.srcCompD[compsL[0][0]]
                    
                # Loop over the destination compositions (copies of srcComp)
                for jsonComp in compsL:

                    # convert each composition to a dict_item (this is not a true dict)
                    dct = jsonComp[1].__dict__.items()

                    # Convert the dict_item to a true dict (oterhwise UpdateDict does not work)
                    jsonCompD[jsonComp[0]] = dict ( list (dct) )

                    # Create a dict from the corresponding source compositon (srcComp)
                    srcD = dict ( list ( srcCompD[jsonComp[0]].__dict__.items() ) )
                    
                    # Combine the json object dict and the srccomp dict
                    mainD = UpdateDict( jsonCompD[jsonComp[0]], srcD)

                    # drop compid
                    mainD.pop('compid')
                    
                    # Complement with any default settings (i.e. the parameters required for a destination composition)
                    AssembleDstComp(mainD,jsonComp[0],srcD)
                

    def _TransferComp(self,session):
        ''' For processes where the destination comp is transfered from the source comp
        '''
        
        if not session:
            
            ### NOT YET IMPLEMENTED
            
            return

        # Select the composition(s) from the database
        queryD = {'subprocid':self.process.processid,
                  'parent': 'process', 'element': 'parameters', 'paramid': 'copycomp'}

        paramL =['defaultvalue', 'required']

        copyRec = session._MultiSearch(queryD, paramL, 'process','processparams')

        if len(copyRec) == 0:

            return

        for c in copyRec:

            if hasattr(self.process.parameters,'copycomp'):

                copyComp = self.process.parameters.copycomp

            elif not c[1]:

                copyComp = c[0]

            if copyComp.lower() == '1to1' and len(self.srcCompD) == 1 and len(self.dstCompD) == 1:

                locus = self.srcLocations.locusL[0]

                comp = list(self.srcCompD.keys())[0]

                datum = self.srcLayerDateExistD[locus][comp][0]

                srcCompD = {comp: self.srcLayerD[locus][datum][comp]._RetrieveLayerComp(session)}

                # Combine the jsonparams and the src composition

                dstCompD = {}

                if not isinstance(self.process.dstcomp, list):

                    exitstr = 'dstcomp is not a list'

                    exit(exitstr)

                for c in self.process.dstcomp:

                    # Convert the jsonComposition from Struct to dict_items
                    dct = c.__dict__.items()

                    # Loop over the dcit_items (only contains 1 item)
                    for item in dct:

                        # convert to ordinary dict
                        dstCompD[item[0]] = dict ( list(item[1].__dict__.items() ) )

                # Combine
                for comp in self.srcCompD:

                    for key in dstCompD[comp]:

                        if dstCompD[comp][key] == 'src':

                            dstCompD[comp][key] = self.srccompD[comp][key]

                # Complete the dstCompD with the db retrieve srcComdD
                self.dstCompD[comp] = UpdateDict(dstCompD[comp], srcCompD[comp])

                # Reset the compostion

                self.dstCompD[comp] = Composition(self.dstCompD[comp], self.process.parameters, self.procsys.dstsystem, self.procsys.dstdivision, self.dstPath)

            else:

                exitstr = 'ADD alternative to paramsjson.ProcessParams._TransferComp: %s' %(copyComp.lower())

        # Reset the dstLayers

        self.dstLayerD = {}
        self._SetDstLayers()

    def _Location(self, session):
        ''' Set source and destination location
        '''

        self.srcLocations = Location(self.process.parameters, self.process.processid, self.defregion, self.procsys, session, True)

        self.dstLocations = Location(self.process.parameters, self.process.processid, self.defregion, self.procsys, session, False)

    def _SetTimeStep(self, timestep):
        ''' Set the timestep for the source and destination of this process
        '''
        if hasattr(self.process, 'srcperiod'):

            self.srcPeriod = TimeSteps(self.process.srcperiod, self.verbose)

        else:

            self.srcPeriod = timestep

        if hasattr(self.process, 'dstperiod'):

            self.dstPeriod = TimeSteps(self.process.dstperiod)

        else:

            self.dstPeriod = timestep

    def _GetProcessSystem(self, session, jsonParams):
        ''' The process system components is predefined
        '''
        
        if session:
            
            queryD = {'subprocid': self.process.processid,'system':jsonParams['userproject']['system']}
    
            paramL = ('srcsystem', 'dstsystem', 'srcdivision', 'dstdivision', 'srcepsg', 'dstepsg')
    
            record = session._SelectProcessSystem(queryD, paramL)
    
            if record == None:
    
                exitstr = 'EXITING - the system userproject system setting %(system)s is not support the process %(subprocid)s' %queryD
    
                exit(exitstr)
                
        else:
            
            record = jsonParams['procsys']
            

        self.procsys = Struct(record)

    def _GetRootProcess(self, session):
        ''' Get the root process to which the sub process belongs
        '''

        queryD = {'subprocid': self.process.processid}

        record = session._SelectRootProcess(queryD)

        if record == None:

            exitstr = 'Exiting - The process "%(subprocid)s" is not registered in the Framework db' %queryD

            exit(exitstr)

        self.rootprocid, self.processStratum = record

    def _AssembleSrcRaw(self, session):
        ''' Set the raw source data
        '''

        srcRawD = {}; self.srcRawD = {}

        for item in self.process.srcraw:

            itemL = list( item.__dict__.items() )
            
            for src in itemL:

                srcRawD[src[0]] = src[1]

        # The srcraw must correspond to the destination compositions
        for key in self.dstCompD:

            if key not in srcRawD:

                exitstr = 'EXITING - the composition %s has no corresponding srcraw element (params.paramsjson._AssembleSrcRaw)' %(key)

                exit(exitstr)

            queryD = {'subprocid':self.process.processid, 'parent':'srcraw',
                      'element': '*'}

            paramL =['paramid', 'defaultvalue', 'required']

            paramRecs = session._AncillaryMultiSearch(queryD, paramL, 'process','processparams')

            # Create a dict with non-required parameters
            defaultD  = dict( [ (i[0],i[1]) for i in paramRecs if not i[2] ] )

            # Create a dict with compulsory parameters
            compulsD = dict( [ (i[0],i[1]) for i in paramRecs if i[2] ] )

            # Check that all compulsory parameters are included
            for c in compulsD:

                if not hasattr(srcRawD[key], c):

                    exitstr = 'EXITING - srcraw element %s missing for composition %s' %(c,key)

                    exit(exitstr)

            # Update the parameters
            jsonParams = dict( list( srcRawD[key].__dict__.items() ) )

            #newDict = UpdateDict(srcRawD[key], defaultD)
            rawDict = UpdateDict(jsonParams, defaultD)

            #self.srcRawD[key] = Struct(srcRawD[key])
            self.srcRawD[key] = Struct(rawDict)

    def _UserStratumRights(self):
        ''' Check if the logged in use has the right to the process and region
        '''
        pass
        '''
        stratum
        print ('hej',self.processStratum)
        print ('hej',self.processStratum)
        '''

class Struct(object):
    ''' Recursive class for building project objects
    '''
    def __init__(self, data):
        for name, value in data.items():
            setattr(self, name, self._wrap(value))

    def _wrap(self, value):
        if isinstance(value, (tuple, list, set, frozenset)):
            return type(value)([self._wrap(v) for v in value])
        else:
            return Struct(value) if isinstance(value, dict) else value

class JsonParams ():
    '''
    classdocs
    '''
    def __init__(self, session=False):
        '''
        '''

        if session:
            
            self.session = session
            
        else:
            
            self.session = False

    def _JsonObj(self,jsonFPN):

        # Set the process list
        processD = {}

        # Read the initial (default) parameters
        defaultjsonFPN = '/Users/thomasgumbricht/Documents/geoimagine_default_thomasg.json'

        # Get the default periodicity
        queryD = {'subprocid': 'Periodicity', 'element':'period'}

        paramL = ['paramid','defaultvalue','paramtyp']
        
        if self.session:
            
            periodL = self.session._MultiSearch(queryD,paramL, 'process','processparams')
            
            pL = [ (i[0],int( i[1]) ) for i in periodL if i[2][0:3] == 'int' ]

            pL.extend([ (i[0],i[1]) for i in periodL if i[2][0:3] != 'int' ] )
            
            periodD = dict(pL)
            
        else:
            
            periodL = ['0']
            
            periodD = {'0': '0'}

            

        

        iniParams = self._JsonParams (defaultjsonFPN)

        self.jsonParams = self._JsonParams (jsonFPN)

        # update userproject, fill in any default setting from the project default json
        self.jsonParams = UpdateDict(self.jsonParams, iniParams)

        # update period, fill in any default database setting
        if 'period' in self.jsonParams:

            self.jsonParams['period'] = UpdateDict(self.jsonParams['period'], periodD)

        else:

            self.jsonParams['period'] = periodD

        # update processes, fill in any default setting from the project default json
        if 'process' in self.jsonParams:

            pnr = -1

            for p in self.jsonParams['process']:

                pnr += 1

                processD[pnr] = {'json': p}

                # This seems to have no effect
                p = UpdateDict(p, iniParams['process'][0], jsonFPN)
                
                if p['dryrun']:
                    
                    p['verbose'] = 3
                
                

        # Convert jasonParams for class attributes
        self.params = Struct(self.jsonParams)

        if self.session:
            
            # Get the default region
            self._GetDefRegion(self.session)

        # Get the overall timestep
        if hasattr(self.params, 'period'):

            timestep = TimeSteps(self.params.period, p['verbose'])

        else:

            # When set to False timestep becomes static
            timestep = TimeSteps( False, p['verbose'] )

        # Loop over the processes defined in this object
        pnr = -1

        for p in self.params.process:

            pnr += 1

            # Set the core information
            PP = ProcessParams(p,pnr,jsonFPN)

            # Check and assemble process parameters
            PP._AssembleParameters(self.session)

            # Said the level of text response
            PP._Verbose()

            if self.session:
                
                # Set default region
                PP._SetDefRegion(self.defregion)
    
                # Set postgres database as Struct
                PP._SetDb(self.jsonParams['postgresdb'])
    
                # Set the user project as a Struct
                PP._SetUserProject(self.jsonParams['userproject'])
    
                # Get the process rootprocess
                PP._GetRootProcess(self.session)
    
                
    
                # Set the system for the postgres session
    
                self.session._SetSystem(self.jsonParams['userproject']['system'])
                
            else:
                
                PP.defregion = self.jsonParams['userproject']['tractid']
                
            # Get the processing system
            PP._GetProcessSystem(self.session, self.jsonParams)

            # Set the locations
            PP._Location(self.session)

            # Get the processing system
            PP._SetTimeStep(timestep)

            # Set the paths
            PP._SetPaths(self.session)

            # Set the compositions
            PP._SetCompositions(self.session, self.jsonParams['process'][0])

            PP._CopyCompositions(self.session, self.jsonParams['process'][0])

            PP._SetLayers()

            PP._TransferComp(self.session)

            PP._UserStratumRights()

            processD[pnr]['PP'] = PP

            processD[pnr]['p'] = p

        return processD

    def _UpdateProjectOld(self, mainD, defaultD):
        '''
        '''

        d = {key: defaultD.get(key, mainD[key]) for key in mainD}

        for key in defaultD:

            if key not in d:

                mainD[key] = defaultD[key]

    def _JsonParams(self,path):
        '''
        '''

        # Opening JSON file
        f = open(path, "r")

        # returns JSON object
        return json.load(f)

    def _GetDict(self):
        '''
        '''
        return self.jsonParams

    def _GetDefRegion(self, session):
        ''' Set default region
        '''

        if self.params.userproject.tractid:

            queryD = {'tract': self.params.userproject.tractid}

            rec = session._SelectTractDefRegion(queryD)

            if rec == None:

                exitstr = 'EXITING - the tractid %(tract)s does not exist in the database' %queryD

                exit(exitstr)

            elif rec[1] == 'D':

                exitstr = 'EXITING - the tractid %(tract)s is set to a default region' %queryD

                exit(exitstr)

            self.defregion = rec[0]

        elif self.self.params.userproject.siteid:

            exit('add site defregion')

        if self.defregion == 'globe':

            exit('EXITING - your defined region "%s" is based on the default region "globe", which is not allowed' %(self.params.userproject.tractid))

        #print ("TGTODO Check if user has the right to this region")

if __name__ == '__main__':

    jsonFPN = '/Users/thomasgumbricht/GitHub/geoimagine-setup_db/doc/xmlsql/general_schema_v80_sql.json'

    params = Params(jsonFPN)
