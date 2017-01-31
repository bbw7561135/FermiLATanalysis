import os

# Initialize Fermi and HEASOFT software
print('Initializing Fermi and HEASOFT software.')
os.system('. $HEADAS/headas-init.sh && source $FERMI_DIR/fermi-init.sh')

# Import modules
print 'Importing Python modules.\n'
import sys
sys.path.append('/Users/qifeng/soft/PyFITS')
import math
from gt_apps import *
import pyfits
import pyLikelihood
from decimal import Decimal
from UnbinnedAnalysis import *
from UpperLimits import UpperLimits
from make3FGLxml import *

class analyze3FGL:
    def __init__(self, name='3FGL J2202.7+4217', id_string='BLLac', ra=330.68, dec=42.2778,
                 work_dir='.', d_dir = '.', e_dir = 'events/', emin=100, emax=300000,
                 ext_dir='/Users/qifeng/Data/Blazars/BLLac/Fermi/Extended_archive_v15/Templates/',
                 evlistname = '@BLLac_events.list', scfilename = 'L1701271500537C62B87191_SC00.fits',
                 do_diffuse=True, ev_file_tag = '10deg2016', fname_DT = 'BLLac_like_results_2016.txt',
                 model='BLLac_3FGL_modelPL.xml', modeloutbase = 'BLLac_3FGL_modelPL_out', rad=10, srcrad=20):
        self.name = name
        self.id_string = id_string
        self.ra = ra
        self.dec = dec
        self.rad = rad
        self.srcrad = srcrad
        self.emin = emin
        self.emax = emax
        # Define work and download directories
        self.work_dir = work_dir
        self.d_dir = d_dir
        self.e_dir = e_dir
        self.scfilename = scfilename
        self.evlistname = evlistname
        self.do_diffuse = do_diffuse
        self.ev_file_tag = ev_file_tag
        self.fname_DT = fname_DT  # file to write results of detection
        try:
            # Opening file stream for light curve detection results
            print('Opening file stream for writing results.\n')
            fileDT = open(fname_DT, 'a')
        except IOError:
            print('There was an error opening file {0}'.format(fname_DT))
            sys.exit()

        self.model = model
        self.modeloutbase = modeloutbase
        #dir for extended source model fits files
        self.ext_dir = ext_dir
        #defaults:
        self.zmax = 90
        self.roicut = 'no'
        self.nenergies = 37

    def get_logistics(self):
        name         = self.name
        id_string    = self.id_string
        ra           = self.ra
        dec          = self.dec
        rad          = self.rad
        emin          = self.emin
        emax          = self.emax
        work_dir     = self.work_dir
        d_dir        = self.d_dir
        e_dir        = self.e_dir
        scfilename   = self.scfilename
        evlistname   = self.evlistname
        do_diffuse   = self.do_diffuse
        ev_file_tag  = self.ev_file_tag
        fname_DT     = self.fname_DT  # file to write results of detection
        model        = self.model
        modeloutbase = self.modeloutbase
        return name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase

    def run_gtselect(self, t_start=495936004.000, t_stop=499392004.000):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        if t_stop==None:
            # Recovers TSTOP keyword value from spacecraft.fits file using pyfits and appends value to time.txt
            hdulist = pyfits.open(str(work_dir)+'/'+str(scfilename))
            t_stop = hdulist[1].header['TSTOP'] # TSTOP time from spacecraft file
            print('Setting the stop time to MET value = '+str(t_stop)+' from spacecraft file '+str(scfilename)+'.\n')
        if t_stop==None:
            hdulist = pyfits.open(str(work_dir)+'/'+str(scfilename))
            t_start = hdulist[1].header['TSTART'] # TSTART time from spacecraft file
            print('Setting the start time to MET value = '+str(t_start)+' from spacecraft file '+str(scfilename)+'.\n')

        try:
            if not os.path.exists(str(work_dir)+'/'+str(e_dir)+'/'+str(id_string)+'_'+str(ev_file_tag)+'filtered_events.fits'):
                print('* Running gtselect *\n')
                filter['evclass'] = 128 #new evclass 128 for source and 256 for clean
                filter['evtype'] = 3 #new evtype 3 for front+back, recommended
                filter['ra'] = ra
                filter['dec'] = dec
                filter['rad'] = rad
                filter['emin'] = emin
                filter['emax'] = emax
                filter['zmax'] = self.zmax # new value, changed from 105 deg
                filter['tmin'] = t_start
                filter['tmax'] = t_stop
                filter['convtype'] = -1
                filter['infile'] = str(evlistname)
                filter['outfile'] = str(work_dir)+'/'+str(e_dir)+''+str(id_string)+'_'+str(ev_file_tag)+'filtered_events.fits'
                filter.run()
            if not os.path.exists(str(work_dir)+'/'+str(e_dir)+'/'+str(id_string)+'_'+str(ev_file_tag)+'filtered_gti_events.fits'):
                print('* Running gtmktime *\n')
                maketime['scfile'] = str(work_dir)+'/'+str(scfilename)
                maketime['filter'] = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
                maketime['roicut'] = self.roicut
                maketime['evfile'] = str(work_dir)+'/'+str(e_dir)+'/'+str(id_string)+'_'+str(ev_file_tag)+'filtered_events.fits'
                maketime['outfile'] = str(work_dir)+'/'+str(e_dir)+'/'+str(id_string)+'_'+str(ev_file_tag)+'filtered_gti_events.fits'
                maketime.run()
        except RuntimeError:
            print('*** error in gtmktime ***\n')

    def run_gtexpcube(self):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        print('* Running gtltcube *\n')
        try:
          if not os.path.exists(str(work_dir)+'/'+str(e_dir)+'/expCube_'+str(ev_file_tag)+''+str(id_string)+'.fits'):
            expCube['evfile'] = str(work_dir)+'/'+str(e_dir)+'/'+str(id_string)+'_'+str(ev_file_tag)+'filtered_gti_events.fits'
            expCube['scfile'] = str(work_dir)+'/'+str(scfilename)
            expCube['outfile'] = str(work_dir)+'/'+str(e_dir)+'/expCube_'+str(ev_file_tag)+str(id_string)+'.fits'
            expCube['zmax'] = 90
            expCube['dcostheta'] = 0.025
            expCube['binsz'] = 1
            expCube.run()
        except RuntimeError:
            print('*** error in gtltcube ***\n')

    def run_gtexpmap(self):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        print('* Running gtexpmap *\n')
        if not os.path.exists(str(work_dir)+'/'+str(e_dir)+'/expMap_'+str(ev_file_tag)+''+str(id_string)+'.fits'):
          expMap['evfile'] = str(work_dir)+'/'+str(e_dir)+'/'+str(id_string)+'_'+str(ev_file_tag)+'filtered_gti_events.fits'
          expMap['scfile'] = str(work_dir)+'/'+str(scfilename)
          expMap['expcube'] = str(work_dir)+'/'+str(e_dir)+'/expCube_'+str(ev_file_tag)+str(id_string)+'.fits'
          expMap['outfile'] = str(work_dir)+'/'+str(e_dir)+'/expMap_'+str(ev_file_tag)+str(id_string)+'.fits'
          expMap['irfs'] = 'CALDB'  # update this regularly
          expMap['srcrad'] = self.srcrad
          expMap['nlong'] = 120
          expMap['nlat'] = 120
          expMap['nenergies'] = self.nenergies
          expMap.run()

    def make_xml(self):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        print('* Generating source model file *\n')
        if not os.path.exists(str(work_dir)+'/'+str(model)):
            mymodel = srcList('gll_psc_v16.fit',str(work_dir)+'/'+str(e_dir)+'/'+str(id_string)+'_'+str(ev_file_tag)+'filtered_gti_events.fits',str(work_dir)+'/'+str(model))
            #mymodel.makeModel('gll_iem_v05.fit', 'gll_iem_v05', 'iso_source_v05.txt', 'iso_source_v05')
            mymodel.makeModel('gll_iem_v06.fits', 'gll_iem_v06', 'iso_P8R2_SOURCE_V6_v06.txt', 'iso_P8R2_SOURCE_V6_v06', extDir=self.ext_dir)

    def run_gtdiff(self):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        print('* Running gtdiffrsp *\n')
        diffResps['evfile'] = str(work_dir)+'/'+str(e_dir)+'/'+str(id_string)+'_'+str(ev_file_tag)+'filtered_gti_events.fits'
        diffResps['scfile'] = str(work_dir)+'/'+str(scfilename)
        diffResps['srcmdl'] = str(work_dir)+'/'+str(model)
        #diffResps['irfs'] = 'P7REP_SOURCE_V15'
        diffResps['irfs'] = 'CALDB'
        diffResps.run()

    def run_gtlike(self):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        print('* Running unbinned likelihood analysis *\n')
        obs = UnbinnedObs(str(work_dir)+'/'+str(e_dir)+'/'+str(id_string)+'_'+str(ev_file_tag)+'filtered_gti_events.fits',
                          str(work_dir)+'/'+str(scfilename),
                          expMap=str(work_dir)+'/'+str(e_dir)+'/expMap_'+str(ev_file_tag)+str(id_string)+'.fits',
                          expCube=str(work_dir)+'/'+str(e_dir)+'/expCube_'+str(ev_file_tag)+str(id_string)+'.fits',
                          irfs='CALDB')
        # DRMNGB is faster for finding initial value
        like = UnbinnedAnalysis(obs,str(work_dir)+'/'+str(model),optimizer='DRMNGB')
        TS_old = -10000
        TS_new = like.Ts(str(name)) # TS value of the source
        n = 0 # TS counter
        # Fitting as long as TS value changes by 0.01
        #likeobj = pyLike.NewMinuit(like.logLike)
        like.tol = 0.1
        #like.fit(verbosity=0)
        likeobj = pyLike.NewMinuit(like.logLike)
        like.fit(verbosity=0,covar=True,optObject=likeobj)
        print("this will be zero if fit converged: {0}".format(likeobj.getRetCode()))
        #pyLike.DRMNGB(like.logLike).getRetCode()
        like.logLike.writeXml(str(work_dir)+'/'+str(modeloutbase)+'_'+str(n)+'.xml')
        flag_converge = 0
        while (abs(TS_old - TS_new) > 0.01 or TS_new is '-inf' or flag_converge==0): # Fitting procedure
          try:
            like2 = UnbinnedAnalysis(obs,str(work_dir)+'/'+str(modeloutbase)+'_'+str(n)+'.xml' ,optimizer='NewMinuit')
            #like2.fit(verbosity = 0, covar = True)
            like2obj = pyLike.NewMinuit(like2.logLike)
            like2.fit(verbosity=0,covar=True,optObject=like2obj)
            TS_old = TS_new
            TS_new = like.Ts(str(name))
            n += 1
            print(str(n)+'\n')
            print("this will be zero if fit converged: {0}".format(like2obj.getRetCode()))
            if like2obj.getRetCode() >0:
                 print("***Error: fit did not converge!!!***")
            sourceDetails = {}
            for source in like2.sourceNames():
                sourceDetails[source] = like2.Ts(source)
            print(sourceDetails)
            for source,TS in sourceDetails.iteritems():
                print("source {0}: TS={1}".format(source, TS))
                if (TS < 1.0):
                    #print "Deleting..."
                    like2.deleteSource(source)
                    #print "freezing parameters for weak sources..."
                    #print "not really, deleteing!..."
                    #print like2[source]
                    #pars = like2.freePars(source)
                    #for par in pars:
                    #    like2.setFreeFlag(source, par, 0)
                    #print "now source model", like2[source]
            if like2obj.getRetCode() == 0:
                 flag_converge=1
            like2.logLike.writeXml(str(work_dir)+'/'+str(modeloutbase)+'_'+str(n)+'.xml')
            #pyLike.NewMinuit(like2.logLike).getRetCode()
            # Checking if fit converges in 10 rounds
            if (n > 10):
              print(str(n)+' is over limit.\n')
              break
          except RuntimeError:
              print('*** error in likelihood fit ***\n')

        if (abs(TS_old - TS_new) < 0.1): # If the fit has converged in 10 rounds to difference of 0.1
            if(TS_new < 1): # Upper limit, 1 sigma limit is fine for each bin if the source exists
              UL = UpperLimits(like)
              print('Calculating upper limit.\n')
              try:
               Limit = UL[str(name)].compute() # Calculating upper limit
              except RuntimeError:
                print('***  error in unbinned likelihood analysis ***\n')
        else: # detection, sigma is greater than 1 for the bin considered
            print('Calculating photon and energy fluxes with errors.\n')
            try:
              Flux = like2.flux(str(name), emin = filter['emin'], emax = filter['emax']) # Checking the energy interval.
              Ferror = like2.fluxError(str(name), emin = filter['emin'], emax = filter['emax'])
              Eflux = like2.energyFlux(str(name), emin = filter['emin'], emax = filter['emax'])
              EfluxError = like2.energyFluxError(str(name), emin = filter['emin'], emax = filter['emax'])
              print('Writing exposure time.\n')
              hdulist = pyfits.open(str(work_dir)+'/'+str(e_dir)+'/'+str(id_string)+'_'+str(ev_file_tag)+'filtered_gti_events.fits')
              ontime = hdulist[2].header['ONTIME']
              MJD = ((t_start + t_stop)/(2*86400) + 51910 + 7.428703703703703e-4) # Converting MET to MJD
              print('* Writing detection to file. *\n')
              fileDT.write(str(this_bin)+'\t'+str("%.2f" % MJD)+'\t'+str(TS_new)+'\t'+str(Flux)+'\t'+str(Ferror)+'\t'+str(Eflux)+'\t'+str(EfluxError)+'\t'+str(ontime)+'\n')
              fileDT.close()
            except RuntimeError:
              print('*** error in unbinned likelihood analysis ***\n')

    def run(self, t_start=495936004.000, t_stop=499392004.000):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        self.run_gtselect(t_start= t_start, t_stop=t_stop)
        self.run_gtexpcube()
        self.run_gtexpmap()
        self.make_xml()
        if do_diffuse:
            self.run_gtdiff()
        self.run_gtlike()
        print('* DONE!!! *')

    def make_TSmap(self, outname='J0648_3FGL_modelPL_residual.fits',
                   xmlname = 'J0648_3FGL_modelPL_fix_out_2.xml'):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        try:
            print('* Running Tsmap *\n')
            TsMap['statistic'] = "UNBINNED"
            TsMap['scfile'] = str(work_dir)+'/'+str(scfilename)
            TsMap['evfile'] = str(work_dir)+'/'+str(e_dir)+'/'+str(id_string)+'_'+str(ev_file_tag)+'filtered_gti_events.fits'
            TsMap['expcube'] = str(work_dir)+'/'+str(e_dir)+'/expCube_'+str(ev_file_tag)+str(id_string)+'.fits'
            TsMap['expmap'] = str(work_dir)+'/'+str(e_dir)+'/expMap_'+str(ev_file_tag)+str(id_string)+'.fits'
            TsMap['srcmdl'] = str(xmlname)
            TsMap['irfs'] = 'CALDB'  # update this regularly
            TsMap['optimizer'] = "NEWMINUIT"
            TsMap['outfile'] = str(outname)
            TsMap['nxpix'] = 25
            TsMap['nypix'] = 25
            TsMap['binsz'] = 0.5
            TsMap['coordsys'] = "CEL"
            TsMap['xref'] = ra
            TsMap['yref'] = dec
            TsMap['proj'] = 'STG'
            TsMap.run()
        except RuntimeError:
            print('*** error in making TS map ***\n')


if __name__ == '__main__':
    # Defining working directory and source name for analysis
    print('Defining file and source names for LAT analysis.\n')
    bllac3fgl = analyze3FGL()
    #t_start = hdulist[1].header['TSTART'] # TSTART time from spacecraft file
    #t_stop = hdulist[1].header['TSTOP'] # TSTOP time from spacecraft file
    #2011-01-01 00:00:00 315532802.000
    #2012-01-01 00:00:00 347068802.000
    #2013-01-01 00:00:00 378691203.000
    #2014-01-01 00:00:00 410227203.000
    #2015-01-01 00:00:00 441763203.000

    t_start = 495936004.000
    #t_stop = 410227203.000
    #381369603.000 is 2013-02-01
    t_stop = 499392004.000

    bllac3fgl.run(t_start=t_start, t_stop=t_stop)
