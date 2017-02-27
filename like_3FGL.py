import os

# Initialize Fermi and HEASOFT software
print('Initializing Fermi and HEASOFT software.')
os.system('. $HEADAS/headas-init.sh && source $FERMI_DIR/fermi-init.sh')

# Import modules
print 'Importing Python modules.\n'
import sys

import math
from gt_apps import *
try:
    import pyfits
except:
    sys.path.append('/Users/qifeng/soft/PyFITS')
    import pyfits

import pyLikelihood
from decimal import Decimal
from UnbinnedAnalysis import *
from UpperLimits import UpperLimits
from make3FGLxml import *

from optparse import OptionParser
import matplotlib.pyplot as plt
import numpy as np


class analyze3FGL:
    def __init__(self, name='3FGL J2202.7+4217', id_string='BLLac', ra=330.68, dec=42.2778,
                 work_dir='.', d_dir='.', e_dir='events/', emin=100, emax=300000,
                 ext_dir='/Users/qifeng/Data/Blazars/BLLac/Fermi/Extended_archive_v15/Templates/',
                 evlistname='@BLLac_events.list', scfilename='L1701271500537C62B87191_SC00.fits',
                 do_diffuse=True, ev_file_tag='10deg2016', fname_DT='BLLac_like_results_2016.txt',
                 model='BLLac_3FGL_modelPL.xml', modeloutbase='BLLac_3FGL_modelPL_out', rad=10, srcrad=20, TS_lower=5.0):
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
            self.fileDT = open(fname_DT, 'a')
        except IOError:
            print('There was an error opening file {0}'.format(fname_DT))
            sys.exit()

        self.model = model
        self.modeloutbase = modeloutbase
        # dir for extended source model fits files
        self.ext_dir = ext_dir
        # defaults:
        self.zmax = 90
        self.roicut = 'yes'
        self.nenergies = 37
        self.TS_lower = TS_lower

    def get_logistics(self):
        name = self.name
        id_string = self.id_string
        ra = self.ra
        dec = self.dec
        rad = self.rad
        emin = self.emin
        emax = self.emax
        work_dir = self.work_dir
        d_dir = self.d_dir
        e_dir = self.e_dir
        scfilename = self.scfilename
        evlistname = self.evlistname
        do_diffuse = self.do_diffuse
        ev_file_tag = self.ev_file_tag
        fname_DT = self.fname_DT  # file to write results of detection
        model = self.model
        modeloutbase = self.modeloutbase
        return name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase

    def run_gtselect(self, t_start=495936004.000, t_stop=499392004.000):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        if t_stop == None:
            # Recovers TSTOP keyword value from spacecraft.fits file using pyfits and appends value to time.txt
            hdulist = pyfits.open(str(work_dir) + '/' + str(scfilename))
            t_stop = hdulist[1].header['TSTOP']  # TSTOP time from spacecraft file
            print(
            'Setting the stop time to MET value = ' + str(t_stop) + ' from spacecraft file ' + str(scfilename) + '.\n')
        if t_stop == None:
            hdulist = pyfits.open(str(work_dir) + '/' + str(scfilename))
            t_start = hdulist[1].header['TSTART']  # TSTART time from spacecraft file
            print('Setting the start time to MET value = ' + str(t_start) + ' from spacecraft file ' + str(
                scfilename) + '.\n')

        try:
            if not os.path.exists(str(work_dir) + '/' + str(e_dir) + '/' + str(id_string) + '_' + str(
                    ev_file_tag) + 'filtered_events.fits'):
                print('* Running gtselect *\n')
                filter['evclass'] = 128  # new evclass 128 for source and 256 for clean
                filter['evtype'] = 3  # new evtype 3 for front+back, recommended
                filter['ra'] = ra
                filter['dec'] = dec
                filter['rad'] = rad
                filter['emin'] = emin
                filter['emax'] = emax
                filter['zmax'] = self.zmax  # new value, changed from 105 deg
                filter['tmin'] = t_start
                filter['tmax'] = t_stop
                filter['convtype'] = -1
                filter['infile'] = str(evlistname)
                filter['outfile'] = str(work_dir) + '/' + str(e_dir) + '' + str(id_string) + '_' + str(
                    ev_file_tag) + 'filtered_events.fits'
                filter.run()
            if not os.path.exists(str(work_dir) + '/' + str(e_dir) + '/' + str(id_string) + '_' + str(
                    ev_file_tag) + 'filtered_gti_events.fits'):
                print('* Running gtmktime *\n')
                maketime['scfile'] = str(work_dir) + '/' + str(scfilename)
                maketime['filter'] = '(DATA_QUAL>0)&&(LAT_CONFIG==1)'
                maketime['roicut'] = self.roicut
                maketime['evfile'] = str(work_dir) + '/' + str(e_dir) + '/' + str(id_string) + '_' + str(
                    ev_file_tag) + 'filtered_events.fits'
                maketime['outfile'] = str(work_dir) + '/' + str(e_dir) + '/' + str(id_string) + '_' + str(
                    ev_file_tag) + 'filtered_gti_events.fits'
                maketime.run()
        except RuntimeError:
            print('*** error in gtmktime ***\n')

    def run_gtexpcube(self):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        print('* Running gtltcube *\n')
        try:
            if not os.path.exists(str(work_dir) + '/' + str(e_dir) + '/expCube_' + str(ev_file_tag) + '' + str(
                    id_string) + '.fits'):
                expCube['evfile'] = str(work_dir) + '/' + str(e_dir) + '/' + str(id_string) + '_' + str(
                    ev_file_tag) + 'filtered_gti_events.fits'
                expCube['scfile'] = str(work_dir) + '/' + str(scfilename)
                expCube['outfile'] = str(work_dir) + '/' + str(e_dir) + '/expCube_' + str(ev_file_tag) + str(
                    id_string) + '.fits'
                expCube['zmax'] = 90
                expCube['dcostheta'] = 0.025
                expCube['binsz'] = 1
                expCube.run()
        except RuntimeError:
            print('*** error in gtltcube ***\n')

    def run_gtexpmap(self):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        print('* Running gtexpmap *\n')
        if not os.path.exists(str(work_dir) + '/' + str(e_dir) + '/expMap_' + str(ev_file_tag) + '' + str(id_string) + '.fits'):
            expMap['evfile'] = str(work_dir) + '/' + str(e_dir) + '/' + str(id_string) + '_' + str(ev_file_tag) + 'filtered_gti_events.fits'
            expMap['scfile'] = str(work_dir) + '/' + str(scfilename)
            expMap['expcube'] = str(work_dir) + '/' + str(e_dir) + '/expCube_' + str(ev_file_tag) + str(id_string) + '.fits'
            expMap['outfile'] = str(work_dir) + '/' + str(e_dir) + '/expMap_' + str(ev_file_tag) + str(id_string) + '.fits'
            expMap['irfs'] = 'CALDB'  # update this regularly
            expMap['srcrad'] = self.srcrad
            expMap['nlong'] = 120
            expMap['nlat'] = 120
            expMap['nenergies'] = self.nenergies
            expMap.run()

    def make_xml(self):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        print('* Generating source model file *\n')
        if not os.path.exists(str(work_dir) + '/' + str(model)):
            mymodel = srcList('gll_psc_v16.fit', str(work_dir) + '/' + str(e_dir) + '/' + str(id_string) + '_' + str(ev_file_tag) + 'filtered_gti_events.fits', str(work_dir) + '/' + str(model))
            # mymodel.makeModel('gll_iem_v05.fit', 'gll_iem_v05', 'iso_source_v05.txt', 'iso_source_v05')
            mymodel.makeModel('gll_iem_v06.fits', 'gll_iem_v06', 'iso_P8R2_SOURCE_V6_v06.txt', 'iso_P8R2_SOURCE_V6_v06',extDir=self.ext_dir)

    def run_gtdiff(self):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        print('* Running gtdiffrsp *\n')
        diffResps['evfile'] = str(work_dir) + '/' + str(e_dir) + '/' + str(id_string) + '_' + str(ev_file_tag) + 'filtered_gti_events.fits'
        diffResps['scfile'] = str(work_dir) + '/' + str(scfilename)
        diffResps['srcmdl'] = str(work_dir) + '/' + str(model)
        # diffResps['irfs'] = 'P7REP_SOURCE_V15'
        diffResps['irfs'] = 'CALDB'
        diffResps.run()

    def run_gtlike(self, t_start=495936004.000, t_stop=499392004.000):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        print('* Running unbinned likelihood analysis *\n')
        obs = UnbinnedObs(str(work_dir) + '/' + str(e_dir) + '/' + str(id_string) + '_' + str(
            ev_file_tag) + 'filtered_gti_events.fits',
                          str(work_dir) + '/' + str(scfilename),
                          expMap=str(work_dir) + '/' + str(e_dir) + '/expMap_' + str(ev_file_tag) + str(
                              id_string) + '.fits',
                          expCube=str(work_dir) + '/' + str(e_dir) + '/expCube_' + str(ev_file_tag) + str(
                              id_string) + '.fits',
                          irfs='CALDB')
        # DRMNGB is faster for finding initial value
        like = UnbinnedAnalysis(obs, str(work_dir) + '/' + str(model), optimizer='DRMNGB')
        TS_old = -10000
        TS_new = like.Ts(str(name))  # TS value of the source
        n = 0  # TS counter
        # Fitting as long as TS value changes by 0.01
        # likeobj = pyLike.NewMinuit(like.logLike)
        like.tol = 0.1
        # like.fit(verbosity=0)
        likeobj = pyLike.NewMinuit(like.logLike)
        like.fit(verbosity=0, covar=True, optObject=likeobj)
        print("this will be zero if fit converged: {0}".format(likeobj.getRetCode()))
        # pyLike.DRMNGB(like.logLike).getRetCode()
        like.logLike.writeXml(str(work_dir) + '/' + str(modeloutbase) + '_' + str(n) + '.xml')
        flag_converge = 0
        while (abs(TS_old - TS_new) > 0.01 or TS_new is '-inf' or flag_converge == 0):  # Fitting procedure
            try:
                like2 = UnbinnedAnalysis(obs, str(work_dir) + '/' + str(modeloutbase) + '_' + str(n) + '.xml',
                                         optimizer='NewMinuit')
                # like2.fit(verbosity = 0, covar = True)
                like2obj = pyLike.NewMinuit(like2.logLike)
                like2.fit(verbosity=0, covar=True, optObject=like2obj)
                TS_old = TS_new
                TS_new = like.Ts(str(name))
                n += 1
                print(str(n) + '\n')
                print("this will be zero if fit converged: {0}".format(like2obj.getRetCode()))
                if like2obj.getRetCode() > 0:
                    print("***Error: fit did not converge!!!***")
                    sourceDetails = {}
                    for source in like2.sourceNames():
                        sourceDetails[source] = like2.Ts(source)
                    print(sourceDetails)
                    for source, TS in sourceDetails.iteritems():
                        print("source {0}: TS={1}".format(source, TS))
                        if (TS < self.TS_lower):
                            print("Deleting sources with a TS value smaller than {0:.1f}...".format(self.TS_lower))
                            like2.deleteSource(source)
                elif like2obj.getRetCode() == 0:
                    flag_converge = 1
                like2.logLike.writeXml(str(work_dir) + '/' + str(modeloutbase) + '_' + str(n) + '.xml')
                # pyLike.NewMinuit(like2.logLike).getRetCode()
                # Checking if fit converges in 10 rounds
                if (n > 5):
                    print(str(n) + ' iterations in likelihood fitting is over limit.\n')
                    print("*** Quitting the likelihood procedure! ***")
                    break
            except RuntimeError:
                print('*** error in likelihood fit ***\n')
                raise

        if (abs(TS_old - TS_new) < 0.1):  # If the fit has converged in 10 rounds to difference of 0.1
            if (TS_new < 1):  # Upper limit, 1 sigma limit is fine for each bin if the source exists
                UL = UpperLimits(like)
                print('Calculating upper limit.\n')
                try:
                    Limit = UL[str(name)].compute()  # Calculating upper limit
                    MJD = ((t_start + t_stop) / (2 * 86400) + 51910 + 7.428703703703703e-4)  # Converting MET to MJD
                    self.fileDT.write(
                        str("%.2f" % MJD) + '\t' + str(TS_new) + '\t' + str(Limit) + '\t' + str('NA') + '\t' + str(
                            'NA') + '\t' + str('NA') + '\t' + str("NA") + '\n')
                    self.fileDT.close()
                except RuntimeError:
                    print('***  error in unbinned likelihood analysis ***\n')
            else:  # detection, sigma is greater than 1 for the bin considered
                print('Calculating photon and energy fluxes along with uncertainties.\n')
                try:
                    Flux = like2.flux(str(name), emin=filter['emin'],
                                      emax=filter['emax'])  # Checking the energy interval.
                    Ferror = like2.fluxError(str(name), emin=filter['emin'], emax=filter['emax'])
                    Eflux = like2.energyFlux(str(name), emin=filter['emin'], emax=filter['emax'])
                    EfluxError = like2.energyFluxError(str(name), emin=filter['emin'], emax=filter['emax'])
                    print('Writing exposure time.\n')
                    hdulist = pyfits.open(str(work_dir) + '/' + str(e_dir) + '/' + str(id_string) + '_' + str(
                        ev_file_tag) + 'filtered_gti_events.fits')
                    ontime = hdulist[2].header['ONTIME']
                    MJD = ((t_start + t_stop) / (2 * 86400) + 51910 + 7.428703703703703e-4)  # Converting MET to MJD
                    print('* Writing detection to file. *\n')
                    self.fileDT.write(
                        str("%.2f" % MJD) + '\t' + str(TS_new) + '\t' + str(Flux) + '\t' + str(Ferror) + '\t' + str(
                            Eflux) + '\t' + str(EfluxError) + '\t' + str(ontime) + '\n')
                    self.fileDT.close()
                except RuntimeError:
                    print('*** error in unbinned likelihood analysis ***\n')
        self.like = like2
        self.obs = obs

    def plot_like(self, outfile=None, show=False):
        if not hasattr(self, 'like'):
            print("Run run_gtlike first")
        like=self.like
        E = (like.energies[:-1] + like.energies[1:]) / 2.
        plt.figure(figsize=(9, 9))
        #plt.ylim((0.4, 1e4))
        #plt.xlim((200, 300000))
        sum_model = np.zeros_like(like._srcCnts(like.sourceNames()[0]))
        for sourceName in like.sourceNames():
            sum_model = sum_model + like._srcCnts(sourceName)
            if sourceName == self.name:
                ls_ = '-'
            else:
                ls_= '--'
            plt.loglog(E, like._srcCnts(sourceName), label=sourceName[1:], ls=ls_)
        plt.loglog(E, sum_model, label='Total Model')
        plt.errorbar(E, like._Nobs(), yerr=np.sqrt(like._Nobs()), fmt='o', label='Counts')
        #plt.legend()
        lgd = plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
        plt.xlabel('E (MeV)')
        #plt.tight_layout()
        if outfile is not None:
            plt.savefig(outfile, dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
        if show:
            plt.show()

    def run(self, t_start=495936004.000, t_stop=499392004.000, diagnostic_file=None):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        self.run_gtselect(t_start=t_start, t_stop=t_stop)
        self.run_gtexpcube()
        self.run_gtexpmap()
        self.make_xml()
        if do_diffuse:
            self.run_gtdiff()
        self.run_gtlike(t_start=t_start, t_stop=t_stop)
        if diagnostic_file is not None:
            self.plot_like(outfile=diagnostic_file, show=False)
        print('* DONE!!! *')

    def make_TSmap(self, outname='BLLac_3FGL_modelPL_residual.fits',
                   xmlname='BLLac_3FGL_modelPL_freezeIndices_out_11.xml'):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        try:
            print('* Running Tsmap *\n')
            TsMap['statistic'] = "UNBINNED"
            TsMap['scfile'] = str(work_dir) + '/' + str(scfilename)
            TsMap['evfile'] = str(work_dir) + '/' + str(e_dir) + '/' + str(id_string) + '_' + str(
                ev_file_tag) + 'filtered_gti_events.fits'
            TsMap['expcube'] = str(work_dir) + '/' + str(e_dir) + '/expCube_' + str(ev_file_tag) + str(
                id_string) + '.fits'
            TsMap['expmap'] = str(work_dir) + '/' + str(e_dir) + '/expMap_' + str(ev_file_tag) + str(
                id_string) + '.fits'
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

    def make_ctsmap(self, outname="ctsmap.fits"):
        name, id_string, ra, dec, rad, emin, emax, work_dir, d_dir, e_dir, scfilename, evlistname, do_diffuse, ev_file_tag, fname_DT, model, modeloutbase = self.get_logistics()
        try:
            counts_map['evfile'] = str(work_dir) + '/' + str(e_dir) + '/' + str(id_string) + '_' + str(
                ev_file_tag) + 'filtered_gti_events.fits'
            counts_map['outfile'] = str(outname)
            counts_map['scfile'] = str(work_dir) + '/' + str(scfilename)
            counts_map['nxpix'] = 25
            counts_map['nypix'] = 25
            counts_map['binsz'] = 0.25
            counts_map['coordsys'] = "CEL"
            counts_map['xref'] = ra
            counts_map['yref'] = dec
            counts_map['proj'] = 'AIT'
            counts_map.run()
        except RuntimeError:
            print('*** error in making counts map ***\n')

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-n", "--name", dest="name", default='3FGL J2202.7+4217')
    parser.add_option("-s", "--source_name", dest="source_name", default='BLLac')
    parser.add_option("-r", "--ra", dest="ra", default=330.68)
    parser.add_option("-d", "--dec", dest="dec", default=42.2778)
    parser.add_option("--rad", dest="rad", default=10)
    parser.add_option("--srcrad", dest="srcrad", default=20)
    parser.add_option("-l", "--list", dest="evlistname", default='@BLLac_events.list')
    parser.add_option("--scfile", dest="scfilename", default='L1701271500537C62B87191_SC00.fits')
    parser.add_option("-w", "--work_dir", dest="work_dir", default='.')
    parser.add_option("-e", "--e_dir", dest="e_dir", default='events/')
    parser.add_option("--data_dir", dest="d_dir", default='.')
    parser.add_option("--e_min", dest="e_min", default=100)
    parser.add_option("--e_max", dest="e_max", default=300000)
    parser.add_option("--ext_dir", dest="ext_dir",
                      default='/Users/qifeng/Data/Blazars/BLLac/Fermi/Extended_archive_v15/Templates/')
    # parser.add_option("--do_diffuse",dest="do_diffuse",default=True)
    parser.add_option("--skip_diffuse", action="store_false", dest="do_diffuse", default=True, help="don't run gtdiff")
    parser.add_option("-m", "--model", dest="model", default='BLLac_3FGL_modelPL.xml')
    parser.add_option("--modeloutbase", dest="modeloutbase", default='BLLac_3FGL_modelPL_out')
    parser.add_option("-o", "--outfile", dest="fname_DT", default='BLLac_like_results_2016.txt')
    parser.add_option("-t", "--tag", dest="ev_file_tag", default='10deg2016')
    parser.add_option("--t_start", dest="t_start", default=495936004.000)
    parser.add_option("--t_stop", dest="t_stop", default=499392004.000)
    parser.add_option("--do_diagnostic", action="store_true", dest="do_diagnostic", default=False, help="save likelihood fit raw spectrum")
    parser.add_option("--TS_lower", dest="TS_lower", default=3.0)

    (options, args) = parser.parse_args()

    # Defining working directory and source name for analysis
    print('Defining file and source names for LAT analysis.\n')
    bllac3fgl = analyze3FGL(name=options.name, id_string=options.source_name, ra=options.ra, dec=options.dec,
                            work_dir=options.work_dir, d_dir=options.d_dir, e_dir=options.e_dir, emin=options.e_min,
                            emax=options.e_max,
                            ext_dir=options.ext_dir,
                            evlistname=options.evlistname, scfilename=options.scfilename,
                            do_diffuse=options.do_diffuse, ev_file_tag=options.ev_file_tag, fname_DT=options.fname_DT,
                            model=options.model, modeloutbase=options.modeloutbase, rad=options.rad,
                            srcrad=options.srcrad, TS_lower=float(options.TS_lower))
    # t_start = hdulist[1].header['TSTART'] # TSTART time from spacecraft file
    # t_stop = hdulist[1].header['TSTOP'] # TSTOP time from spacecraft file
    # 2011-01-01 00:00:00 315532802.000
    # 2012-01-01 00:00:00 347068802.000
    # 2013-01-01 00:00:00 378691203.000
    # 2014-01-01 00:00:00 410227203.000
    # 2015-01-01 00:00:00 441763203.000
    # t_stop = 410227203.000
    # 381369603.000 is 2013-02-01

    if options.do_diagnostic:
        bllac3fgl.run(t_start=float(options.t_start), t_stop=float(options.t_stop),
                      diagnostic_file="Diagnostic_plot" + str("{0}".format(options.ev_file_tag)) + ".png")
    else:
        bllac3fgl.run(t_start=float(options.t_start), t_stop=float(options.t_stop))
