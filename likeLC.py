from like_3FGL import *
import numpy as np
from optparse import OptionParser


if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-n","--name",dest="name", default='3FGL J2202.7+4217')
    parser.add_option("-s","--source_name",dest="source_name", default='BLLac')
    parser.add_option("-r","--ra",dest="ra", default=330.68)
    parser.add_option("-d","--dec",dest="dec", default=42.2778)
    parser.add_option("--rad",dest="rad", default=10)
    parser.add_option("--srcrad",dest="srcrad", default=20)
    parser.add_option("-l","--list",dest="evlistname", default='@BLLac_events.list')
    parser.add_option("--scfile",dest="scfilename", default='L1701271500537C62B87191_SC00.fits')
    parser.add_option("-w","--work_dir",dest="work_dir",default='.')
    parser.add_option("-e","--e_dir",dest="e_dir",default='events/')
    parser.add_option("--data_dir",dest="d_dir",default='.')
    parser.add_option("--e_min",dest="e_min",default=100)
    parser.add_option("--e_max",dest="e_max",default=300000)
    parser.add_option("--ext_dir",dest="ext_dir",default='/Users/qifeng/Data/Blazars/BLLac/Fermi/Extended_archive_v15/Templates/')
    #parser.add_option("--do_diffuse",dest="do_diffuse",default=True)
    parser.add_option("--skip_diffuse", action="store_false", dest="do_diffuse", default=True, help="don't run gtdiff")
    parser.add_option("-m", "--model",dest="model",default='BLLac_3FGL_modelPL.xml')
    parser.add_option("--modeloutbase",dest="modeloutbase",default='BLLac_3FGL_modelPL_out')
    parser.add_option("-o", "--outfile",dest="fname_DT",default='BLLac_like_results_2016.txt')
    parser.add_option("-t", "--tag",dest="ev_file_tag",default='10deg2016')
    parser.add_option("--t_start",dest="t_start",default=495936004.000)
    parser.add_option("--t_stop",dest="t_stop",default=499392004.000)
    parser.add_option("--t_bin",dest="t_bin",default=86400.000)
    parser.add_option("--do_diagnostic", action="store_true", dest="do_diagnostic", default=False, help="save likelihood fit raw spectrum")
    parser.add_option("--TS_lower", dest="TS_lower", default=3.0)


    (options, args) = parser.parse_args()

    # Defining working directory and source name for analysis
    print('Defining file and source names for LAT analysis.\n')
    for t in np.arange(float(options.t_start), float(options.t_stop), float(options.t_bin)):

        bllac3fgl = analyze3FGL(name=options.name, id_string=options.source_name, ra=options.ra, dec=options.dec,
                 work_dir=options.work_dir, d_dir=options.d_dir, e_dir=options.e_dir, emin=options.e_min, emax=options.e_max,
                 ext_dir=options.ext_dir,
                 evlistname=options.evlistname, scfilename=options.scfilename,
                 do_diffuse=options.do_diffuse, ev_file_tag=str(t)+options.ev_file_tag, fname_DT=options.fname_DT,
                 model=options.model, modeloutbase=options.modeloutbase, rad=options.rad, srcrad=options.srcrad, TS_lower=float(options.TS_lower))
        if options.do_diagnostic:
            bllac3fgl.run(t_start=t, t_stop=t + options.t_bin, diagnostic_file="Diagnostic_plot"+str("{0:.0f}".format(t))+".png")
        else:
            bllac3fgl.run(t_start=t, t_stop=t+options.t_bin)



