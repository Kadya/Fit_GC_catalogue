from __future__ import print_function
import matplotlib
import pPXF_MUSE.ppxf_MUSE as ppxf_MUSE
import numpy as np
import os
import GC_support as sup
import glob
import warnings
warnings.filterwarnings("ignore")
matplotlib.use('Agg')


def mkdir(path):
    if not os.path.isdir(path):
        os.mkdir(path)


def isfile(path):
    if not os.path.isfile(path):
        print('File {0} not found!'.format(path))
        return 0
    else:
        return 1


def do_gal(galaxy, ID_lim=0):
    GC_cat_file = './catalogs/{0}_GC_cat_to_fit.fits'.format(galaxy)
    out_dir = './ppxf_output/{0}/'.format(galaxy)
    mkdir(out_dir)

    n = 100
    GCs = sup.initialize_catalog(GC_cat_file)
    kin_only = True
    alpha_fit = False

    for GCi in GCs:
        if GCi.ID >= ID_lim:
            print('Fitting GC {0} of {1}'.format(int(GCi.ID), len(GCs)))
            savebase = '{0}_GC_{1}'.format(galaxy, int(GCi.ID))
            if GCi.SNR >= 8:
                kin_only = False
                age_lim = 10
                if GCi.SNR >= 10:
                    age_lim = 1
            else:
                kin_only = True
                age_lim = 12
            if galaxy == 'FCCB1241':
                line_mask = 'line_mask_AO.dat'
            else:
                line_mask = 'line_mask.dat'

            ppxf_MUSE.ppxf_MUSE_MC(GCi.spec, GCi.wave, galaxy=galaxy, out_dir=out_dir, vel=1200,
                                   filebase=savebase, n=n, kin_only=kin_only, cores=6, lam_range=[4700, 8900],
                                   save_plot=True, age_lim=age_lim, mask_file=line_mask, gas_fit=False)


if __name__ == "__main__":  # only executed when the program is run by itself and not imported

    #galaxies = ['FCC202W', 'FCC202E', 'VCC1226', 'VCC579']
    galaxies = ['VCC990', 'VCC2019']
    ID_lims = [3, 0]

    for i, gal in enumerate(galaxies):
        do_gal(gal, ID_lim=ID_lims[i])
