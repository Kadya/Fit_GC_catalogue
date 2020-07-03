from __future__ import print_function

import pPXF_MUSE.ppxf_MUSE as ppxf_MUSE
import numpy as np
import matplotlib.pyplot as plt
import os
from . import GC_support as sup
import glob
import warnings
warnings.filterwarnings("ignore")


def hists(vs, v, dv):
    bins = np.linspace(v - 3 * dv, v+3*dv, 50)
    fig, ax = plt.subplots()
    ax.hist(vs)


def fit_for_plots(GCs, SNR_threshold=[3, 8], galaxy='FCC47', outdir_plots='./ppxf_plots/'):
    plot_out = outdir_plots + '/{0}/'.format(galaxy)
    for GCi in GCs:
        if GCi.SNR > SNR_threshold[0]:
            print('Getting plot for GC {}'.format(int(GCi.ID)))
            plot_file = '{0}_GC_{1}_ppxf_fit.png'.format(galaxy, GCi.ID)
            ppxf_MUSE.ppxf_MUSE(GCi.spec, GCi.wave, galaxy=galaxy, kin_only=True, lam_range=[
                                4700, 8900], plot_out=plot_out, plot_kin_title=plot_file, save_plots=True)


def fit_single_spec(wave, spec, save=False, filebase='Spec_fit', galaxy='FCC47', plot=True, regul=0, kin_only=False, lam_range=[4700, 8900]):
    plot_pop_title = filebase + '_pop_fit.png'
    pp, miles = ppxf_MUSE.ppxf_MUSE(spec, wave, kin_only=kin_only, regul=regul, save_plots=save,
                                    plot_pop_title=plot_pop_title,
                                    return_pp=True, plot=plot, galaxy=galaxy, lam_range=lam_range)
    return pp, miles


def do_the_fit(ID, output_specs='./Outputs/Specs/FCC47', galaxy='FCC47',  n=300, out_dir='./Outputs/pPXF_Output/FCC47/', kin_only=True, prefix='', alpha_fit=False):
    # print('Doing GC # {0} {1} times'.format(ID, n))
    f = np.loadtxt(output_specs + '{0}_GC_spec_{1}.dat'.format(galaxy, ID), unpack=1)
    wave, spec, spec_bg = f[0], f[1], f[2]
    savebase = '{0}_GC_{1}_{2}'.format(galaxy, ID, prefix)

    print('Doing GC {0}'.format(ID))
    ppxf_MUSE.ppxf_MUSE_MC(spec-spec_bg, wave, galaxy=galaxy, out_dir=out_dir,
                           filebase=savebase, n=n, kin_only=kin_only, cores=7, lam_range=[4700, 8900], alpha_fit=alpha_fit)


def fit_GCs(GCs, SNR_threshold=[3, 8], galaxy='FCC47', n=3, out_dir='./ppxf_output/', fit_new=False, vsys=1444.1):
    alpha_fit = False
    updated_GCs = []
    for GCi in GCs:
        print('Fitting GC {0}'.format(int(GCi.ID)))
        savebase = '{0}_GC_{1}'.format(galaxy, int(GCi.ID))
        print('This GC has a SNR = {0}'.format(np.round(GCi.SNR, 2)))
        if GCi.SNR > SNR_threshold[0]:
            kin_only = True
            if GCi.SNR >= SNR_threshold[1]:
                # print('Doing pop fit for this GC!')
                kin_only = False
        else:
            # print('No fitting for this GC!')
            continue
        if fit_new:
            ppxf_MUSE.ppxf_MUSE_MC(GCi.spec, GCi.wave, galaxy=galaxy, out_dir=out_dir, savetxt=True,
                                   filebase=savebase, n=n, kin_only=kin_only, cores=7, lam_range=[4700, 8900], alpha_fit=alpha_fit, save_plot=True)
        # print(out_dir + savebase + '_'+'{0}_runs*.dat'.format(n))
        file = glob.glob(out_dir + savebase + '_'+'*{0}*runs*.dat'.format(n))[0]
        f = np.loadtxt(file, unpack=1)
        v, dv = sup.get_statistics_from_dist(f[0])

        #bins = np.linspace(v - 3 * dv, v+3*dv, 50)
        #fig, ax = plt.subplots()
        #ax.hist(f[0], bins=bins)
        # plt.show()

        GCi.set_v(v - vsys, dv)
        if not kin_only:
            m, dm = sup.get_statistics_from_dist(f[3])
            print(GCi.ID, m, dm)
        else:
            m, dm = 0, 0
        GCi.set_m(m, dm)
        updated_GCs.append(GCi)
    return updated_GCs


def get_balmer_metallicity(GCs, n=100, galaxy='FCC47', out_dir='./ppxf_output/', vsys=1444.1):
    updated_GCs = []
    for GCi in GCs:
        if GCi.SNR >= 9.8:
            print('SNR = {0}'.format(np.round(GCi.SNR, 2)))
            print('Getting [Fe/H] for GC {0}'.format(GCi.ID))
            savebase = '{0}_GC_{1}_Balmer'.format(galaxy, int(GCi.ID))
            file = glob.glob(out_dir + savebase + '_'+'*{0}*runs*.dat'.format(n))[0]
            f = np.loadtxt(file, unpack=1)
            v, dv = sup.get_statistics_from_dist(f[0])
            print('Without Balmer fit:')
            print('v = {0} +- {1} km/s'.format(np.round(v - vsys, 2), np.round(dv, 2)))
            print('EMILES:')
            print('v = {0} +- {1} km/s'.format(GCi.v, GCi.dv))
            m, dm = sup.get_statistics_from_dist(f[3])
            GCi.m_b = m
            GCi.dm_b = dm
            print('Without Balmer fit:')
            print('m = {0} +- {1} dex'.format(np.round(m, 2), np.round(dm, 2)))
            print('EMILES:')
            print('m = {0} +- {1} dex'.format(GCi.m, GCi.dm))
        updated_GCs.append(GCi)
    return updated_GCs


def get_ages(GCs, n=100, galaxy='FCC47', out_dir='./ppxf_output/', vsys=1444.1):
    updated_GCs = []
    for GCi in GCs:
        if GCi.SNR >= 9.8:
            print('SNR = {0}'.format(np.round(GCi.SNR, 2)))
            print('Getting age for GC {0}'.format(GCi.ID))
            savebase = '{0}_GC_{1}_age'.format(galaxy, int(GCi.ID))
            file = glob.glob(out_dir + savebase + '_'+'*{0}*runs*.dat'.format(n))[0]
            f = np.loadtxt(file, unpack=1)
            v, dv = sup.get_statistics_from_dist(f[0])
            m, dm = sup.get_statistics_from_dist(f[3])
            age, dage = sup.get_statistics_from_dist(f[2])
            GCi.m_age = m
            GCi.dm_age = dm
            GCi.age = age
            GCi.dage = dage
        updated_GCs.append(GCi)
    return updated_GCs


def get_MILES_metallicity(GCs, n=100, galaxy='FCC47', out_dir='./ppxf_output/', vsys=1444.1):
    updated_GCs = []
    for GCi in GCs:
        if GCi.SNR >= 8:
            print('SNR = {0}'.format(np.round(GCi.SNR, 2)))
            print('Getting [Fe/H] for GC {0}'.format(GCi.ID))
            savebase = '{0}_GC_{1}_MILES'.format(galaxy, int(GCi.ID))
            file = glob.glob(out_dir + savebase + '_'+'*{0}*runs*.dat'.format(n))[0]
            f = np.loadtxt(file, unpack=1)
            v, dv = sup.get_statistics_from_dist(f[0])
            print('MILES fit:')
            print('v = {0} +- {1} km/s'.format(np.round(v - vsys, 2), np.round(dv, 2)))
            print('EMILES:')
            print('v = {0} +- {1} km/s'.format(GCi.v, GCi.dv))
            m, dm = sup.get_statistics_from_dist(f[3])
            GCi.m_miles = m
            GCi.dm_miles = dm
            print('MILES fit:')
            print('m = {0} +- {1} dex'.format(np.round(m, 2), np.round(dm, 2)))
            print('EMILES:')
            print('m = {0} +- {1} dex'.format(GCi.m, GCi.dm))
        updated_GCs.append(GCi)
    return updated_GCs


def get_aMILES_metallicity(GCs, n=100, galaxy='FCC47', out_dir='./ppxf_output/', vsys=1444.1):
    updated_GCs = []
    for GCi in GCs:
        if GCi.SNR >= 8:
            print('SNR = {0}'.format(np.round(GCi.SNR, 2)))
            print('Getting [Fe/H] for GC {0}'.format(GCi.ID))
            savebase = '{0}_GC_{1}_aMILES'.format(galaxy, int(GCi.ID))
            file = glob.glob(out_dir + savebase + '_'+'*{0}*runs*.dat'.format(n))[0]
            f = np.loadtxt(file, unpack=1)
            v, dv = sup.get_statistics_from_dist(f[0])
            print('aMILES fit:')
            print('v = {0} +- {1} km/s'.format(np.round(v - vsys, 2), np.round(dv, 2)))
            print('EMILES:')
            print('v = {0} +- {1} km/s'.format(GCi.v, GCi.dv))
            m, dm = sup.get_statistics_from_dist(f[3])
            GCi.m_amiles = m
            GCi.dm_amiles = dm
            print('aMILES fit:')
            print('m = {0} +- {1} dex'.format(np.round(m, 2), np.round(dm, 2)))
            print('EMILES:')
            print('m = {0} +- {1} dex'.format(GCi.m, GCi.dm))
        updated_GCs.append(GCi)
    return updated_GCs


def get_best_SSP_metallicity(GCs, n=100, galaxy='FCC47', out_dir='./ppxf_output/', vsys=1444.1):
    updated_GCs = []
    for GCi in GCs:
        if GCi.SNR >= 10:
            print('SNR = {0}'.format(np.round(GCi.SNR, 2)))
            print('Getting best SSP metallicity for GC {0}'.format(GCi.ID))
            savebase = '{0}_GC_{1}_best'.format(galaxy, int(GCi.ID))
            file = glob.glob(out_dir + savebase + '_'+'*{0}*runs*.dat'.format(n))[0]
            f = np.loadtxt(file, unpack=1)
            v, dv = sup.get_statistics_from_dist(f[0])
            print('best SSP fit:')
            print('v = {0} +- {1} km/s'.format(np.round(v - vsys, 2), np.round(dv, 2)))
            print('EMILES:')
            print('v = {0} +- {1} km/s'.format(GCi.v, GCi.dv))
            m, dm = sup.get_statistics_from_dist(f[3])
            GCi.m_ssp = m
            GCi.dm_ssp = dm
            print('best SSP fit:')
            print('m = {0} +- {1} dex'.format(np.round(m, 2), np.round(dm, 2)))
            print('EMILES:')
            print('m = {0} +- {1} dex'.format(GCi.m, GCi.dm))
        updated_GCs.append(GCi)
    return updated_GCs


def get_CaT_metallicity(GCs, n=100, galaxy='FCC47', out_dir='./ppxf_output/', vsys=1444.1):
    updated_GCs = []
    fig, ax = plt.subplots()
    for GCi in GCs:
        if GCi.SNR >= 8:
            print('SNR = {0}'.format(np.round(GCi.SNR, 2)))
            print('Getting [Fe/H] for GC {0}'.format(GCi.ID))
            savebase = '{0}_GC_{1}_CaT'.format(galaxy, int(GCi.ID))
            file = glob.glob(out_dir + savebase + '_'+'*{0}*runs*.dat'.format(n))[0]
            f = np.loadtxt(file, unpack=1)
            v, dv = sup.get_statistics_from_dist(f[0])
            print('CaT fit:')
            print('v = {0} +- {1} km/s'.format(np.round(v - vsys, 2), np.round(dv, 2)))
            print('EMILES:')
            print('v = {0} +- {1} km/s'.format(GCi.v, GCi.dv))
            m, dm = sup.get_statistics_from_dist(f[3])
            if (np.abs(GCi.v - (v-vsys)) < 30) & (dv < 50):
                GCi.m_CaT = m
                GCi.dm_CaT = dm
                ax.errorbar(GCi.v, (v-vsys), xerr=GCi.dv,
                            yerr=np.sqrt(dv**2), c='k', fmt='o')
            print('CaT fit:')
            print('m = {0} +- {1} dex'.format(np.round(m, 2), np.round(dm, 2)))
            print('EMILES:')
            print('m = {0} +- {1} dex'.format(GCi.m, GCi.dm))

        updated_GCs.append(GCi)
    return updated_GCs


if __name__ == "__main__":  # only executed when the program is run by itself and not imported
    GC_cat_file = '../GC_Finder_MUSE/outputs/FCC161_GC_catalog.fits'
    galaxy = 'FCC161'
    out_dir = './ppxf_output/'
    n = 5
    GCs = sup.initialize_catalog(GC_cat_file)
    kin_only = True
    alpha_fit = False

    for GCi in GCs:
        print('Fitting GC {0}'.format(int(GCi.ID)))
        savebase = '{0}_GC_{1}'.format(galaxy, int(GCi.ID))
        ppxf_MUSE.ppxf_MUSE_MC(GCi.spec, GCi.wave, careful_masking=True, galaxy=galaxy, out_dir=out_dir,
                               filebase=savebase, n=n, kin_only=kin_only, cores=7, lam_range=[4700, 8900], alpha_fit=alpha_fit, save_plot=True)
    # ppxf_MUSE.ppxf_MUSE(spec - spec_bg, wave, plot=True, lam_range = [4700, 8900], galaxy = galaxy, quiet = False)
