
from __future__ import print_function
import matplotlib.pyplot as plt
import os
from matplotlib.colors import LogNorm
import numpy as np
from copy import deepcopy as copy
from filter_routines import process_spec
from astropy.io import fits
import warnings
warnings.filterwarnings("ignore")


class GC():
    def __init__(self, ID, x, y, ra, dec, SNR=0, g=np.nan, z=np.nan, wave=None, spec=None, bg_spec=None, v=0, dv=0, m=0, dm=0, r=0, galaxy='FCC47', g_MUSE=0,
                 z_MUSE=0, age=0, dage=0, m_age=0, dm_age=0, m_miles=0, dm_miles=0, m_amiles=0, dm_amiles=0, m_CaT=0, dm_CaT=0, m_ssp=0, dm_ssp=0, dm_b=0, m_b=0):
        self.ID = ID
        self.x = x
        self.y = y
        self.ra = ra
        self.dec = dec
        self.SNR = SNR
        self.g = g
        self.z = z
        self.r = r
        self.wave = wave
        spec[np.isinf(spec)] = 0.0
        spec[np.isnan(spec)] = 0.0
        self.spec = spec
        bg_spec[np.isinf(bg_spec)] = 0.0
        bg_spec[np.isnan(bg_spec)] = 0.0
        self.bg_spec = bg_spec
        self.v = v
        self.dv = dv
        self.m = m
        self.dm = dm
        self.g_MUSE = g_MUSE  # process_spec(self.wave, self.spec, 'F475W')
        self.z_MUSE = z_MUSE  # process_spec(self.wave, self.spec, 'F850LP')
        self.galaxy = galaxy
        self.age = age
        self.dage = dage
        self.m_age = m_age
        self.dm_age = dm_age
        self.m_miles = m_miles
        self.dm_miles = dm_miles
        self.m_amiles = m_amiles
        self.dm_amiles = dm_amiles
        self.m_CaT = m_CaT
        self.dm_CaT = dm_CaT
        self.m_b = m_b
        self.dm_b = dm_b
        self.m_ssp = m_CaT
        self.dm_ssp = dm_CaT

    def set_v(self, v=0, dv=0):
        self.v = v
        self.dv = dv

    def set_ra_dec(self, ra=0, dec=0):
        self.ra = np.round(ra, 7)
        self.dec = np.round(dec, 7)

    def set_m(self, m=0, dm=0):
        self.m = m
        self.dm = dm

    def set_ID(self, ID=0):
        self.ID = int(ID)

    def get_MUSE_mag(self):
        self.g_MUSE = process_spec(self.wave, self.spec, 'F475W')
        self.z_MUSE = process_spec(self.wave, self.spec, 'F850LP')

    def set_r(self, r=0):
        self.r = r


def get_rad(x_pix, y_pix, xc=207, yc=387):
    '''
    Simple routine to calculate the projected distance in arcsec
    '''
    rs = np.zeros(len(x_pix))
    for i in range(len(x_pix)):
        rs[i] = np.sqrt((x_pix[i] - xc)**2 + (y_pix[i] - yc)**2)  # in pixel coords
    return rs * 0.2  # in arcsec


def create_catalog(GCs, output_directory='./', galaxy='FCC47'):
    '''
    Creates the MUSE catalog as a fits file in the output directory
    '''
    # the spectra
    wave = GCs[0].wave

    filename = '{0}_GC_catalog.fits'.format(galaxy)
    print("- Writing: " + output_directory + filename)
    cols = []
    cols.append(fits.Column(name='GC_ID',  format='D', array=[GCs[i].ID for i in range(len(GCs))]))
    cols.append(fits.Column(name='X_PIX',    format='D', array=[GCs[i].x for i in range(len(GCs))]))
    cols.append(fits.Column(name='Y_PIX',    format='D', array=[GCs[i].y for i in range(len(GCs))]))
    cols.append(fits.Column(name='SNR',    format='D', array=[GCs[i].SNR for i in range(len(GCs))]))
    cols.append(fits.Column(name='RA',    format='D', array=[GCs[i].ra for i in range(len(GCs))]))
    cols.append(fits.Column(name='DEC',    format='D', array=[GCs[i].dec for i in range(len(GCs))]))
    cols.append(fits.Column(name='g',    format='D', array=[GCs[i].g for i in range(len(GCs))]))
    cols.append(fits.Column(name='z',    format='D', array=[GCs[i].z for i in range(len(GCs))]))
    cols.append(fits.Column(name='v',    format='D', array=[GCs[i].v for i in range(len(GCs))]))
    cols.append(fits.Column(name='dv',    format='D', array=[GCs[i].dv for i in range(len(GCs))]))
    cols.append(fits.Column(name='m',    format='D', array=[GCs[i].m for i in range(len(GCs))]))
    cols.append(fits.Column(name='dm',    format='D', array=[GCs[i].dm for i in range(len(GCs))]))
    cols.append(fits.Column(name='r', format='D', array=[GCs[i].r for i in range(len(GCs))]))
    cols.append(fits.Column(name='age', format='D', array=[GCs[i].age for i in range(len(GCs))]))
    cols.append(fits.Column(name='g_MUSE',    format='D',
                            array=[GCs[i].g_MUSE for i in range(len(GCs))]))
    cols.append(fits.Column(name='z_MUSE',    format='D',
                            array=[GCs[i].z_MUSE for i in range(len(GCs))]))
    cols.append(fits.Column(name='m_miles',    format='D',
                            array=[GCs[i].m_miles for i in range(len(GCs))]))
    cols.append(fits.Column(name='dm_miles',    format='D',
                            array=[GCs[i].dm_miles for i in range(len(GCs))]))
    cols.append(fits.Column(name='m_amiles',    format='D',
                            array=[GCs[i].m_amiles for i in range(len(GCs))]))
    cols.append(fits.Column(name='dm_amiles',    format='D',
                            array=[GCs[i].dm_amiles for i in range(len(GCs))]))
    cols.append(fits.Column(name='m_CaT',    format='D',
                            array=[GCs[i].m_CaT for i in range(len(GCs))]))
    cols.append(fits.Column(name='dm_CaT',    format='D',
                            array=[GCs[i].dm_CaT for i in range(len(GCs))]))
    cols.append(fits.Column(name='m_b',    format='D',
                            array=[GCs[i].m_b for i in range(len(GCs))]))
    cols.append(fits.Column(name='dm_b',    format='D',
                            array=[GCs[i].dm_b for i in range(len(GCs))]))
    cols.append(fits.Column(name='m_ssp',    format='D',
                            array=[GCs[i].m_ssp for i in range(len(GCs))]))
    cols.append(fits.Column(name='dm_ssp',    format='D',
                            array=[GCs[i].dm_ssp for i in range(len(GCs))]))
    cols.append(fits.Column(name='galaxy',    format='A6',
                            array=[GCs[i].galaxy for i in range(len(GCs))]))
    tbhdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols))

    dwave = wave[1] - wave[0]
    wave0 = wave[0]
    wave1 = wave[-1]
    hdr = fits.Header()
    hdr['wave0'] = np.round(wave0, 4)
    hdr['wave1'] = np.round(wave1 + dwave, 4)
    hdr['dwave'] = dwave
    hdr['galaxy'] = galaxy

    hdu = fits.PrimaryHDU(header=hdr)
    hdu1 = fits.ImageHDU([GCs[i].spec + GCs[i].bg_spec for i in range(len(GCs))])
    hdu2 = fits.ImageHDU([GCs[i].bg_spec for i in range(len(GCs))])
    hdulist = fits.HDUList([hdu, tbhdu, hdu1, hdu2])
    if os.path.isfile(output_directory + filename):
        os.remove(output_directory + filename)
    hdulist.writeto(output_directory + filename)


def read_fits(file, i=1):
    hdulist = fits.open(file)
    data = hdulist[i].data
    # print(hdulist[i].header)
    return data, hdulist[i].header


def initialize_catalog(file, already_fitted=False, galaxy='FCC47'):
    GC_cat_table = read_fits(file, i=1)[0]
    GC_specs_unsubtracted = read_fits(file, i=2)[0]
    GC_specs_bg = read_fits(file, i=3)[0]
    header = read_fits(file, i=0)[1]
    wave = np.arange(header['wave0'], np.round(header['wave1'], 4), header['dwave'])
    if len(wave) == len(GC_specs_unsubtracted[0]) + 1:
        wave = wave[:-1]
    GCs = []
    miles = False
    amiles = False
    CaT = False
    get_ages = False
    best_ssp = False
    # print(GC_cat_table.names)
    if 'm_miles' in GC_cat_table.names:
        miles = True
    if 'm_amiles' in GC_cat_table.names:
        amiles = True
    if 'm_CaT' in GC_cat_table.names:
        CaT = True
    if 'm_ssp' in GC_cat_table.names:
        best_ssp = True
    if 'm_b' in GC_cat_table.names:
        balmer_met = True
    else:
        balmer_met = False
    if 'm_age' in GC_cat_table.names:
        get_ages = True
    for i, ID in enumerate(GC_cat_table.field('GC_ID')):
        xpix, ypix = GC_cat_table.field('X_PIX')[i], GC_cat_table.field('Y_PIX')[i]
        ra, dec = GC_cat_table.field('RA')[i], GC_cat_table.field('DEC')[i]
        SNR = GC_cat_table.field('SNR')[i]
        g = GC_cat_table.field('g')[i]
        z = GC_cat_table.field('z')[i]
        spec = GC_specs_unsubtracted[i] - GC_specs_bg[i]
        bg_spec = GC_specs_bg[i]
        galaxy = GC_cat_table.field('galaxy')[i]
        if already_fitted:
            v = GC_cat_table.field('v')[i]
            dv = GC_cat_table.field('dv')[i]
            m = GC_cat_table.field('m')[i]
            dm = GC_cat_table.field('dm')[i]
            r = GC_cat_table.field('r')[i]
            g_MUSE = GC_cat_table.field('g_MUSE')[i]
            z_MUSE = GC_cat_table.field('z_MUSE')[i]
            if miles:
                m_miles = GC_cat_table.field('m_miles')[i]
                dm_miles = GC_cat_table.field('dm_miles')[i]
            else:
                m_miles = 0
                dm_miles = 0
            if amiles:
                m_amiles = GC_cat_table.field('m_amiles')[i]
                dm_amiles = GC_cat_table.field('dm_amiles')[i]
            else:
                m_amiles = 0
                dm_amiles = 0
            if CaT:
                m_CaT = GC_cat_table.field('m_CaT')[i]
                dm_CaT = GC_cat_table.field('dm_CaT')[i]
            else:
                m_CaT = 0
                dm_CaT = 0
            if best_ssp:
                m_ssp = GC_cat_table.field('m_ssp')[i]
                dm_ssp = GC_cat_table.field('dm_ssp')[i]
            else:
                m_ssp = 0
                dm_ssp = 0
            if balmer_met:
                m_b = GC_cat_table.field('m_b')[i]
                dm_b = GC_cat_table.field('dm_b')[i]
            else:
                m_b = 0
                dm_b = 0
            if get_ages:
                age = GC_cat_table.field('age')[i]
                dage = GC_cat_table.field('dage')[i]
                m_age = GC_cat_table.field('m_age')[i]
                dm_age = GC_cat_table.field('dm_age')[i]
            else:
                age, dage, m_age, dm_age = 0, 0, 0, 0
            GCs.append(GC(ID, xpix, ypix, ra, dec, SNR=SNR, g=g, z=z, wave=wave, spec=spec, bg_spec=bg_spec, galaxy=galaxy,
                          v=v, dv=dv, m=m, dm=dm, r=r, g_MUSE=g_MUSE, z_MUSE=z_MUSE, m_miles=m_miles, dm_miles=dm_miles,
                          m_CaT=m_CaT, dm_CaT=dm_CaT, m_ssp=m_ssp, dm_ssp=dm_ssp, m_amiles=m_amiles, dm_amiles=dm_amiles, m_b=m_b, dm_b=dm_b,
                          age=age, dage=dage, dm_age=dm_age, m_age=m_age))
        else:
            GCs.append(GC(ID, xpix, ypix, ra, dec, SNR=SNR, g=g, z=z,
                          wave=wave, spec=spec, bg_spec=bg_spec, galaxy=galaxy))
    return GCs


def remove_from_cat(remove_file, GCs):
    IDs_to_remove = np.loadtxt(remove_file)

    new_GCs = []
    for i, GCi in enumerate(GCs):
        if GCi.ID in IDs_to_remove:
            print('Will remove GC {0}'.format(i))
        else:
            new_GCs.append(GCi)
    print('The new catalog contains {} GCs!'.format(len(new_GCs)))
    return new_GCs


def get_statistics_from_dist(array):
    med = np.percentile(array, 50)
    err_p = np.percentile(array, 100-16.84) - med
    err_m = med - np.percentile(array, 16.84)
    err = 0.5 * (np.round(err_p, 2) + np.round(err_m, 2))
    return np.round(med, 2), err


def der_snr(flux):
    """
    REFERENCES  * ST-ECF Newsletter, Issue #42:
                www.spacetelescope.org/about/further_information/newsletters/html/newsletter_42.html
                * Software:
                www.stecf.org/software/ASTROsoft/DER_SNR/
    AUTHOR      Felix Stoehr, ST-ECF
                24.05.2007, fst, initial import
                01.01.2007, fst, added more help text
                28.04.2010, fst, return value is a float now instead of a numpy.float64
    """
    flux = np.array(flux[np.isfinite(flux)])
    # Values that are exactly zero (padded) or NaN are skipped
    flux = flux[(flux != 0.0) & np.isfinite(flux)]
    n = len(flux)
    # For spectra shorter than this, no value can be returned
    if n > 4:
        signal = np.median(flux)
        noise = 0.6052697 * np.median(np.abs(2.0 * flux[2:n-2] - flux[0:n-4] - flux[4:n]))
        return float(signal / noise)
    else:
        return 0.0


class host_galaxy():
    def __init__(self, name, xc=0, yc=0, PA=0, eps=0, vsys=0, Reff=0, mB=0, xoffset=0, yoffset=0, fwhm=0.8, mc=0):
        self.name = name
        self.xc = xc
        self.yc = yc
        self.PA = PA
        self.eps = eps
        self.vsys = vsys
        self.Reff = Reff
        self.mB = mB
        self.xoffset = xoffset
        self.yoffset = yoffset
        self.fwhm = fwhm/0.2
        self.mc = mc


def dist_circle(xc, yc, s):
    """
    Returns an array in which the value of each element is its distance from
    a specified center. Useful for masking inside a circular aperture.

    The (xc, yc) coordinates are the ones one can read on the figure axes
    e.g. when plotting the result of my find_galaxy() procedure.

    FROM MGEFIT
    """
    x, y = np.ogrid[:s[0], :s[1]] - np.array([yc, xc])  # note yc before xc
    rad = np.sqrt(x**2 + y**2)
    return rad


def get_galaxy_info(galaxy, file='/Users/kfahrion/Documents/Data_clean/MUSE_Data/galaxies.dat'):
    f = open(file)
    for line in f.readlines()[1:]:  # skip header:
        columns = line.strip().split()
        gal_name = columns[0]
        if gal_name == galaxy:
            xc = int(columns[1])
            yc = int(columns[2])
            PA = float(columns[3])
            eps = float(columns[4])
            vsys = float(columns[5])
            Reff = float(columns[6])
            mB = float(columns[7])
            xoffset = float(columns[8])
            yoffset = float(columns[9])
            fwhm = float(columns[10])
            m_c = float(columns[11])
    try:
        return host_galaxy(galaxy, xc=xc, yc=yc, PA=PA, eps=eps, vsys=vsys, Reff=Reff, mB=mB, xoffset=xoffset, yoffset=yoffset, fwhm=fwhm, mc=m_c)
    except:
        print('{} not found in galaxies.dat!'.format(galaxy))
        return host_galaxy(galaxy)
