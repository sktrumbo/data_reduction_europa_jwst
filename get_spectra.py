import pickle
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy import interpolate
from get_ephemerides import get_ephemerides
from scipy.io import readsav
import scipy.signal
from scipy.optimize import curve_fit

def gaussfit(x, a0, a1, a2, a3, a4):
    z = (x - a1)
    y = a0 * np.exp(-z ** 2 / (2 * (a2 ** 2))) + a3 + a4 * x
    return y

def spec_extract_with_headers(filename, path=''):

    specs = {}
    hdulist = fits.open(filename)
    new_cubes = readsav('europa.sav') #read in IDL save file with cubes

    if 'europa_g395h_1' in filename:
        cube = new_cubes['l1']
    elif 'europa_g395h_2' in filename:
        cube =new_cubes['l2']
    elif 'europa_g235h_1' in filename:
        cube = new_cubes['m1']
    elif 'europa_g235h_2' in filename:
        cube =new_cubes['m2']
    elif 'europa_g140h_1' in filename:
        cube =new_cubes['s1']
    elif 'europa_g140h_2' in filename:
        cube =new_cubes['s2']

    gstars = readsav('gstars.sav')
    savename = path + 'specs_' + filename[:-5] + '_new.p'
    pa_aper = hdulist[1].header['PA_APER']
    plate_scale = 0.1
    crval3 = hdulist[1].header['CRVAL3']
    cdelt3 = hdulist[1].header['CDELT3']
    grating = hdulist[0].header['GRATING']
    detector = hdulist[0].header['DETECTOR']
    target = hdulist[0].header['TARGNAME']

    if target == 'GANYMEDE':
        target = 'Ganymede'
    elif target == 'EUROPA':
        target = 'Europa'

    wvl = cdelt3 * np.arange(hdulist[1].data.shape[0]) + crval3

    if ((grating == 'G395H') & (detector == 'NRS1')):
        star = gstars['s395h1']
        wvl_star = gstars['w395h1']
        filter_left = 2.925
        filter_right = 4.04
    elif ((grating == 'G395H') & (detector == 'NRS2')):
        star = gstars['s395h2']
        wvl_star = gstars['w395h2']
        filter_left = 4.18
        filter_right = 5.18
    if ((grating == 'G235H') & (detector == 'NRS1')):
        star = gstars['s235h1']
        wvl_star = gstars['w235h1']
        filter_left = 1.7
        filter_right = 2.4
    elif ((grating == 'G235H') & (detector == 'NRS2')):
        star = gstars['s235h2']
        wvl_star = gstars['w235h2']
        filter_left = 2.46
        filter_right = 3.15
    if ((grating == 'G140H') & (detector == 'NRS1')):
        star = gstars['s140h1']
        wvl_star = gstars['w140h1']
        filter_left = 0.98
        filter_right = 1.44
    elif ((grating == 'G140H') & (detector == 'NRS2')):
        star = gstars['s140h2']
        wvl_star = gstars['w140h2']
        filter_left = 1.46
        filter_right = 1.87

    tstart = hdulist[0].header['DATE-OBS'] + ' ' + hdulist[0].header['TIME-OBS'][:-4]
    tend = tstart + '.99'
    horizons_data = get_ephemerides(target, tstart, tend)
    npa = float(horizons_data[0][6])
    ang_diam = float(horizons_data[0][3])
    central_lon = float(horizons_data[0][4])
    central_lat = float(horizons_data[0][5])
    heliocentric_v = float(horizons_data[0][11])
    geocentric_v = float(horizons_data[0][13])
    pa_target = (pa_aper - npa) - 90
    pa_target = pa_target * np.pi / 180.
    central_lon = central_lon * np.pi / 180
    central_lat = central_lat * np.pi / 180
    central_lon = -central_lon
    ang_rad = ang_diam / 2.
    pixel_rad = ang_rad / plate_scale

    #find center by maximizing flux inside circle with radius pixel_rad
    flux_square = np.nanmedian(cube, axis=0)
    flux_square[np.isnan(flux_square)] = 0
    x = np.arange(flux_square.shape[1])
    y = np.arange(flux_square.shape[0])
    scale_factor = 3
    x_new = np.linspace(0, flux_square.shape[1], scale_factor * flux_square.shape[1])
    y_new = np.linspace(0, flux_square.shape[0], scale_factor * flux_square.shape[0])
    f = interpolate.interp2d(x, y, flux_square)
    flux_square_subsampled = f(x_new, y_new)

    xx, yy = np.meshgrid(x_new, y_new)
    max_flux = -10000
    x_center = -1000
    y_center = -1000
    x_estimate = np.argmax(np.nanmedian(flux_square, axis=0))
    y_estimate = np.argmax(np.nanmedian(flux_square, axis=1))
    print('Estimate Center:', (x_estimate, y_estimate))

    for i in range(len(x_new)):
        for j in range(len(y_new)):
            dd = (xx - x_new[i])**2 + (yy - y_new[j])**2
            flux = np.nanmean(flux_square_subsampled[dd < pixel_rad**2])
            if flux > max_flux:
                max_flux = flux
                x_center = x_new[i]
                y_center = y_new[j]

    print('Center:', (x_center, y_center))
    plt.figure()
    plt.imshow(flux_square, origin='lower')
    plt.scatter(x_center, y_center, c='red')
    circle = plt.Circle((x_center, y_center), pixel_rad, fill=False, color='red', linewidth=2.0)
    ax = plt.gca()
    ax.add_patch(circle)
    plt.show()

    x_mat, y_mat = np.meshgrid(x, y)

    dist = (x_mat - x_center)**2 + (y_mat - y_center)**2
    off_target = dist > pixel_rad**2

    x_prime = ((x_mat - x_center) * plate_scale) / ang_rad
    y_prime = ((y_mat - y_center) * plate_scale) / ang_rad
    x_targ = x_prime * np.cos(pa_target) - y_prime * np.sin(pa_target)
    y_targ = x_prime * np.sin(pa_target) + y_prime * np.cos(pa_target)
    x1_prime = x_prime - (plate_scale / ang_diam)
    x2_prime = x_prime + (plate_scale / ang_diam)
    y1_prime = y_prime + (plate_scale / ang_diam)
    y2_prime = y_prime - (plate_scale / ang_diam)
    x11 = x1_prime * np.cos(pa_target) - y1_prime * np.sin(pa_target)
    x12 = x1_prime * np.cos(pa_target) - y2_prime * np.sin(pa_target)
    x21 = x2_prime * np.cos(pa_target) - y1_prime * np.sin(pa_target)
    x22 = x2_prime * np.cos(pa_target) - y2_prime * np.sin(pa_target)
    y11 = x1_prime * np.sin(pa_target) + y1_prime * np.cos(pa_target)
    y12 = x1_prime * np.sin(pa_target) + y2_prime * np.cos(pa_target)
    y21 = x2_prime * np.sin(pa_target) + y1_prime * np.cos(pa_target)
    y22 = x2_prime * np.sin(pa_target) + y2_prime * np.cos(pa_target)

    # transform into lat, lon coordinates.
    rho = np.sqrt(x_targ ** 2. + y_targ ** 2.)
    c = np.arcsin(rho)
    lats = np.arcsin((np.cos(c)) * (np.sin(central_lat)) + (y_targ * (np.cos(central_lat))))
    lons = central_lon + np.arctan2(x_targ * np.sin(c), (rho * (np.cos(c)) * (np.cos(central_lat)) - y_targ * (np.sin(c)) * (np.sin(central_lat))))
    lons = -lons

    rho11 = np.sqrt(x11 ** 2. + y11 ** 2.)
    c11 = np.arcsin(rho11)
    lat11 = np.arcsin((np.cos(c11)) * (np.sin(central_lat)) + (y11 * (np.cos(central_lat))))
    lon11 = central_lon + np.arctan2(x11 * np.sin(c11), (rho11 * (np.cos(c11)) * (np.cos(central_lat)) - y11 * (np.sin(c11)) * (np.sin(central_lat))))

    rho12 = np.sqrt(x12 ** 2. + y12 ** 2.)
    c12 = np.arcsin(rho12)
    lat12 = np.arcsin((np.cos(c12)) * (np.sin(central_lat)) + (y12 * (np.cos(central_lat))))
    lon12 = central_lon + np.arctan2(x12 * np.sin(c12), (rho12 * (np.cos(c12)) * (np.cos(central_lat)) - y12 * (np.sin(c12)) * (np.sin(central_lat))))

    rho21 = np.sqrt(x21 ** 2. + y21 ** 2.)
    c21 = np.arcsin(rho21)
    lat21 = np.arcsin((np.cos(c21)) * (np.sin(central_lat)) + (y21 * (np.cos(central_lat))))
    lon21 = central_lon + np.arctan2(x21 * np.sin(c21), (rho21 * (np.cos(c21)) * (np.cos(central_lat)) - y21 * (np.sin(c21)) * (np.sin(central_lat))))

    rho22 = np.sqrt(x22 ** 2. + y22 ** 2.)
    c22 = np.arcsin(rho22)
    lat22 = np.arcsin((np.cos(c22)) * (np.sin(central_lat)) + (y22 * (np.cos(central_lat))))
    lon22 = central_lon + np.arctan2(x22 * np.sin(c22), (rho22 * (np.cos(c22)) * (np.cos(central_lat)) - y22 * (np.sin(c22)) * (np.sin(central_lat))))

    corner_lons = (np.array((lon11, lon12, lon21, lon22))) * (-1.)  # west longitude
    corner_lats = np.array((lat11, lat12, lat21, lat22))
    corners = np.dstack((corner_lons, corner_lats))

    scale_factor_star = 5
    wvl_new = np.linspace(wvl_star[0], wvl_star[-1], scale_factor_star * len(wvl_star))

    # load in solar model
    solar = np.loadtxt('solar_kurucz.txt')
    solar_wvl = solar[:, 0] / 1000.
    solar_spec = solar[:, 1]
    solar_spec = solar_spec[((solar_wvl >= wvl[0] - .25) & (solar_wvl <= wvl[-1] + .25))]
    solar_wvl = solar_wvl[((solar_wvl > wvl[0] - .25) & (solar_wvl < wvl[-1] + .25))]

    # set up vector to hold the solar model binned to JWST exact resolution
    if (grating == 'G395H'):
        r = fits.open('jwst_nirspec_g395h_disp.fits')
    elif (grating == 'G235H'):
        r = fits.open('jwst_nirspec_g235h_disp.fits')
    elif (grating == 'G140H'):
        r = fits.open('jwst_nirspec_g140h_disp.fits')

    if ((grating == 'G395H') & (detector == 'NRS1')):
        bin_res = .000012
    elif ((grating == 'G395H') & (detector == 'NRS2')):
        bin_res = .00006
    if ((grating == 'G235H') & (detector == 'NRS1')):
        bin_res = .000012
    elif ((grating == 'G235H') & (detector == 'NRS2')):
        bin_res = .000012
    if ((grating == 'G140H') & (detector == 'NRS1')):
        bin_res = .000012
    elif ((grating == 'G140H') & (detector == 'NRS2')):
        bin_res = .000012

    c = 2.99792458e5  # speed of light in km/s

    # shift the solar spectrum to match Europa's frame
    solar_wvl_binned = np.arange(wvl[0] - .25, wvl[-1] + .25, bin_res)
    solar_binned = np.zeros(len(solar_wvl_binned))
    solar_jwst = np.zeros(len(wvl))
    solar_smoothed = np.zeros(len(solar_wvl_binned))
    r = np.vstack((r[1].data['WAVELENGTH'], r[1].data['R']))
    r_interp = np.interp(solar_wvl_binned, r[0, :], r[1, :])

    for j in range(1, len(solar_wvl_binned) - 1):
        wvl_ref = (solar_wvl_binned[j] + solar_wvl_binned[j + 1]) / 2.
        shift_solar = (wvl_ref - wvl_ref / np.sqrt((1 + heliocentric_v / c) / (1 - heliocentric_v / c)))
        shift_obs = (wvl_ref - wvl_ref / np.sqrt((1 + geocentric_v / c) / (1 - geocentric_v / c)))
        total_shift = shift_solar + shift_obs
        solar_binned[j] = np.mean(solar_spec[((solar_wvl >= solar_wvl_binned[j] - total_shift) & (
                    solar_wvl <= solar_wvl_binned[j + 1] - total_shift))])

    for i in range(0, len(solar_wvl_binned)):
        fwhm = solar_wvl_binned[i]/r_interp[i]
        sigma = fwhm/2.355
        wvl_trunc = solar_wvl_binned[((solar_wvl_binned > solar_wvl_binned[i] - 3 * sigma) & (solar_wvl_binned < solar_wvl_binned[i] + 3 * sigma))]
        solar_spec_trunc = solar_binned[((solar_wvl_binned > solar_wvl_binned[i] - 3 * sigma) & (solar_wvl_binned < solar_wvl_binned[i] + 3 * sigma))]
        gauss_filter = 1. / np.sqrt(2. * np.pi * (sigma**2.)) * np.exp(-((wvl_trunc - solar_wvl_binned[i])**2.) / (2. * (sigma**2.)))
        solar_smoothed[i] = np.trapz(gauss_filter * solar_spec_trunc, wvl_trunc)
        if np.mod(i, int(len(solar_wvl_binned)/10.)) == 0:
            print('Done ', int(i / len(solar_wvl_binned) * 100), 'percent')

    for j in range(1, len(wvl) - 1):
        solar_jwst[j] = np.mean(solar_smoothed[(solar_wvl_binned >= ((wvl[j] + wvl[j - 1]) / 2.)) & (solar_wvl_binned < ((wvl[j] + wvl[j + 1]) / 2.))])
    solar_jwst[-1] = solar_jwst[-2]
    solar_jwst[0] = solar_jwst[1]

    #correct for any small error in JWST wavelengths
    cube_norm = cube / np.nanmean(cube, axis = 0)
    spec = np.nanmean(cube_norm[:, ~off_target], axis=1)

    spec_upsampled = np.interp(wvl_new, wvl, spec)
    solar_upsampled = np.interp(wvl_new, wvl, solar_jwst)
    spec_upsampled_clipped = spec_upsampled[((~np.isnan(spec_upsampled)) & (~np.isnan(solar_upsampled)))]
    solar_upsampled_clipped = solar_upsampled[((~np.isnan(spec_upsampled)) & (~np.isnan(solar_upsampled)))]
    wvl_new_clipped = wvl_new[((~np.isnan(spec_upsampled)) & (~np.isnan(solar_upsampled)))]

    # filter out low frequencies
    b, a = scipy.signal.butter(3, 0.05, 'highpass')
    spec_filtered = scipy.signal.filtfilt(b, a, spec_upsampled_clipped)
    solar_filtered = scipy.signal.filtfilt(b, a, solar_upsampled_clipped)
    spec_filtered = spec_filtered[((wvl_new_clipped > filter_left) & (wvl_new_clipped < filter_right))]
    solar_filtered = solar_filtered[((wvl_new_clipped > filter_left) & (wvl_new_clipped < filter_right))]

    cc = np.correlate(spec_filtered / np.std(spec_filtered), solar_filtered / np.std(solar_filtered), 'full')
    shifts = np.arange(len(cc))
    shifts = shifts - np.floor((len(cc) / 2))
    cc = cc - np.nanmin(cc)
    cc = cc / np.nanmax(cc)
    cc_max = cc[np.nanargmax(cc) - 9: np.nanargmax(cc) + 10]
    shifts_max = shifts[np.nanargmax(cc) - 9: np.nanargmax(cc) + 10]

    parameters, covariance = curve_fit(gaussfit, shifts_max, cc_max)
    gaussfitdata = gaussfit(shifts_max, *parameters)
    solar_shift = shifts_max[np.nanargmax(gaussfitdata)]

    wvl_new_shifted = np.roll(wvl_new, np.int(solar_shift))
    if solar_shift >= 0:
        wvl_new_shifted[:np.int(solar_shift)] = np.nan
    else:
        wvl_new_shifted[np.int(solar_shift):] = np.nan
    wvl_shifted = np.interp(wvl, wvl_new, wvl_new_shifted)

    star_upsampled = np.interp(wvl_new, wvl_star, star)
    spec_upsampled_clipped = spec_upsampled[((np.isnan(spec_upsampled) == False) & (np.isnan(star_upsampled) == False))]
    star_upsampled_clipped = star_upsampled[((np.isnan(spec_upsampled) == False) & (np.isnan(star_upsampled) == False))]
    wvl_new_clipped = wvl_new[((np.isnan(spec_upsampled) == False) & (np.isnan(star_upsampled) == False))]

    #filter out low frequencies
    b, a = scipy.signal.butter(3, 0.05, 'highpass')
    spec_filtered = scipy.signal.filtfilt(b, a, spec_upsampled_clipped)
    star_filtered = scipy.signal.filtfilt(b, a, star_upsampled_clipped)

    spec_filtered = spec_filtered[((wvl_new_clipped > filter_left) & (wvl_new_clipped < filter_right))]
    star_filtered = star_filtered[((wvl_new_clipped > filter_left) & (wvl_new_clipped < filter_right))]

    cc = np.correlate(spec_filtered/np.std(spec_filtered), star_filtered/np.std(star_filtered), 'full')
    shifts = np.arange(len(cc))
    shifts = shifts - np.floor((len(cc)/2))
    cc = cc - np.nanmin(cc)
    cc = cc / np.nanmax(cc)
    cc_max = cc[np.nanargmax(cc) - 9: np.nanargmax(cc) + 10]
    shifts_max = shifts[np.nanargmax(cc) - 9: np.nanargmax(cc) + 10]

    parameters, covariance = curve_fit(gaussfit, shifts_max, cc_max)
    gaussfitdata = gaussfit(shifts_max, *parameters)
    star_shift = shifts_max[np.nanargmax(gaussfitdata)]
    star_upsampled_shifted = np.roll(star_upsampled, np.int(star_shift))
    if star_shift > 0:
        star_upsampled_shifted[:np.int(star_shift)] = np.nan
    else:
        star_upsampled_shifted[np.int(star_shift):] = np.nan
    star_shifted = np.interp(wvl, wvl_new, star_upsampled_shifted)
    cube_div = cube / np.reshape(star_shifted, (len(star_shifted), 1, 1))

    solar_jwst_for_div = np.zeros(len(wvl_shifted))
    for j in range(1, len(wvl_shifted) - 1):
        solar_jwst_for_div[j] = np.mean(solar_smoothed[(solar_wvl_binned >= ((wvl_shifted[j] + wvl_shifted[j - 1]) / 2.)) & (solar_wvl_binned < ((wvl_shifted[j] + wvl_shifted[j + 1]) / 2.))])
    solar_jwst_for_div[-1] = solar_jwst_for_div[-2]
    solar_jwst_for_div[0] = solar_jwst_for_div[1]
    cube_div_model = cube / np.reshape(solar_jwst, (len(solar_jwst), 1, 1))

    specs['cube'] = cube
    specs['wavelengths'] = wvl_shifted
    specs['wavelengths_original'] = wvl
    specs['cube_div'] = cube_div
    specs['cube_div_model'] = cube_div_model
    specs['star_shifted'] = star_shifted
    specs['solar_jwst'] = solar_jwst_for_div
    specs['star_shift'] = star_shift/scale_factor_star
    specs['solar_shift'] = solar_shift/scale_factor_star
    specs['lats'] = lats * 180. / np.pi
    specs['lons'] = lons * 180. / np.pi
    specs['corners'] = corners * 180. / np.pi
    specs['central_lat'] = central_lat * 180. / np.pi
    specs['central_lon'] = -central_lon * 180. / np.pi
    specs['pa_target'] = pa_target * 180. / np.pi
    specs['angular_diameters'] = ang_diam
    specs['off_target'] = off_target
    specs['x_center'] = x_center
    specs['y_center'] = y_center
    specs['pixel_radius'] = pixel_rad
    specs['npa'] = npa
    pickle.dump(specs, open(savename, "wb"))


