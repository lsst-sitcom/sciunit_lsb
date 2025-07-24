import requests
import yaml
from scipy.interpolate import make_interp_spline
import os
import astropy.units as u

from .YBSC_to_DES import *


def get_colorterm_spline(colorterm_file_string, band):
    """
    Get the colorterm spline for a specific band.

    This function retrieves the colorterm spline from a specified file.

    Parameters
    ----------
    colorterm_file_string : `str`
        The base string of the colorterm file name.
    band : `str`
        The band for which to get the colorterm spline
        (e.g., 'g', 'r', 'i').

    Returns
    -------
    colorterm_spline : `dict`
        Dict containing the colorterm spline nodes and values, plus
        additional info.
    """
    colorterm_path = "https://raw.githubusercontent.com/lsst-dm/the_monster/refs/heads/main/colorterms"

    colorterm_url = os.path.join(
                colorterm_path,
                colorterm_file_string+f'_{band}.yaml',
            )

    print(f"Fetching {colorterm_url}")
    # Retrieve the file content from the URL
    response = requests.get(colorterm_url, allow_redirects=True)
    # Convert bytes to string
    content = response.content.decode("utf-8")
    assert content != '404: Not Found', f"File {colorterm_url} not found."

    # Load the yaml
    colorterms_dict = yaml.safe_load(content)

    return colorterms_dict

def apply_colorterms(source_color_flux_1, source_color_flux_2, source_flux, colorterms_dict):
    """Apply the color term spline model.
    Parameters
    ----------
    source_color_flux_1 : `np.ndarray` (N,)
        Array of source fluxes used for color (1).
    source_color_flux_2 : `np.ndarray` (N,)
        Array of source fluxes used for color (2).
    source_flux : `np.ndarray` (N,)
        Array of source fluxes to convert.
    colorterms_dict : `dict`
        Dict containing colorterm spline nodes and values
    Returns
    -------
    model_flux : `np.ndarray` (N,)
        Array of fluxes converted to target system.
    """
    mag_1 = (np.array(source_color_flux_1)*u.nJy).to_value(u.ABmag)
    mag_2 = (np.array(source_color_flux_2)*u.nJy).to_value(u.ABmag)
    mag_color = mag_1 - mag_2
    spl = make_interp_spline(colorterms_dict['nodes'], colorterms_dict['spline_values'])
    model_flux = np.array(source_flux) * np.array(spl(mag_color))
    model_flux -= colorterms_dict['flux_offset']

    # Check that things are in range: colors out of range simply should
    # not be corrected.
    bad = ((mag_color < colorterms_dict['nodes'][0]) |\
           (mag_color > colorterms_dict['nodes'][-1]) |\
           (~np.isfinite(mag_color)))
    model_flux[bad] = np.nan

    return model_flux

def ab_mag_to_njy(mags):
    """Convert AB magnitudes to flux in nanoJansky
    Parameters
    ----------
    mags : `np.ndarray`
        Array of mags to convert
    Returns
    -------
    fluxes : `np.ndarray`
        Array of fluxes in nJy
    """
    return (mags * u.ABmag).to(u.nJy).value

def njy_to_ab_mag(flux):

    #return -2.5*np.log10(flux/(3631*1e9))
    return (flux * u.nJy).to(u.ABmag).value
    
def des_to_lsst(vmag, bmv, band):

    match band:
        case 'g':
            #print('g')
            mag1 = get_des_gmag(vmag, bmv)
            mag2 = get_des_imag(vmag, bmv)
            mag_source = get_des_gmag(vmag, bmv)
        case 'r':
            #print('r')
            mag1 = get_des_gmag(vmag, bmv)
            mag2 = get_des_imag(vmag, bmv)
            mag_source = get_des_rmag(vmag, bmv)
        case 'i':
            #print('i')
            mag1 = get_des_gmag(vmag, bmv)
            mag2 = get_des_imag(vmag, bmv)
            mag_source = get_des_imag(vmag, bmv)
        case 'z':
            #print('z')
            mag1 = get_des_imag(vmag, bmv)
            mag2 = get_des_zmag(vmag, bmv)
            mag_source = get_des_zmag(vmag, bmv)
        case 'y':
            #print('y')
            mag1 = get_des_imag(vmag, bmv)
            mag2 = get_des_zmag(vmag, bmv)
            mag_source = get_des_ymag(vmag, bmv)
        case 'u':
            #print('u mags are not available yet')
            vflux = ab_mag_to_njy(vmag)
            return njy_to_ab_mag(vflux), vflux
            
    
    colorterm_file_string = "Monster_to_ComCam_band" #ComCam SynthLSST
    
    colorterms_dict = get_colorterm_spline(colorterm_file_string, band)

    flux1 = ab_mag_to_njy(mag1)
    flux2 = ab_mag_to_njy(mag2)
    flux_source = ab_mag_to_njy(mag_source)
    
    lsst_flux = apply_colorterms(flux1, flux2, flux_source, colorterms_dict)
    #print(lsst_flux)
    return njy_to_ab_mag(lsst_flux), lsst_flux
    