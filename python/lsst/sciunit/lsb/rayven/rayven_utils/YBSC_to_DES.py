import numpy as np

def get_ybsc_rc(vmag, bmv):
    rmag = vmag - 0.59*(bmv) + 0.01
    return rmag

def get_ybsc_ic(vmag, bmv):
    imag = vmag - 1.244*(bmv) + 0.382
    return imag
    
def get_des_gmag(vmag, bmv):
    """
    Vectorized DES g-band magnitude calculation from Vmag and B-V color.
    
    Parameters
    ----------
    vmag : array_like
        V-band magnitudes
    bmv : array_like
        B-V colors
        
    Returns
    -------
    array_like
        DES g-band magnitudes
    """
    # Convert inputs to numpy arrays
    vmag = np.asarray(vmag)
    bmv = np.asarray(bmv)
    
    # Initialize output array with NaNs
    gmag = np.full_like(vmag, np.nan)
    
    # Create masks for each condition
    mask1 = (bmv > -0.2) & (bmv <= 0.4)
    mask2 = (bmv > 0.4) & (bmv <= 2.2)
    
    # Apply transformations
    gmag[mask1] = vmag[mask1] + 0.552 * bmv[mask1] - 0.099
    gmag[mask2] = vmag[mask2] + 0.493 * bmv[mask2] - 0.067
    
    return gmag

def get_des_rmag(vmag, bmv):
    """
    Vectorized DES R-band magnitude calculation from Rmag and R-I color.
    
    Parameters
    ----------
    vmag : array_like
        V-band magnitudes
    bmv : array_like
        B-V colors
        
    Returns
    -------
    array_like
        DES r-band magnitudes
    """
    # Convert inputs to numpy arrays

    #rmag = vmag - 0.59*(bmv) + 0.01  ##from Caldwell 1993
    
    # rmag = np.asarray(rmag)
    vmag = np.asarray(vmag)
    bmv = np.asarray(bmv)
    
    # Initialize output array with NaNs
    Rmag = np.full_like(vmag, np.nan)
    
    # Create masks for each condition
    mask1 = (bmv > -0.2) & (bmv <= 2.2)
    #mask2 = (bmv > 0.7) & (bmv <= 2.)
    
    # Apply transformations
    Rmag[mask1] = vmag[mask1] -0.543 * bmv[mask1] + 0.128
    #Rmag[mask2] = rmag[mask2] + 0.127 * bmv[mask2] + 0.113
    
    return Rmag

def get_des_imag(vmag, bmv):
    """
    Vectorized DES I-band magnitude calculation from Imag and R-I color.
    
    Parameters
    ----------
    imag : array_like
        I-band magnitudes
    rmv : array_like
        R-I colors
        
    Returns
    -------
    array_like
        DES g-band magnitudes
    """
    #imag = vmag - 1.244*(bmv) + 0.382 ##Caldwell 1993
    
    # Convert inputs to numpy arrays
    vmag = np.asarray(vmag)
    bmv = np.asarray(bmv)
    
    # Initialize output array with NaNs
    Imag = np.full_like(vmag, np.nan)
    
    # Create masks for each condition
    mask1 = (bmv > -0.2) & (bmv <= 2.2)
    #mask2 = (bmv > 0.7) & (bmv <= 2.)
    
    # Apply transformations
    Imag[mask1] = vmag[mask1] -1.04 * bmv[mask1] + 0.312
    #Imag[mask2] = imag[mask2] + 0.049 * bmv[mask2] + 0.416
    
    return Imag

def get_des_zmag(vmag, bmv):
    """
    Vectorized DES Z-band magnitude calculation from Imag and R-I color.
    
    Parameters
    ----------
    imag : array_like
        I-band magnitudes
    bmv : array_like
        R-I colors
        
    Returns
    -------
    array_like
        DES g-band magnitudes
    """
    #imag = vmag - 1.244*(bmv) + 0.382 ##Caldwell 1993
    
    # Convert inputs to numpy arrays
    vmag = np.asarray(vmag)
    bmv = np.asarray(bmv)
    
    # Initialize output array with NaNs
    zmag = np.full_like(vmag, np.nan)
    
    # Create masks for each condition
    mask1 = (bmv > -0.2) & (bmv <= 2.2)
    # mask2 = (bmv > 0.2) & (bmv <= 0.7)
    # mask3 = (bmv > 0.7) & (bmv <= 2.)
    
    # Apply transformations
    zmag[mask1] = vmag[mask1] - 1.302 * bmv[mask1] + 0.417
    
    return zmag

def get_des_ymag(vmag, bmv):
    """
    Vectorized DES Y-band magnitude calculation from Imag and R-I color.
    
    Parameters
    ----------
    imag : array_like
        I-band magnitudes
    rmv : array_like
        R-I colors
        
    Returns
    -------
    array_like
        DES g-band magnitudes
    """
    #imag = vmag - 1.244*(bmv) + 0.382 ##Caldwell 1993
    
    # Convert inputs to numpy arrays
    vmag = np.asarray(vmag)
    bmv = np.asarray(bmv)
    
    # Initialize output array with NaNs
    ymag = np.full_like(vmag, np.nan)
    
    # Create masks for each condition
    mask1 = (bmv > -0.2) & (bmv <= 2.2)
    # mask2 = (bmv > 0.2) & (bmv <= 0.7)
    # mask3 = (bmv > 0.7) & (bmv <= 2.)
    
    # Apply transformations
    ymag[mask1] = vmag[mask1] - 1.416 * bmv[mask1] + 0.504
    # ymag[mask2] = imag[mask2] - 0.831 * bmv[mask2] + 0.698
    # ymag[mask3] = imag[mask3] - 0.437 * bmv[mask2] + 0.451
    return ymag