import os
import numpy as np
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, QTable

from rayven_utils.DES_to_LSST import des_to_lsst

REPO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
CAT_PATH = os.path.join(REPO_DIR, 'data', 'ybsc_v5.fits')

class BrightStarCatalog:

    def __init__(self, 
                 ra, 
                 dec, 
                 band, 
                 zeropoint, 
                 photocalib=None, 
                 ybsc=None,
                 table=None, 
                 base_path=CAT_PATH):
        
        self._path = base_path
        self.ra = ra
        self.dec = dec
        self.band = band
        self.zeropoint = zeropoint
        self.photocalib = photocalib

        if ybsc is None:
            self.ybsc = self.get_ybsc()
            self.filter_ybsc()
        else:
            self.ybsc = ybsc
            self.filter_ybsc()
            
        if table is None:
            ra, dec, mag, instflux = self.get_bright_stars()
            self.table = self.make_table(ra, dec, mag, instflux)
        else:
            self._validate_bright_star_table(table)
            self.table = table

    
    def _validate_bright_star_table(self, table):
        if not isinstance(table, QTable):
            raise TypeError(f"bright star catalog must be an astropy.table.QTable object")

        required_cols = {'ra':u.deg, 'dec':u.deg, 'mag':u.mag, 'flux':u.ct}
        
        if not set(required_cols.keys()).issubset(table.colnames):
            missing = required_cols - set(table.colnames)
            raise ValueError(f"bright star catalog table is missing required column(s): {', '.join(missing)}")

        for column, unit in required_cols.items():
            if table[column].unit is None:
                table[column].unit = unit

    def get_ybsc(self):
        ybsc = fits.open(self._path)
        ybsc = ybsc[1].data
        return ybsc

    def filter_ybsc(self):
        mask = np.isfinite(self.ybsc['coord_ra']) & np.isfinite(self.ybsc['coord_dec'])
        self.ybsc = self.ybsc[mask]

    def get_bright_stars(self, flux_threshold=1e7):
        boresight = SkyCoord(ra=[self.ra]*u.degree, dec=[self.dec]*u.degree)
        yale_catalog = SkyCoord(ra=self.ybsc['coord_ra']*u.degree, dec=self.ybsc['coord_dec']*u.degree)
                
        ybsc_idx, _, d2d, d3d = boresight.search_around_sky(yale_catalog, 1.9*u.degree)
        
        vmag, bmv = self.ybsc[ybsc_idx]['Vmag'], self.ybsc[ybsc_idx]['B-V']
        lsst_mag, lsst_flux = des_to_lsst(vmag, bmv, self.band)

        nan_mask = np.isnan(lsst_flux)
        lsst_mag[nan_mask] = vmag[nan_mask]
        lsst_flux[nan_mask] = 10**(-(vmag[nan_mask]-self.zeropoint)/2.5)

        mask = (lsst_flux>flux_threshold) & (np.isfinite(lsst_flux))
        ra = self.ybsc[ybsc_idx]['coord_ra'][mask]
        dec = self.ybsc[ybsc_idx]['coord_dec'][mask]
        flux = lsst_flux[mask]
        mag = lsst_mag[mask]

        if self.photocalib is not None:
            instflux = [self.photocalib.magnitudeToInstFlux(m) for m in mag]
        else:
            instflux = [10**(-(m-self.zeropoint)/2.5) for m in mag]

        return ra, dec, mag, instflux
            
        

    def make_table(self, ra, dec, mag, instflux):
        
        columns = [ra, dec, mag, instflux]
        colnames = ['ra', 'dec', 'mag', 'flux']
        units = [u.deg, u.deg, u.mag, u.ct]

        table = QTable(data=columns, names=colnames, units=units)

        return table

    