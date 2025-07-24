from lsst.summit.utils import ConsDbClient
import numpy as np
import os

class ObservationParameters:

    def __init__(self, visit=None, ra=None, dec=None, band=None, zeropoint=None, exposure_catalog=None, **kwargs):
        
        if visit is not None:
            self.visit = visit
            
            if exposure_catalog is None:
                exposure_catalog = self._load_exposure_catalog()
                
            self.ra, self.dec, self.band, self.zeropoint = self._load_from_exposure_catalog(visit, exposure_catalog)
    
        elif None not in (ra, dec, band, zeropoint):
            self.visit = None
            self.ra = ra
            self.dec = dec
            self.band = band
            self.zeropoint = zeropoint
    
        else:
            raise ValueError("Either visit or all of (ra, dec, band, zeropoint) must be provided.")

    def _load_exposure_catalog(self):
        os.environ["no_proxy"] += ",.consdb"
        query = """
            SELECT v.visit_id as visit, v.band, v.exp_time, v.s_ra as ra, v.s_dec as dec, v.sky_rotation, v.img_type, v.physical_filter, q.zero_point_median
            FROM cdb_lsstcam.visit1 v, cdb_lsstcam.visit1_quicklook q
            WHERE v.visit_id = q.visit_id and q.zero_point_median IS NOT NULL"""

        client = ConsDbClient('http://consdb-pq.consdb:8080/consdb/')
        table = client.query(query)
        # outfile = f'catalogs/LSSTcam_exposure_list.csv'
        #table.write(outfile, overwrite=True)
        return table

    def _load_from_exposure_catalog(self, visit, exposure_catalog):

        row = exposure_catalog[exposure_catalog['visit']==visit]
   
        ra = row['ra'].data[0]
        dec = row['dec'].data[0]
        band = row['band'].data[0]
        zeropoint = row['zero_point_median'].data[0]

        return ra, dec, band, zeropoint

    