from .observation_parameters import ObservationParameters
from .constants import LSSTCamConstants
from .reflectance import Reflectance
from .bright_star_catalog import BrightStarCatalog
from .camera_geometry import CameraGeometry
from .batoid_simulator import BatoidSimulator
#from .ghost_data import Ghost

from lsst.daf.butler import Butler
from lsst.obs.lsst import LsstCam

from astropy.table import join

class GhostTool:

    def __init__(self, 
                 obs_params=None, 
                 reflectance=None,
                 bright_star_catalog=None,
                 **kwargs
                ):
        
        if obs_params is not None:
            self.obs_params = obs_params
            
        else:
            self.obs_params = ObservationParameters(
                visit=kwargs.get('visit', None),
                ra=kwargs.get('ra', None),
                dec=kwargs.get('dec', None),
                band=kwargs.get('band', None),
                zeropoint=kwargs.get('zeropoint', None),
                exposure_catalog=kwargs.get('exposure_catalog', None)
            )
        
        ######## INITIALIZING DATA PRODUCTS ########

        self.data_products = {}
        self.initialize_data_products()

        ############################################
        
        if reflectance is not None:
            self.reflectance = reflectance
            
        else:
            self.reflectance = Reflectance(
                band=self.obs_params.band,
                L1=kwargs.get('L1', None),
                L2=kwargs.get('L2', None),
                L3=kwargs.get('L3', None),
                fil=kwargs.get('fil', None),
                det=kwargs.get('det', None)
            )

        
        if bright_star_catalog is not None:
            self.bright_star_catalog = bright_star_catalog

        else:
            self.bright_star_catalog = BrightStarCatalog(
                ra=self.obs_params.ra,
                dec=self.obs_params.dec,
                band=self.obs_params.band,
                zeropoint=self.obs_params.zeropoint,
                photocalib=self.data_products['photocalib'],
                ybsc=kwargs.get('yale_catalog', None),
                table=kwargs.get('bright_star_table', None)
            )
        
        self.camera_geometry = CameraGeometry(
            camera=self.data_products['camera'],
            bright_star_catalog=self.bright_star_catalog,
            wcs=self.data_products['wcs'],
            coord_transform_table=kwargs.get('coord_transform_table', None)
        )

        ######## MAKING STAR TABLE ########
        
        self.star_table = self.make_star_table(self.bright_star_catalog, self.camera_geometry)
        
        ###################################

        self.batoid_simulator = BatoidSimulator(
            obs_params=self.obs_params,
            reflectance=self.reflectance,
            star_table=self.star_table,
            scaling=kwargs.get('batoid_scaling', 'flux'),
            verbose=False
        )

    def get_data_product(self, dataset_name, detector=94):
        visit = self.obs_params.visit
        data_id = {"visit": visit, "exposure": visit, "detector": detector}
        
        try:
            return self.butler.get(dataset_name, data_id)
        except Exception as e:
            print(f"Warning: Failed to retrieve '{dataset_name}' - {e}")
            return None
            
    def initialize_data_products(self, **kwargs):

        # camera
        self.data_products['camera'] = LsstCam.getCamera()

        if self.obs_params.visit is not None:
            # butler
            butler_dict = kwargs.get('butler_dict', {'repo': "/repo/embargo_new", 
                                                     'collections':['LSSTCam/raw/all',
                                                                    "LSSTCam/runs/nightlyValidation"],
                                                     'instrument':"LSSTCam"})
    
            self.butler = Butler(butler_dict['repo'], collections=butler_dict['collections'], 
                                            instrument=butler_dict['instrument'])
            
            # world coordinate system (WCS)
            self.data_products['wcs'] = self.get_data_product('preliminary_visit_image.wcs')
    
            # photoCalib
            self.data_products['photocalib'] = self.get_data_product('preliminary_visit_image.photoCalib')
        else:
            self.butler = None
            self.data_products['wcs'] = None
            self.data_products['photocalib'] = None
            
    def make_star_table(self, bright_star_catalog, camera_geometry):

        if not len(bright_star_catalog.table) == len(camera_geometry.coord_transform_table):
            raise ValueError(f"Table size mismatch: bright_star_catalog.table has {len(bright_star_catalog.table)} rows, camera_geometry.coord_transform_table has {len(camera_geometry.coord_transform_table)} rows.")

        table = join(bright_star_catalog.table, camera_geometry.coord_transform_table, join_type='inner')
        
        return table

        


        
        