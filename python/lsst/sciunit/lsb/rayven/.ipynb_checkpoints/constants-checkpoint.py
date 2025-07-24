import astropy.units as u
from dataclasses import dataclass

from lsst.obs.lsst import LsstCam

@dataclass(frozen=True)
class LSSTCamConstants:
    telescope_components: tuple = ('L1_entrance', 'L1_exit', 'L2_entrance', 'L2_exit', 'Filter_entrance', 'Filter_exit', 'L3_entrance', 'L3_exit', 'Detector')

    telescope_component_labels: tuple = ('L11', 'L12', 'L21', 'L22', 'F1', 'F2',
                                         'L31', 'L32', 'D')
                                         
    bands: tuple = ('u', 'g', 'r', 'i', 'z', 'y')

    median_wavelengths = dict(u=372, g=481, r=622, i=756, z=868, y=975) ##nm
    
    pixel_to_arcsec: float = 0.2 * u.arcsec/u.pixel # arcsec/pixel
    mm_to_pixel: float =  1e-3/(10*1e-6) * u.pixel/u.mm # 1 px = 10 um
    
    focal_plane_bbox = LsstCam.getCamera().getFpBBox()
    fp_min_x, fp_max_x = focal_plane_bbox.getMinX() * u.mm, focal_plane_bbox.getMaxX() * u.mm
    fp_min_y, fp_max_y = focal_plane_bbox.getMinY() * u.mm, focal_plane_bbox.getMaxY() * u.mm
    fp_width, fp_height = focal_plane_bbox.getWidth() * u.mm, focal_plane_bbox.getHeight() * u.mm

    ## top 10 brightest ghosts: the only ones seen in lsstcam visits so far
    _top10: tuple = ('F2-F1', 'D-L32', 'L32-L31', 'L31-F2', 'L31-F1', 'D-L31', 'L32-F2', 'L32-F1', 'D-F2', 'D-F1')

    
