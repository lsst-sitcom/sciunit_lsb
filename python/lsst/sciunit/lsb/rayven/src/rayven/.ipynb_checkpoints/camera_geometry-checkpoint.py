import numpy as np
import astropy.units as u
from astropy.table import Table, QTable
import lsst.afw.cameraGeom as afwCameraGeom
from lsst.geom import Point2D

class CameraGeometry:

    def __init__(self, camera, bright_star_catalog, wcs=None, coord_transform_table=None):

        self.camera = camera
        self.bright_star_catalog = bright_star_catalog
        self.wcs = wcs
        
        self.fp_geom = {} #focal plane geometry, units: mm

        self.get_focal_plane_extent()
        
        det_geom_dict = self.get_detector_extent()
        self.det_geometry_table = self.make_det_geometry_table(det_geom_dict=det_geom_dict)

        if coord_transform_table is not None:
            self._validate_coord_transform_table(coord_transform_table)
            self.coord_transform_table = coord_transform_table
        elif coord_transform_table is None and wcs is not None:
            coord_transform_dict = self.get_star_coordinate_transforms()
            self.coord_transform_table = self.make_coord_transform_table(coord_transform_dict=coord_transform_dict)
        else:
            raise ValueError("Either 'wcs' or 'coord_transform_table' must be provided to CameraGeometry.")

        self.match_star_to_detector()
            
    def get_focal_plane_extent(self):
        focal_plane_bbox = self.camera.getFpBBox()

        self.fp_geom['min_x'], self.fp_geom['max_x'] = focal_plane_bbox.getMinX(), focal_plane_bbox.getMaxX()
        self.fp_geom['min_y'], self.fp_geom['max_y'] = focal_plane_bbox.getMinY(), focal_plane_bbox.getMaxY()
        self.fp_geom['width'], self.fp_geom['height'] = focal_plane_bbox.getWidth(), focal_plane_bbox.getHeight()

    def get_detector_extent(self):
        detector_ids = []
        detector_types = []
        center_fp_x, center_fp_y = [], []
        center_fa_x, center_fa_y = [], []
        minimum_x, minimum_y = [], []
        maximum_x, maximum_y = [], []
        for detector in self.camera:
            min_x, min_y = float('inf'), float('inf')
            max_x, max_y = float('-inf'), float('-inf')
            
            corners = detector.getCorners(afwCameraGeom.FOCAL_PLANE)
            for corner in corners:
                min_x = min(min_x, corner.getX())
                min_y = min(min_y, corner.getY())
                max_x = max(max_x, corner.getX())
                max_y = max(max_y, corner.getY())

            detector_ids.append(detector.getId())
            detector_types.append(detector.getPhysicalType())
            center_fa_x.append(detector.getCenter(afwCameraGeom.FIELD_ANGLE).getX())
            center_fa_y.append(detector.getCenter(afwCameraGeom.FIELD_ANGLE).getY())
            center_fp_x.append(detector.getCenter(afwCameraGeom.FOCAL_PLANE).getX())
            center_fp_y.append(detector.getCenter(afwCameraGeom.FOCAL_PLANE).getY())
            minimum_x.append(min_x)
            minimum_y.append(min_y)
            maximum_x.append(max_x)
            maximum_y.append(max_y)

        return {'detector': detector_ids,
                'detector_type': detector_types,
                'center_fa_x': center_fa_x,
                'center_fa_y': center_fa_y,
                'center_fp_x': center_fp_x,
                'center_fp_y': center_fp_y,
                'min_x': minimum_x,
                'min_y': minimum_y,
                'max_x': maximum_x,
                'max_y': maximum_y
               }
        
    def make_det_geometry_table(self, det_geom_dict):
        
        columns = [det_geom_dict[key] for key in det_geom_dict.keys()]
        colnames = list(det_geom_dict.keys())
        units = [None, None, u.deg, u.deg, u.mm, u.mm, u.mm, u.mm, u.mm, u.mm]

        table = QTable(data=columns, names=colnames, units=units)

        return table

    def get_star_coordinate_transforms(self):
        ## transform from ra, dec to pixel coordinates
        ra, dec = self.bright_star_catalog.table['ra'].value, self.bright_star_catalog.table['dec'].value

        px_x, px_y = self.wcs.skyToPixelArray(ra, 
                                              dec, 
                                              degrees=True)
        
        central_detector = self.camera[94] 

        ## Transform from pixel coordinates to focal plane coords in mm
        tx_fp = central_detector.getTransform(afwCameraGeom.PIXELS, afwCameraGeom.FOCAL_PLANE)
        fpx, fpy = tx_fp.getMapping().applyForward(np.vstack((px_x, px_y)))
        fp_x, fp_y = fpx.ravel(), fpy.ravel()
        
        ## Transform from focal plane coords to RA, DEC offset from boresight
        tx_fa = central_detector.getTransform(afwCameraGeom.PIXELS, afwCameraGeom.FIELD_ANGLE)
        fax, fay = tx_fa.getMapping().applyForward(np.vstack((px_x, px_y)))
        fa_x,fa_y = fax.ravel(), fay.ravel()

        return {'ra': ra, 
                'dec': dec, 
                'px_x': px_x, 
                'px_y': px_y, 
                'fp_x': fp_x, 
                'fp_y': fp_y, 
                'fa_x': fa_x, 
                'fa_y': fa_y
               }
    
    def make_coord_transform_table(self, coord_transform_dict):

        columns = [coord_transform_dict[key] for key in coord_transform_dict.keys()]
        colnames = list(coord_transform_dict.keys())
        units = [u.deg, u.deg, u.pixel, u.pixel, u.mm, u.mm, u.rad, u.rad]
        
        table = QTable(data=columns, names=colnames, units=units)

        return table

    def _validate_coord_transform_table(self, table):
        if not isinstance(table, QTable):
            raise TypeError(f"camera geometry coordinate transform table must be an astropy.table.QTable object")
            
        required_cols = {'ra':u.deg, 'dec':u.deg, 'fa_x':u.rad, 'fa_y':u.rad}
        if not set(required_cols.keys()).issubset(table.colnames):
            missing = required_cols - set(table.colnames)
            raise ValueError(f"camera geometry coordinate transform table is missing required column(s): {', '.join(missing)}")

        for column, unit in required_cols.items():
            if table[column].unit is None:
                table[column].unit = unit

    def match_star_to_detector(self):
        detector, detector_type = [], []
        for fa_x, fa_y in self.coord_transform_table['fa_x', 'fa_y']:
            
            star = Point2D(fa_x.value, fa_y.value)
            
            min_dist = np.inf
            min_detid = np.nan
            min_dettype = np.nan
            
            for det in self.camera:
                if det.getId() > 188: ##do not consider AOS detectors
                    continue
                    
                det_center = det.getCenter(afwCameraGeom.FIELD_ANGLE)
                dist = star.distanceSquared(det_center)
    
                if dist < min_dist:
                    min_dist = dist
                    min_detid = det.getId()
                    min_dettype = det.getPhysicalType()
                    
            detector.append(min_detid)
            detector_type.append(min_dettype)
        
        self.coord_transform_table['detector_id'] = detector
        self.coord_transform_table['detector_type'] = detector_type
        