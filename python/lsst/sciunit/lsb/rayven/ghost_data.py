from astropy import units as u
from dataclasses import dataclass, field
from batoid import RayVector
from typing import List

from scipy.stats import binned_statistic_2d
import numpy as np

from .constants import LSSTCamConstants

@dataclass
class Ghost:
    name: str
    ray: RayVector
    x: List[float] = field(default_factory=list)
    y: List[float] = field(default_factory=list)
    flux: List[float] = field(default_factory=list)
    #star_flux: float
    
    @property
    def x_size(self):
        return np.nanmax(self.x)-np.nanmin(self.x)

    @property
    def y_size(self):
        return np.nanmax(self.y)-np.nanmin(self.y)
        
    def bin(self, bins):
        binned_ghost, _, _, _ = binned_statistic_2d(self.x, self.y, 
                                                    range=[[LSSTCamConstants.fp_min_x.value, 
                                                            LSSTCamConstants.fp_max_x.value],
                                                           [LSSTCamConstants.fp_min_y.value, 
                                                            LSSTCamConstants.fp_max_y.value]], 
                                                    values=self.flux, statistic='sum', bins=bins)
        return np.flipud(binned_ghost.T)


    def calculate_area(self, bins=1000, units='mm'):
        
        binned_ghost = self.bin(bins = bins)
        
        x_binsize, y_binsize = LSSTCamConstants.fp_width/bins, LSSTCamConstants.fp_height/bins
        bin_area = np.count_nonzero(binned_ghost) * x_binsize * y_binsize
        
        match units:
            case 'mm':
                self.area = bin_area
            case 'pixel':
                self.area = bin_area * LSSTCamConstants.mm_to_pixel**2
            case 'arcsec':
                self.area = bin_area * (LSSTCamConstants.mm_to_pixel * LSSTCamConstants.pixel_to_arcsec)**2
            case 'deg':
                self.area = (bin_area * (LSSTCamConstants.mm_to_pixel * LSSTCamConstants.pixel_to_arcsec)**2).to(u.deg**2)
                
        
        
@dataclass
class StarGhostSet:
    ghosts: List[Ghost] = field(default_factory=list)
    
    def __getitem__(self, index):
        """
        Allow indexing by integer (position) or string (ghost name).
        """
        if isinstance(index, int):
            return self.ghosts[index]
        elif isinstance(index, str):
            for ghost in self.ghosts:
                if ghost.name == index:
                    return ghost
            raise KeyError(f"Ghost with name '{index}' not found.")
        else:
            raise TypeError("Index must be an integer or a string (ghost name).")

    def __len__(self):
        return len(self.ghosts)

    @property
    def labels(self):
        return [ghost.name for ghost in self.ghosts]
        
    @property
    def x(self):
        return np.concatenate([ghost.x for ghost in self.ghosts])

    @property
    def y(self):
        return np.concatenate([ghost.y for ghost in self.ghosts])
        
    @property
    def flux(self):
        return np.concatenate([ghost.flux for ghost in self.ghosts])

    @property
    def total_flux(self):
        return np.sum(self.flux)
    
    def index(self, name):
        return [ghost.name for ghost in self.ghosts].index(name)

    def append(self, ghost: Ghost):
        self.ghosts.append(ghost)
        
    def bin(self, bins):
        binned_ghosts, _, _, _ = binned_statistic_2d(self.x, self.y, 
                                                    range=[[LSSTCamConstants.fp_min_x.value, 
                                                            LSSTCamConstants.fp_max_x.value],
                                                           [LSSTCamConstants.fp_min_y.value, 
                                                            LSSTCamConstants.fp_max_y.value]], 
                                                    values=self.flux, statistic='sum', bins=bins)
        return np.flipud(binned_ghosts.T)
    

@dataclass
class FieldGhostSet:
    star_ghost_sets: List[StarGhostSet] = field(default_factory=list)

    def __getitem__(self, index):
        return self.star_ghost_sets[index]

    def __len__(self):
        return len(self.star_ghost_sets)

    @property
    def ghosts(self) -> List[Ghost]:
        """Flattened list of all Ghosts across all StarGhostSets."""
        return [ghost for star_set in self.star_ghost_sets for ghost in star_set.ghosts]

    @property
    def x(self):
        return np.concatenate([ghost.x for ghost in self.ghosts])

    @property
    def y(self):
        return np.concatenate([ghost.y for ghost in self.ghosts])

    @property
    def flux(self):
        return np.concatenate([ghost.flux for ghost in self.ghosts])

    @property
    def total_flux(self):
        return np.sum(self.flux)

    def append(self, star_ghost_set: StarGhostSet):
        self.star_ghost_sets.append(star_ghost_set)

    def bin(self, bins):
        binned_ghosts, _, _, _ = binned_statistic_2d(self.x, self.y, 
                                                    range=[[LSSTCamConstants.fp_min_x.value, 
                                                            LSSTCamConstants.fp_max_x.value],
                                                           [LSSTCamConstants.fp_min_y.value, 
                                                            LSSTCamConstants.fp_max_y.value]], 
                                                    values=self.flux, statistic='sum', bins=bins)
        return np.flipud(binned_ghosts.T)
    