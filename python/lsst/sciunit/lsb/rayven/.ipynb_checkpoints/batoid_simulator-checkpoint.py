from .constants import LSSTCamConstants
from .ghost_data import Ghost, StarGhostSet, FieldGhostSet

import batoid
import numpy as np
from astropy.table import QTable

class BatoidSimulator:

    def __init__(self, obs_params, reflectance, star_table, scaling, telescope=None, verbose=False):

        self._verbose = verbose
        
        self.obs_params = obs_params

        if telescope is None:
            self.telescope = batoid.Optic.fromYaml(f"LSST_{self.obs_params.band}.yaml")
        else:
            self.telescope = telescope
        
        self.reflectance = reflectance
        
        self._validate_star_table(star_table)
        self.star_table = star_table
        self.num_stars = len(star_table)
        
        self.scaling = self.set_scaling(scaling)
        

    def _validate_star_table(self, table):
        if not isinstance(table, QTable):
            raise TypeError(f"star table must be an astropy.table.QTable object")
            
        required_cols = {'ra', 'dec', 'mag', 'flux', 'fa_x', 'fa_y'}
        if not required_cols.issubset(table.colnames):
            missing = required_cols - set(table.colnames)
            raise ValueError(f"camera geometry coordinate transform table is missing required column(s): {', '.join(missing)}")

    def set_scaling(self, scaling):

        match scaling:
            case 'constant':
                return [1]*self.num_stars
            case 'flux':
                return self.star_table['flux'].value
            case 'mag':
                return self.star_table['mag'].value
            case _:
                raise ValueError(f"Batoid scaling must be 'flux', 'mag' or 'constant', currently: {scaling}")

    def set_optic_reflectance(self):

        for surface in self.telescope.itemDict.values():

            # L1, L2, L3 reflectances
            if surface.name in self.reflectance.values.keys()and isinstance(surface, batoid.RefractiveInterface):
                surface.forwardCoating = batoid.SimpleCoating(self.reflectance.values[surface.name], 
                                                              1-self.reflectance.values[surface.name])
                
                surface.reverseCoating = batoid.SimpleCoating(self.reflectance.values[surface.name], 
                                                              1-self.reflectance.values[surface.name])
            # Filter reflectances
            elif 'Filter' in surface.name and isinstance(surface, batoid.RefractiveInterface):
                surface.forwardCoating = batoid.SimpleCoating(
                    self.reflectance.values[f'{self.obs_params.band}{surface.name[6:]}'], 
                    1-self.reflectance.values[f'{self.obs_params.band}{surface.name[6:]}'])
                
                surface.reverseCoating = batoid.SimpleCoating(
                    self.reflectance.values[f'{self.obs_params.band}{surface.name[6:]}'], 
                    1-self.reflectance.values[f'{self.obs_params.band}{surface.name[6:]}'])
                
            # other optics besides detector
            elif 'Detector' in surface.name:
                continue
            else:
                surface.forwardCoating = batoid.SimpleCoating(0.0, 1.0)
                surface.reverseCoating = batoid.SimpleCoating(0.0, 1.0)

    def set_detector_reflectance(self, dettype):

        for surface in self.telescope.itemDict.values():

            if isinstance(surface, batoid.Detector) and 'Detector' in surface.name:
                surface.forwardCoating = batoid.SimpleCoating(self.reflectance.values[dettype], 
                                                              1-self.reflectance.values[dettype])

    def simulate_single_star(self, star_index=None, scaling=None, fa_x=None, fa_y=None, dettype=None):

        # option to provide just the star index and then retrieving everything from the star_table
        if star_index is not None and None in (fa_x, fa_y, dettype):
            fa_x, fa_y, dettype = self.star_table['fa_x', 'fa_y', 'detector_type'][star_index]
            scaling = self.scaling[star_index]
            
        # option to provide all information manually
        elif None not in (fa_x, fa_y, dettype, scaling):
            star_index = None

        else:
            raise ValueError("Either star_index or all of (fa_x, fa_y, dettype, scaling) must be provided.")

        # set reflectances of the telescope optics and detector
        self.set_optic_reflectance()
        self.set_detector_reflectance(dettype)

        # set wavelength to simulate at
        wavelength = LSSTCamConstants.median_wavelengths[self.obs_params.band] * 1e-9

        rays = batoid.RayVector.asPolar(
            optic=self.telescope, wavelength=wavelength,
            theta_x=(fa_x), theta_y=(fa_y),
            nrad=300, naz=2000
        )
    
        rForward, rReverse = self.telescope.traceSplit(rays, minFlux=1e-4, _verbose=self._verbose)
        forwardFlux = np.sum([np.sum(rr.flux) for rr in rForward])
        reverseFlux = np.sum([np.sum(rr.flux) for rr in rReverse])
    
        x_arr = [rr.x for rr in rForward] 
        y_arr = [rr.y for rr in rForward]
        
        #calculate total batoid flux produced by star
        tot_flux = np.sum(np.concatenate([rr.flux for rr in rForward]))
        
        labels, flux_arr = [], []
        for ray in rForward:
            #normalize flux using scaling 
            normed_flux = self.normalize_flux(flux = ray.flux,
                                              total_flux = tot_flux,
                                              scaling = scaling)
            flux_arr.append(normed_flux)
            #label each ghost
            l = self.label_ghost(ray)
            labels.append(l)
        
        gb = self.store_simulated_data(labels, rForward, x_arr, y_arr, flux_arr, tot_flux)
        
        return gb

    def label_ghost(self, ray):
        comp_labels = LSSTCamConstants.telescope_component_labels
        comp = LSSTCamConstants.telescope_components
        
        ghost_generating_optics = []

        # find optical elements that repeat with one other element in between
        # indicating the bounce of a photon
        for i in range(len(ray.path)-2):
            if ray.path[i] == ray.path[i+2]:
                ghost_generating_optics.append(ray.path[i+1])

        # case when there are no repeated elements
        if np.size(ghost_generating_optics) == 0:
            ghost_generating_optics = ['', '']

        # labelling the ghost with shorter labels: 
        # ['Filter_entrance', 'Filter_exit'] becomes 'F2-F1' etc.
        if '' in ghost_generating_optics: # no bounce case: star flux
            label = 'Star'
        else:
            label = f'{comp_labels[comp.index(ghost_generating_optics[0])]
                }-{
                comp_labels[comp.index(ghost_generating_optics[1])]}'
        return label
        
    def normalize_flux(self, flux, total_flux, scaling):
        
        return flux/total_flux * scaling

    def store_simulated_data(self, label_arr, ray_arr, x_arr, y_arr, flux_arr, tot_flux):
        ghost_arr = []

        for l, r, x, y, f in zip(label_arr, ray_arr, x_arr, y_arr, flux_arr):
            ghost = Ghost(name = l, ray = r, x = x *1e3, y = y *1e3, flux = f) # 1e3 to convert from m to mm
            ghost_arr.append(ghost)

        star_ghost_set = StarGhostSet(ghosts = ghost_arr)
        
        return star_ghost_set
        
    def simulate_fov(self):
        star_ghost_sets = []
        for i in range(self.num_stars):
            
            ghost_set = self.simulate_single_star(star_index=i)

            star_ghost_sets.append(ghost_set)

        field_ghost_set = FieldGhostSet(star_ghost_sets = star_ghost_sets)
        return field_ghost_set
        
    

    

    
        
    
                         

            
        

                
                