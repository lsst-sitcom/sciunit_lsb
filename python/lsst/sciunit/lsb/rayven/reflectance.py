from .constants import LSSTCamConstants

import os, sys
import numpy as np
import astropy.units as u
from astropy.table import Table

REPO_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
DATA_DIR = os.path.join(REPO_DIR, 'data', 'throughputs')

class Reflectance:
    def __init__(
        self,
        band,
        L1=None,
        L2=None,
        L3=None,
        fil=None,
        det=None, #(E2V, ITL)
        base_path=DATA_DIR
    ):
        self._base_path = base_path
        self.band = band
        self.values = {}

        self._process_detector(val=det)
        self._process_optic(val=L1, optic_name='L1')
        self._process_optic(val=L2, optic_name='L2')
        self._process_optic(val=L3, optic_name='L3')
        self._process_optic(val=fil, optic_name=f'{self.band}')


    def _process_detector(self, val, detector_name=['E2V', 'ITL']):
        
        if val is None:
            wav = LSSTCamConstants.median_wavelengths[self.band]

            for name in detector_name:
                throughput_data = Table.read(f'{os.path.join(self._base_path, name)}.dat', 
                                               format='ascii', comment='#', names=('wav', 'throughput'))
                index = np.argmin(np.abs(throughput_data['wav']-wav))
    
                self.values[name] = 1-throughput_data[index]['throughput']
            
        elif isinstance(val, (int, float)):

            self.values['E2V'] = val
            self.values['ITL'] = val
            
        elif isinstance(val, (tuple, list)) and len(val)==2:

            self.values[f'{name}_entrance'] = val[0]
            self.values[f'{name}_exit'] = val[1]

        else:
            raise ValueError("detector input must be None or int/float or tuple/list of length 2")
            
    def _process_optic(self, val, optic_name):
        
        if val is None:

            wav = LSSTCamConstants.median_wavelengths[self.band]

            throughput_data = Table.read(f'{os.path.join(self._base_path, optic_name)}.dat', 
                                           format='ascii', comment='#', names=('wav', 'throughput'))
            index = np.argmin(np.abs(throughput_data['wav']-wav))

            self.values[f'{optic_name}_entrance'] = 1-throughput_data[index]['throughput']
            self.values[f'{optic_name}_exit'] = 1-throughput_data[index]['throughput']
            
        elif isinstance(val, (int, float)):

            self.values[f'{optic_name}_entrance'] = val
            self.values[f'{optic_name}_exit'] = val
            
        elif isinstance(val, (tuple, list)) and len(val)==2:

            self.values[f'{optic_name}_entrance'] = val[0]
            self.values[f'{optic_name}_exit'] = val[1]

        else:
            raise ValueError(f"optic {optic_name} input must be None or int/float or tuple/list of length 2")
        

