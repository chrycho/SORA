from sora.config import input_tests
from sora.body.meta import PhysicalData
from sora.body import Body
from sora.ephem.meta import BaseEphem
from astropy.coordinates import SkyCoord
from .meta import BaseRing
import astropy.units as u
import numpy as np
import warnings

class Ring(BaseRing):    
    def __init__(self, ephem, **kwargs):
        """ 
        Parameters:
        'ring_id',
        'radius', 
        'radius_err',
        'eccentricity', 
        'eccentricity_err',
        'earth_pole_coord',
        'pole_orientation',
        'normal_opacity_err',
        'normal_optical_depth', 
        'normal_optical_depth_err',
        'radial_width', 
        'radial_width_err',
        'equivalent_depth', 
        'equivalent_depth_err',
        'equivalent_width', 
        'equivalent_width_err'
	
                              
        """
        
        allowed_kwargs = ['ring_id',
                          'radius', 'radius_err',
                          'eccentricity', 'eccentricity_err',
                          'earth_pole_coord',
                          'pole_orientation',
                          'normal_opacity', 'normal_opacity_err',
                          'normal_optical_depth', 'normal_optical_depth_err',
                          'radial_width', 'radial_width_err',
                          'equivalent_depth', 'equivalent_depth_err',
                          'equivalent_width', 'equivalent_width_err']
        
            
        input_tests.check_kwargs(kwargs, allowed_kwargs=allowed_kwargs)

        earth_pole_coord = kwargs.get('earth_pole_coord', '12 00 00.00000 +90 00 00.00000')
        earth_pole = SkyCoord(earth_pole_coord, unit=(u.hourangle, u.deg))

        self.ephem = ephem
        
        self.ring_id = kwargs.get('ring_id')      
        
        if not 'pole_orientation' in kwargs:
            raise ValueError(f'Pole orientation is not defined for ring {self.ring_id}') 
        else:
            pole = SkyCoord(kwargs.get('pole_orientation'))        
            self._pole_orientation = pole.icrs

        if not 'radius' in kwargs:
            raise ValueError(f'Radius is not defined for ring {self.ring_id}')
        else:            
            self._radius = PhysicalData('Radius', 
                                        kwargs.get('radius'), 
                                        kwargs.get('radius_err', 0.0), unit=u.km)    
        
        self._normal_opacity = PhysicalData('Normal opacity', 
                                            kwargs.get('normal_opacity', None), 
                                            kwargs.get('normal_opacity_err', 0.0))  
        
        self._normal_optical_depth = PhysicalData('Normal optical depth', 
                                                  kwargs.get('normal_optical_depth', None), 
                                                  kwargs.get('normal_optical_depth_err', 0.0), unit=u.km)
        self._radial_width = PhysicalData('Radial width', 
                                          kwargs.get('radial_width', None), 
                                          kwargs.get('radial_width_err', 0.0), unit=u.km)
        
        self._eccentricity = PhysicalData('Eccentricity', 
                                          kwargs.get('eccentricity', None), 
                                          kwargs.get('eccentricity_err', 0.0))
        
        self._equivalent_depth = PhysicalData('Equivalent depth', 
                                              kwargs.get('equivalent_depth', None), 
                                              kwargs.get('equivalent_depth_err', 0.0), unit=u.km)
        
        self._equivalent_width = PhysicalData('Equivalent width', 
                                              kwargs.get('equivalent_width', None), 
                                              kwargs.get('equivalent_width_err', 0.0), unit=u.km)
        

    def get_ring_orientation(self, time, observer='geocenter'):
        pos = self.ephem.get_position(time, observer=observer)
        ring_position_angle = pos.position_angle(self._pole_orientation).rad*u.rad
        ring_opening_angle = np.arcsin(-(np.sin(self._pole_orientation.dec)*np.sin(pos.dec) + 
                                          np.cos(self._pole_orientation.dec)*np.cos(pos.dec)*np.cos(self._pole_orientation.ra-pos.ra)
                                         )
                                       )
        return ring_position_angle.to('deg'), ring_opening_angle.to('deg')
        
    def __str__(self):
        out = []
        out.append('Pole orientation (ICRS):\n    RA: {:.2f}\n    DEC: {:.2f}\n'.format(self._pole_orientation.icrs.ra,
                                                                                        self._pole_orientation.icrs.dec)
                  )
        #out.append('Pole orientation (Barycentric Mean Ecliptic):\n    lat: {:.4f}\n    lon: {:.4f}\n'.format(self.pole_orientation.barycentricmeanecliptic.lat,
        #                                                                                  self.pole_orientation.barycentricmeanecliptic.lon)
        #          )

        out.append(self._radius.__str__())        
        out.append(self._normal_opacity.__str__())
        out.append(self._normal_optical_depth.__str__())
        out.append(self._radial_width.__str__())
        out.append(self._eccentricity.__str__())
        out.append(self._equivalent_width.__str__())  
        out.append(self._equivalent_depth.__str__())
        
        return ''.join(out)
