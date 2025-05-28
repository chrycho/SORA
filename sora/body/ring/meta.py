import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from sora.body.meta import PhysicalData

class BaseRing():
    
    @property
    def radius(self):
        return self._radius    
    @radius.setter
    def radius(self, value):
        if isinstance(value, PhysicalData):
            radius = value
            radius.name = "Radius"
        else:
            radius = PhysicalData('Radius', value)
        if radius < 0:
            raise ValueError("Radius cannot be a negative value")
        self._radius = radius        
        
    @property
    def pole_orientation(self):
        return self._pole_orientation
    @pole_orientation.setter
    def pole_orientation(self, value):
        if value is None:
            self._pole_orientation = SkyCoord(np.nan, np.nan, unit=(u.deg, u.deg))
        else:
            self._pole_orientation = SkyCoord(value, unit=(u.hourangle, u.deg))            
            
    @property
    def normal_opacity(self):
        return self._normal_opacity    
    @normal_opacity.setter
    def normal_opacity(self, value):
        if isinstance(value, PhysicalData):
            normal_opacity = value
            normal_opacity.name = "Normal opacity"
        else:
            normal_opacity = PhysicalData('Normal opacity', value)
            
        if normal_opacity < 0 or normal_opacity > 1:
            raise ValueError("Normal opacity must be betwenn 0 and 1")
            
        self._normal_opacity = normal_opacity

    
    @property
    def normal_optical_depth(self):
        return self._normal_optical_depth    
    @normal_optical_depth.setter
    def normal_optical_depth(self, value):
        if isinstance(value, PhysicalData):
            normal_optical_depth = value
            normal_optical_depth.name = "Normal optical depth"
        else:
            normal_optical_depth = PhysicalData('Normal optical depth', value)
            
        if normal_optical_depth < 0:
            raise ValueError("Normal optical depth must be positive")
            
        self._normal_optical_depth = normal_optical_depth
        
    @property
    def radial_width(self):
        return self._radial_width    
    @radial_width.setter
    def radial_width(self, value):
        if isinstance(value, PhysicalData):
            radial_width = value
            radial_width.name = "Radial width"
        else:
            radial_width = PhysicalData('Radial width', value, unit=u.km)
            
        if radial_width < 0:
            raise ValueError("Radial width must be positive")
            
        self._radial_width = radial_width
        
    @property
    def eccentricity(self):
        return self._eccentricity    
    @eccentricity.setter
    def eccentricity(self, value):
        if isinstance(value, PhysicalData):
            eccentricity = value
            eccentricity.name = "Eccentricity"
        else:
            eccentricity = PhysicalData('Eccentricity', value)
            
        if eccentricity < 0 or eccentricity > 1:
            raise ValueError("Eccentricity must be betwenn 0 and 1")
            
        self._eccentricity = eccentricity

    @property
    def equivalent_depth(self):
        return self._equivalent_depth    
    @equivalent_depth.setter
    def equivalent_depth(self, value):
        if isinstance(value, PhysicalData):
            equivalent_depth = value
            equivalent_depth.name = "Equivalent depth"
        else:
            equivalent_depth = PhysicalData('Equivalent depth', value)
            
        if equivalent_depth < 0:
            raise ValueError("Equivalent depth must be positive")
            
        self._equivalent_depth = equivalent_depth

    @property
    def equivalent_width(self):
        return self._equivalent_width 
    @equivalent_width.setter
    def equivalent_width(self, value):
        if isinstance(value, PhysicalData):
            equivalent_width = value
            equivalent_width.name = "Equivalent width"
        else:
            equivalent_width = PhysicalData('Equivalent width', value)
            
        if equivalent_width < 0:
            raise ValueError("Equivalent width must be positive")
            
        self._equivalent_width = equivalent_width