import astropy.units as u
import numpy as np

def __calc_coef_projecao(self):
    """ Private function that calculate the projection coeficients 
    """ 

    east    = np.array([-np.sin(self.ephem.ra), np.cos(self.ephem.ra), 0])
    north   = np.array([-np.cos(self.ephem.ra)*np.sin(self.ephem.dec), -np.sin(self.ephem.ra)*np.sin(self.ephem.dec), np.cos(self.ephem.dec)]) #vector North
    target  = np.array([np.cos(self.ephem.ra)*np.cos(self.ephem.dec), np.sin(self.ephem.ra)*np.cos(self.ephem.dec), np.sin(self.ephem.dec)])
    pt = np.array([np.cos(self.__earth_pole.ra)*np.cos(self.__earth_pole.dec), np.sin(self.__earth_pole.ra)*np.cos(self.__earth_pole.dec), np.sin(self.__earth_pole.dec)]) #vector Pole Earth

    east_pole = np.cos(self.pole_coord.dec)*np.sin(self.pole_coord.ra - self.ephem.ra)
    north_pole = np.sin(self.pole_coord.dec)*np.cos(self.ephem.dec) - np.cos(self.pole_coord.dec)*np.sin(self.ephem.dec)*np.cos(self.pole_coord.ra - self.ephem.ra)

    a = -east*np.sin(self.B)*np.sin(self.P) - north*np.sin(self.B)*np.cos(self.P) - target*np.cos(self.B)
    b = -east*np.cos(self.P) + north*np.sin(self.P)
    p = +east*np.cos(self.B)*np.sin(self.P) + north*np.cos(self.B)*np.cos(self.P) - target*np.sin(self.B)

    n1 = np.cross(pt,p)
    module_n1 = np.sqrt(np.sum(n1**2))
    n1 = n1/module_n1
    n2 = np.cross(p,n1)

    self.__coef_polo = np.array([np.dot(a,n1),np.dot(b,n1),np.dot(a,n2),np.dot(b,n2)])
    self.__coef = np.array([-east_pole/(np.sin(self.B)*np.cos(self.B)), -north_pole/(np.sin(self.B)*np.cos(self.B)), -north_pole/np.cos(self.B), +east_pole/np.cos(self.B)])

def to_ring_plane(self, ksi, eta, ksi_0=0.0, eta_0=0.0):
    """ Convert between the sky-plane and the (equatorial) ring plane 

    Parameters
    ----------
        ksi : `int, float`
            Ksi position of the ring detection in the sky plane
        eta : `int, float`
            Eta position of the ring detection in the sky plane
        ksi_0 : `int, float`
            Central ksi value, default=0.0
        eta_0  : `int, float`
            Central eta value, default=0.0

    Returns:
        xp : `float` 
            X Position in the ring plane
        yp : `float` 
            Y Position in the ring plane            
    """
    delta_ksi = ksi - ksi_0
    delta_eta = eta - eta_0

    x = delta_ksi*self.__coef[0] + delta_eta*self.__coef[1]
    y = delta_ksi*self.__coef[2] + delta_eta*self.__coef[3]
    xp = x*self.__coef_polo[0] + y*self.__coef_polo[1]
    yp = x*self.__coef_polo[2] + y*self.__coef_polo[3]
    return xp, yp

def to_sky_plane(self, xp, yp, ksi_0=0.0, eta_0=0.0):
    """ Convert between the sky-plane and the (equatorial) ring plane 

    Parameters
    ----------
        xp : `int, float`
            X Position in the ring plane
        yp : `int, float`
            Y Position in the ring plane
        ksi_0 : `int, float`
            Central ksi value, default=0.0
        eta_0 : `int, float`
            Central eta value, default=0.0

    Returns:
        ksi : `float`
            Ksi position of the ring detection in the sky plane
        eta : `float`  
            Eta position of the ring detection in the sky plane             
    """
    x = (self.__coef_polo[3]*xp - self.__coef_polo[1]*yp)/(self.__coef_polo[0]*self.__coef_polo[3] - 
                                                           self.__coef_polo[2]*self.__coef_polo[1])
    y = (yp/self.__coef_polo[3]) - (x*self.__coef_polo[2]/(self.__coef_polo[3]))

    delta_ksi = (self.__coef[3]*x - self.__coef[1]*y)/(self.__coef[0]*self.__coef[3] - 
                                                       self.__coef[2]*self.__coef[1])
    delta_eta = (y/self.__coef[3]) - (delta_ksi*self.__coef[2]/(self.__coef[3]))

    ksi = delta_ksi + ksi_0
    eta = delta_eta + eta_0
    return ksi, eta

def calc_longitude(self, ksi, eta, ksi_0=0.0, eta_0=0.0):
    """ Calculate the True longitude in the J2000 reference frame 

    Parameters
    ----------
        ksi : `int, float`
            Ksi position of the ring detection in the sky plane
        eta : `int, float`
            Eta position of the ring detection in the sky plane
        ksi_0 : `int, float`
            Central ksi value, default=0.0
        eta_0  : `int, float`
            Central eta value, default=0.0


    Returns:
        along (float): True longitude in the J2000 reference frame, in degrees
    """
    delta_ksi = ksi - ksi_0
    delta_eta = eta - eta_0

    x = delta_ksi*self.__coef[0] + delta_eta*self.__coef[1]
    y = delta_ksi*self.__coef[2] + delta_eta*self.__coef[3]
    xp = x*self.__coef_polo[0] + y*self.__coef_polo[1]
    yp = x*self.__coef_polo[2] + y*self.__coef_polo[3]
    along = np.arctan2(yp,xp)*u.rad
    return along.to('deg').value

def calc_equatorial_vel(self, ksi, eta, vel_ksi, vel_eta, ksi_0=0.0, eta_0=0.0):
    """ Calculate the radial velocity in the equatorial (ring) plane

    Parameters
    ----------
        ksi : `int, float`
            Ksi
        eta : `int, float`
            Eta
        vel_ksi : `int, float`
            Velocity Ksi
        vel_eta : `int, float`
            Velocity Eta
        ksi_0 : `int, float`
            Central ksi value, default=0.0
        eta_0 : `int, float`
            Central eta value, default=0.0

    Returns:
        v_rp (float): Radial velocity in the equatorial (ring) plane, in km/s
    """
    delta_ksi = ksi - ksi_0
    delta_eta = eta - eta_0

    x = delta_ksi*self.__coef[0] + delta_eta*self.__coef[1]
    y = delta_ksi*self.__coef[2] + delta_eta*self.__coef[3]
    xp = x*self.__coef_polo[0] + y*self.__coef_polo[1]
    yp = x*self.__coef_polo[2] + y*self.__coef_polo[3]
    rp = np.sqrt(xp**2 + yp**2)

    xf = (delta_ksi+vel_ksi)*self.__coef[0] + (delta_eta+vel_eta)*self.__coef[1]
    yf = (delta_ksi+vel_ksi)*self.__coef[2] + (delta_eta+vel_eta)*self.__coef[3]
    xpf = xf*self.__coef_polo[0] + yf*self.__coef_polo[1]
    ypf = xf*self.__coef_polo[2] + yf*self.__coef_polo[3]
    rpf = np.sqrt(xpf**2 + ypf**2)

    v_rp = rpf - rp
    return v_rp

def print_coef(self):
    print('coeficientes',self.__coef)
    print('coeficientes polo',self.__coef_polo)        