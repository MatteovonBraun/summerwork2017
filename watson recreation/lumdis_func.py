import math as m 
import numpy as np
import astropy.units as u

def lum_dis(z, H0 = 70. , Om0 = .3089, Ode0 = .6911, Ok0 = 0):
 
        
    _inv_efunc_scalar = lambda z, M, L, K : 1/m.sqrt(M*(1.+z)**3+K*(1.+z)**2+L)
    _inv_efunc_scalar_args = (Om0, Ode0, Ok0)
            
    
    def hubble_distance():
            
            _hubble_distance =  299792/ H0      #speed of light in km/s / H0 in km/s/MPc
            
            return _hubble_distance

    def _comoving_distance_z1z2(z1, z2):
    
    
            from scipy.integrate import quad
            f = lambda z1, z2: quad(_inv_efunc_scalar, z1, z2,
                                 args=_inv_efunc_scalar_args)[0]
                
            err = lambda z1, z2: quad(_inv_efunc_scalar, z1, z2,
                                 args=_inv_efunc_scalar_args)[1]
                                 
                             
            return hubble_distance() * vectorize_if_needed(f, z1, z2)
    

    def _comoving_transverse_distance_z1z2(z1, z2):
    
           
            dc = _comoving_distance_z1z2(z1, z2)
            if Ok0 == 0:
                return dc
            sqrtOk0 = m.sqrt(abs(Ok0))
            dh = hubble_distance()
            if Ok0 > 0:
                return dh / sqrtOk0 * np.sinh(sqrtOk0 * dc/ dh)
            else:
                return dh / sqrtOk0 * np.sin(sqrtOk0 * dc/ dh)
    
    
    def comoving_transverse_distance(z):
           
            return _comoving_transverse_distance_z1z2(0, z)
    
    
    def luminosity_distance(z):
           
    
            return (1. + z) * comoving_transverse_distance(z)
            
    return luminosity_distance(z)  *   u.Mpc
    
def vectorize_if_needed(func, *x):
    """ Helper function to vectorize functions on array inputs"""
    if any(map(isiterable, x)):
        return np.vectorize(func)(*x)
    else:
        return func(*x)
        
def isiterable(obj):
    """Returns `True` if the given object is iterable."""

    try:
        iter(obj)
        return True
    except TypeError:
        return False        
        
cosmology = {"s":{'H0':70, "Om0":.3089, 'Ode0':.6911, 'Ok0':0}, 
             "bs0":{'H0':70, "Om0":.1, 'Ode0':.9, 'Ok0':0}, 
             "bs1":{'H0':70, "Om0":0., 'Ode0':0., 'Ok0':1.}, 
             "bs2":{'H0':70, "Om0":.5, 'Ode0':1., 'Ok0':-.5},
             "conversion":{'H0':71, "Om0":.3, 'Ode0':.7, 'Ok0':0}#values gotten from http://adsabs.harvard.edu/abs/2015PASP..127...67B
            } 
            
"""

use the syntax "lum_dis(z, **cosmology[key])"

ex:

lum_dis(1, **cosmology['s'])
            
"""