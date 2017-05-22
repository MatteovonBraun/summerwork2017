from astropy.cosmology import WMAP9 as cosmol
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import WMAP7
from astropy.cosmology import FlatwCDM
from scipy import * 
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import LambdaCDM




class cosmo(object):
    def __init__(self, H, Ok, Om, Ode):
        self.H = H
        self.Ok = Ok
        self.Om = Om
        self.Ode = Ode
        



cosmology = {"s":cosmo(70, 0, .3089, .6911), 
             "bs0":cosmo(70, 0,.1,.9), 
             "bs1":cosmo(70, 1, 0, 0), 
             "bs2":cosmo(70, -.5,.5 , 1), 
            } 




def lumdis(z, key):
    Ho = cosmology[key].H
    mo = cosmology[key].Om
    deo = cosmology[key].Ode
    l = cosmology[key].Ok
    cos = LambdaCDM(H0=Ho, Om0=mo, Ode0=l)    
    return cos.luminosity_distance(z)




z = np.arange(0,4,.001)
norm = lumdis(z, "s")
bs0 = lumdis(z, "bs0")
bs1 = lumdis(z, "bs1")
bs2 = lumdis(z, "bs2")
        
fig, ax = plt.subplots()

ax.loglog(z, norm, 'k--', label="Omega m = .3089 Omega lambda = .6911 Omega k = 0 (Current Best Cosmolgy)")
ax.loglog(z, bs0, 'r:', label="Omega m = .1 Omega lambda = .9 Omega k = 0")
ax.loglog(z, bs1, 'b--', label="Omega m = 0 Omega lambda = 0 Omega k = 1")
ax.loglog(z, bs2, 'g-.', label="Omega m = .5 Omega lambda = 1 Omega k = -.5")

plt.title("Luminosity Distance vs Z For Different Cosmologies", fontsize = 20)
plt.xlabel("Z (Redshift)", fontsize = 14.5)
plt.ylabel("Luminosity Distance (MParsecs)", fontsize = 14.5)


legend = ax.legend(loc="upper left")

for label in legend.get_texts():
    label.set_fontsize('small')
    
plt.show()




#Add quasar data from the email from Jon/Misty

