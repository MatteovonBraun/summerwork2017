from astropy.cosmology import WMAP9 as cosmol
from astropy.cosmology import FlatLambdaCDM

#Change this so that you can have different k (curved or flat)


from astropy.cosmology import WMAP7
from astropy.cosmology import FlatwCDM
from scipy import * 
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np


class cosmo(object):
    def __init__(self, H, Om, Ode):
        self.H = H
        self.Om = Om
        self.Ode = Ode


cosmology = {"s":cosmo(70, .3089, .6911), 
             "up":cosmo(70,.4,.6), 
             "lo":cosmo(69, .7, .2), 
             "bs":cosmo(71, .3, 0), 
            } 


def lumdis(z, key):
    Ho = cosmology[key].H
    mo = cosmology[key].Om
    deo = cosmology[key].Ode
    cos = FlatLambdaCDM(H0=Ho, Om0=mo)    
    return cos.luminosity_distance(z)

z = np.arange(0,4,.001)
normie = lumdis(z, "s")
up = lumdis(z, "up")
lo = lumdis(z, "lo")
bs = lumdis(z, "bs")
    
    
#ax = plt.subplots()
fig, ax = plt.subplots(figsize=(7,7))

ax.loglog(z, normie, 'k-.', label="Omega m = .3089 Omega lambda = .6911")
ax.loglog(z, up, 'r:', label="Omega m = .4 Omega lambda = .6")
ax.loglog(z, lo, 'b--', label="Omega m = .7 Omega lambda = .2")
ax.loglog(z, bs, 'g-.', label="Omega m = .3 Omega lambda = 0")

plt.title("Luminosity Distance vs Z For Different Cosmologies")
plt.xlabel("Z (Redshift)")
plt.ylabel("Luminosity Distance (MParsecs)")

#Make labels larger 


plt.figure

legend = ax.legend(loc="upper left")


for label in legend.get_texts():
    label.set_fontsize('small')
    
plt.show()






