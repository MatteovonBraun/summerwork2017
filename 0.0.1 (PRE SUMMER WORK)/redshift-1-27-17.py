import matplotlib as plt
import math

Dh = 4283
Om = .3
Ol = .7
Ok = 0


def indefinite(z):
    return -0.315798*math.log(3.97906*z**2+2.68048*z+5.70143)+0.631596*math.log(5.27763*z + 12.2776)+1.09396*math.atan(0.870584* z + 0.293233)


def Dm(z):
    return Dh*(indefinite(z)-indefinite(0))

def Dl(z):
    return (1+z)*Dm(z)


plt.plot(Dl,z)
