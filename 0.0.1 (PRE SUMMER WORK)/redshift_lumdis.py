import matplotlib.pyplot as plt
import numpy as np

Dh = 4283
Om = .3
Ol = .7
Ok = 0

def indefinite(z):
    return -0.315798*np.log(3.97906*z**2+2.68048*z+5.70143)+0.631596*np.log(5.27763*z + 12.2776)+1.09396*np.arctan(0.870584* z + 0.293233)

def Dl(z):
    return (1+z)*Dh*(indefinite(z)-indefinite(0))

z = np.arange(0.0, 9, .01)
y = Dl(z)

z_set=[.5,1,6]

plt.plot(z , y)
plt.xlabel("Redshift (z)")
plt.ylabel("Luminosity Distance(D_l)")
plt.title("Luminosity Distance vs Redshift")

for i in z_set:
    plt.plot(i,Dl(i),'o')
    plt.annotate("(%s,%e)"%(i,Dl(i)) , xy = (i,Dl(i)), xytext=(i+.1,Dl(i)-500))
plt.savefig("LumDis v z")    
plt.show()

