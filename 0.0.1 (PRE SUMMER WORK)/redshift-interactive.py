import matplotlib.pyplot as plt
import numpy as np
import math 

Dh = 4283
tH = 4.408*10**17 #in seconds

#what is an afermative response
affirmative = ["yes", "y", "yeah", "yee", "affirmative", "si", "certainly", "quite so", "definetly", "yaeh", "yeha" "yay"]

#what is a negative response
negative = ["no", "n", "nay", "nahh", "nah", "negatory"]

#neumericly integrating with simpson's rule
def simpson(f, a, b, n):
    if n%2==1:
        raise ValueError("n must be even (received n=%w)" % n)
    h = (b-a)/n
    s = f(a) + f(b)
    for i in range(1, n, 2):
        s += 4*f(a+i*h)
    for i in range(2, n-1, 2):
        s += 2*f(a+i*h)
    return s*h/3

def E(Om, Ol, Ok, z):
    return (Om(1+z)**3+Ok(1+z)**2+Ol)**.5

def Dl(Om, Ol, Ok, z):
    if Ok > 0:
        return Dh*Ok**(-.5)*np.sinh(Ok^(.5)*simpson(1/E(Om, Ol, Ok, z), 0, z, 50, z)/Dh)
    elif Ok == 0:
        return (1+z)*Dh*simpson(1/E(Om, Ol, Ok, z), 0, z, 50*z) 
    elif Ok < 0:
        return Dh*(np.absolute(Ok))**(-.5)*np.sin((np.absolute(Ok))^(.5)*simpson(1/E(Om, Ol, Ok, z), 0, z, 50, z)/Dh)
    else:
        print("error: Omega K is non-number")

def tL(Om, Ol, Ok, z):
    E = (Om(1+z)**3+Ok(1+z)**2+Ol)**.5
    return tH*simpson(1/(E(z)*(1+z)), 0, z, 50*z)

#defining the function that produces the graph and plots spicific points
def og(Om, Ol, Ok):
    fig = plt.figure()
    DlPlot = fig.add_subplot(111)
    tLPlot = fig.add_subplot(111)
    length = int(input("How many z values would you like to input? "))
    z_set = [] 
    for num in range(length):
        zed = "z value %a: "%(num+1)
        zquest = input(zed)
        z_set.append(float(zquest))
#This next thing is to have enough space on the graph for the label of a point
    if max(z_set) <= 3:
        z_max = 6
    else:
        z_max = math.ceil(max(z_set)+3)
    z = np.arange(0.0, z_max, .01)
    D = Dl(Om, Ol, Ok, z) 
    t = tL(Om, Ol, Ok, z)
    DlPlot.plot(z , D)
    tLPlot.plot(z , t)
    DlPlot.xlabel("Redshift (z)")
    tLPlot.xlabel("Redshift (z)")
    DlPlot.ylabel("Luminosity Distance (Dl)")
    tLPlot.ylabel("Lookback Time (tL)")
    DlPlot.title("Luminosity Distance vs Redshift")
    tLPlot.title("Lookback Time vs Redshift")
    #change the values added or substracted in xytext if the labels of points are in bad locations.
    for i in z_set: 
        DlPlot.plot(i,Dl(i), "o")
        DlPlot.annotate("(%s,%w)"%(i, Dl(i)), xy = (i, Dl(i)), xytext=(i+.1,Dl(i)-500))
        tLPlot.plot(i,Dl(i), "o")
        tLPlot.annotate("(%s,%w)"%(i, tL(i)), xy = (i, tL(i)), xytext=(i+.1,tL(i)-500))
    Dlplot.saveplot("LumDis, Om="+Om+", Ol="+Ol+", Ok="+Ok+", z="+z_set)
    tLplot.saveplot("LookbackTime, Om="+Om+", Ol="+Ol+", Ok="+Ok+", z="+z_set)
    Dlplot.show()
    tlPlot.show()
    
#asking what the omega value is 
def omegaprompt(a):
    return float(input("What value would you like for " + a +"? "))
    
def main():
    done = 0
    while done < 5:
        response=input("Do you want to specify Omega Values? ")
        if response.lower() in affirmative:
            done = 6
            Om = omegaprompt("Omega M")
            Ol = omegaprompt("Omega L")
            Ok = omegaprompt("Omega K")
            og (Om, Ol, Ok)
        elif response.lower() in negative:
            done = 6
            Om = .3
            Ol = .7
            Ok = 0 
            og (Om, Ol, Ok)
        else:
            print("Could you repeat that?")
# This is in case the while loop freaks out and loops really fast
        done =  done + 1 

main()