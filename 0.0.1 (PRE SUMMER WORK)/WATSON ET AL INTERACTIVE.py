import matplotlib.pyplot as plt
import numpy as np
import math 

Dh = 4283
tH = 4.408*10**17 #in seconds

#what is an afermative response
affirmative = ["yes", "y", "yeah", "yee", "affirmative", "si", "certainly", "quite so", "definetly", "yaeh", "yeha" "yay"]

#what is a negative response
negative = ["no", "n", "nay", "nahh", "nah", "negatory"]


#MAKE A FUNCTION CALLED 



#defining the function that produces the graph and plots spicific points
def og(Om, Ol, Ok):
    fig = plt.figure()
    DlPlot = fig.add_subplot(111)
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
  
    # MAKE A FUNCTION CALLED DL 
    D = Dl(Om, Ol, Ok, z) 
    DlPlot.plot(z , D)
    DlPlot.xlabel("Redshift (z)")
    DlPlot.ylabel("Luminosity Distance (Dl)")
    DlPlot.title("Luminosity Distance vs Redshift")
    #change the values added or substracted in xytext if the labels of points are in bad locations.
    for i in z_set: 
        DlPlot.plot(i,Dl(i), "o")
        DlPlot.annotate("(%s,%w)"%(i, Dl(i)), xy = (i, Dl(i)), xytext=(i+.1,Dl(i)-500))
    Dlplot.saveplot("LumDis, Om="+Om+", Ol="+Ol+", Ok="+Ok+", z="+z_set)
    Dlplot.show()
    
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