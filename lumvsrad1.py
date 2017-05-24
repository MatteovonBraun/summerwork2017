#Code to plot luminosity vs radius

from astropy.io import ascii
from scipy import *
import matplotlib.pyplot as plt
import numpy as np


data = ascii.read("CSVResults.csv")
dic = {
'name':data[0][:][41:113],
'logl':data[1][:][41:113],
'logler':data[2][:][41:113],
'tcent':data[3][:][41:113],
'tcper':data[4][:][41:113],
'tpeak':data[6][:][41:113],
'tpper':data[7][:][41:113],
'tpner':data[8][:][41:113],
'line':data[9][:][41:113]
}

R = []       
L = []
        
def getL():
    for i in dic['logl']:
        L.append(i)
        
print 
def getR():
    c = 29979245800
    for j in dic['tcent']:
        x = float(j)
        y = log(x * c)
        z = str(y)
        R.append(z)

     
def plotLR():
    

    plt.plot(L,R, 'g.')
    plt.axis([41,46.5,24,30.5])
    plt.xlabel('log(L)')
    plt.ylabel('log(R)')
    plt.title('Luminosity vs Radius')


getL()
getR()
plotLR()

#R=ctau
