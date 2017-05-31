from scipy import * 
import matplotlib.pyplot as plt
import numpy as np 
from astropy.io import ascii

def breakuptable(file, rows= [0,10000], columns=[0,10000]):
    raw = ascii.read(file)
    table = {}
    
    names = list(raw.colnames)
    #changes rows if specified
    rstart = rows[0]
    rend = rows[1]
    
    #changes the columns that will be added to the list if specified
    cstart = columns[0]
    cend = columns[1]
    names = names[cstart:cend] 
    
    for i in range(len(names)):
        current_column = list(raw[[names[i]][0]])
        current_name = current_column[0]
        table[current_name] = current_column[rstart:rend]
        
    return table
    

dic = breakuptable('CSVResults.csv',[41,113])
 
names = dic.keys()

#BELOW IS SPICIFIC FOR THE NEEDS OF THE PLOT


L = []
R = []


def getL():
    for i in dic['log_L_agn']:
        L.append(float(i))
        
def getR(): 
    c = 299792458
    for j in dic['t_cent']:
        x = float(j)
        r = log(x * c)
        R.append(r)


def plotLR():
    x = L
    y = R
    
    xgap = (.05*(max(x)-min(x)))
    ygap = (.1*(max(y)-min(y)))
    
    xlow = min(x) - xgap
    xhigh = max(x) + xgap
    ylow = min(y) - ygap
    yhigh = max(y) +ygap
    
    
    plt.plot(L,R, 'g.')
    plt.axis([xlow,xhigh,ylow,yhigh])
    plt.xlabel('log(L)')
    plt.ylabel('log(R)')
    plt.title('Luminosity vs Radius')


getL()
getR()
plotLR()

#R=ctau
#tau = time delay

        

