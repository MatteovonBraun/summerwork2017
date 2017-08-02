'''
=====================================================================================================

Lumdis-lumdis_indicator
  
Version: 1.15

Authors: McDonnell, Matteo; Moore, Marc

=====================================================================================================
'''

import sys
import matplotlib.pyplot as plt
import numpy as np
import lumdis_func as l 
import astropy.units as u
import math as m 
import pandas as pd 
import scipy.odr as o
import datetime as dt



############################################# SETUP ###################################################

"""
TO DO:

- Make sure the second dataset is correctly combining the tau and Luminosity
  (it appears different than the plot Yasaman has)

- Ask Yasaman if code is correct 
 
- Once the simulation file is totally working, import the data as a new dataframe and 
  use it as data in this file. 

- Clean up code? There is a lot of repeated stuff, we can make that cleaner. 
"""



#to stop the commandline spam
pd.options.mode.chained_assignment = None  # default='warn'

#to see how long it takes to run:
startTime = dt.datetime.now()

        
cosmology = {"s":{'H0':70, "Om0":.3089, 'Ode0':.6911, 'Ok0':0}, 
             "bs0":{'H0':70, "Om0":.1, 'Ode0':.9, 'Ok0':0}, 
             "bs1":{'H0':70, "Om0":0., 'Ode0':0., 'Ok0':1.}, 
             "bs2":{'H0':70, "Om0":.5, 'Ode0':1., 'Ok0':-.5},
             "conversion":{'H0':71, "Om0":.3, 'Ode0':.7, 'Ok0':0},
             "test":{'H0':71, "Om0":.0, 'Ode0':1., 'Ok0':0}#values gotten from http://adsabs.harvard.edu/abs/2015PASP..127...67B
            } 


#the following is the pandas setup, see commit link for details. 
#https://stackoverflow.com/questions/15705630/python-getting-the-row-which-has-the-max-value-in-groups-using-groupby
raw = pd.read_csv('CSVResults z.csv')
raw.sort_values(by='name1')


#regular data

trim = pd.DataFrame()
trim = raw[['name1','log_L_agn','log_L_agn_err','t_cent', 't_c_perr', 'z']]

def percenterror(row):
    return (row['t_c_perr']/row['t_cent'])

trim['pct_err_t'] = trim.apply(lambda row: percenterror(row), axis = 1)

idx =  trim.groupby(['name1'])['pct_err_t'].transform(min) == trim['pct_err_t']
full0 = trim[idx]



#combination of table1_ascii and mbh_measurments

raw1 = pd.read_csv('mbh_measurements.csv')
raw2 = pd.read_csv('table1_ascii.csv')

raw1.set_index('RMID')
raw2.set_index('RMID')

trim1 = raw1[['redshift','tau_cent','tau_cent_uperr','tau_cent_loerr']]
trim2 = raw2[['wLw']]

full1 = pd.DataFrame(trim1.join(trim2, how='inner'))
full1 = full1[full1.wLw != -99.0]


full2 = pd.read_csv('du16_with_z.csv')

tables = [full0, full1, full2]

################################################# CALCULATION #######################################################

#using the 'clean2+ExtCorr' values on page 22 of Bentz et al. 2013
alpha = .549
beta =  1.559   #'K' in paper

def lumdis(z, key):
    return l.lum_dis(z, **cosmology[key])
    


#source: https://stackoverflow.com/questions/26886653/pandas-create-new-column-based-on-values-from-other-columns

def log_flux(row,redshift,log_L_agn):
    z = row[redshift]
    log_dis = m.log10(float((lumdis(z, 'conversion'))/u.Mpc)*(3.08567758*10**24)) #conversion to cm
    log_F = row[log_L_agn] - 2*log_dis - m.log10(4*m.pi)
    return log_F

full0['log_F'] = full0.apply(lambda row: log_flux(row, 'z', 'log_L_agn'), axis = 1)
full1['log_F'] = full1.apply(lambda row: log_flux(row, 'redshift', 'wLw'), axis = 1)
full2['log_F'] = full2.apply(lambda row: log_flux(row, 'z', 'log_L_agn'), axis = 1)

def dis(row,t_cent):
    log_d = .5*((m.log10(row[t_cent]) - beta)/alpha + 44 - m.log10(4*m.pi) - row['log_F']) #in units of log10 cm
    if t_cent == 'tau_cent':
        log_d = .5*((m.log10(row[t_cent]*(1+row['redshift'])) - beta)/alpha + 44 - m.log10(4*m.pi) - row['log_F']) #in units of log10 cm
    d = (10**log_d)/(3.08567758*10**24) #converts to linear and cm
    return d

full0['dis_mpc'] = full0.apply(lambda row: dis(row,'t_cent'), axis = 1)
full1['dis_mpc'] = full1.apply(lambda row: dis(row,'tau_cent'), axis = 1)
full2['dis_mpc'] = full2.apply(lambda row: dis(row,'t_cent'), axis = 1)

#error in just Flux and tau
def dis_Ft_err(row, log_L_agn_err,t_c_err): 
    #This is cause there is no L error in one of the data groups
    if log_L_agn_err==0:
        err_F = 0
    else: 
        err_F = row[log_L_agn_err]
    err_r = 1/(m.log(10)*row[t_c_err])
    d = row['dis_mpc']
    
    #https://www.wolframalpha.com/input/?i=derivatibe+of+.5((+r+-b)%2Fa%2B44-log10(4pi)-F)+with+respect+to+F
    ddF = -0.5
    #https://www.wolframalpha.com/input/?i=derivatibe+of+.5((+r+-b)%2Fa%2B44-log10(4pi)-F)+with+respect+to+r
    ddr = 0.5/alpha
    
    dis_err_log = m.sqrt((ddF*err_F)**2+(ddr*err_r)**2)
    #https://www.wolframalpha.com/input/?i=derivative+10%5Ex
    dis_err = m.log(10)*d*dis_err_log #because the error in d is log but d itself is linear, 10^log10(d) just produces d 
    return dis_err
              
full0['dis_Ft_err'] = full0.apply(lambda row: dis_Ft_err(row,'log_L_agn_err','t_c_perr'), axis = 1)
full1['dis_Ft_uperr'] = full1.apply(lambda row: dis_Ft_err(row,0,'tau_cent_uperr'), axis = 1)
full1['dis_Ft_loerr'] = full1.apply(lambda row: dis_Ft_err(row,0,'tau_cent_loerr'), axis = 1)
full2['dis_Ft_uperr'] = full2.apply(lambda row: dis_Ft_err(row,'log_L_agn_err','t_c_perr'), axis = 1)
full2['dis_Ft_loerr'] = full2.apply(lambda row: dis_Ft_err(row,'log_L_agn_err','t_c_nerr'), axis = 1)

#error in flux, tau, alpha, and beta NOW USING SCATTER!
def dis_Ftab_err(row, log_L_agn_err,t_c_err):
    if log_L_agn_err==0:
        err_F = 0
    else: 
        err_F = row[log_L_agn_err]
    err_r = 1/(m.log(10)*row[t_c_err])
    #like aplha and beta, it comes from clean2 + ExtCorr
    scatter = .036
    d = row['dis_mpc']
    
    ddF = -0.5
    ddr = m.sqrt((0.5/alpha)**2 + scatter)
    
    
    dis_err_log = m.sqrt((ddF*err_F)**2+(ddr*err_r)**2)
    dis_err_full = m.log(10)*d*dis_err_log 
    return dis_err_full 


full0['dis_Ftab_err'] = full0.apply(lambda row: dis_Ftab_err(row,'log_L_agn_err','t_c_perr'), axis = 1)
full1['dis_Ftab_uperr'] = full1.apply(lambda row: dis_Ftab_err(row,0,'tau_cent_uperr'), axis = 1)
full1['dis_Ftab_loerr'] = full1.apply(lambda row: dis_Ftab_err(row,0,'tau_cent_loerr'), axis = 1)
full2['dis_Ftab_uperr'] = full2.apply(lambda row: dis_Ftab_err(row,'log_L_agn_err','t_c_perr'), axis = 1)
full2['dis_Ftab_loerr'] = full2.apply(lambda row: dis_Ftab_err(row,'log_L_agn_err','t_c_nerr'), axis = 1)

def dis_FTab_averr(row):
    err = (row['dis_Ftab_uperr']+row['dis_Ftab_loerr'])/2
    return err
    
full1['dis_FTab_averr'] = full1.apply(lambda row: dis_FTab_averr(row), axis = 1)
full2['dis_FTab_averr'] = full2.apply(lambda row: dis_FTab_averr(row), axis = 1)

#only using quasars closer than 2000 MPc

restricted = pd.DataFrame()
restricted = full0.loc[full0['dis_mpc'] <= 1000]

########################################## FIT OF PLOT ######################################################

#like below, made it a function that you can uncomment for easy changing. 

def fit():

    ####closeness of fit to model####
    print 
    print "Calculating Chi Squared Grid"
    print 
    
    #sums the individual values of Chi from the two tables and devides by the number of points
    def x_sum(tables, column):
        data = []
        n = []
        for t in tables:
            data.append(t.loc['X_sum',column])
            n.append(len(t[column]))
        
        return sum(data)/sum(n)
    
    
    #gives one value for x^2 based off of the cosmology you are compairing it with 
    def general_fit(H0, Om, Ol, Ok):
        
        def cosmol_x(row,z,err):
            key = {'H0':H0, "Om0":Om, 'Ode0':Ol, 'Ok0':Ok}
            z = row[z] 
            model = float(l.lum_dis(z, **key) /u.Mpc)
            x = ((row['dis_mpc'] - model)**2)/(row[err]**2)
            return x
        
        #makes a column with a name based off of the cosmology you are testing and produces indivitual values for the sum
        #format: 'x^2_H0_OM_OL_OK'
        current_column = 'x^2_%(H0)f_%(OM)f_%(OL)f_%(OK)f' % {'H0': H0,'OM': Om,'OL': Ol,'OK': Ok}
        
        restricted[current_column] = restricted.apply(lambda row: cosmol_x(row,'z','dis_Ftab_err'), axis = 1) 
        full0[current_column] = full0.apply(lambda row: cosmol_x(row,'z','dis_Ftab_err'), axis = 1) 
        full1[current_column] = full1.apply(lambda row: cosmol_x(row,'redshift','dis_FTab_averr'), axis = 1) 
        full2[current_column] = full2.apply(lambda row: cosmol_x(row,'z','dis_FTab_averr'), axis = 1) 
        
        restricted.loc['X_sum',current_column] = restricted[current_column].sum()
        full0.loc['X_sum',current_column] = full0[current_column].sum()
        full1.loc['X_sum',current_column] = full1[current_column].sum()
        full2.loc['X_sum',current_column] = full2[current_column].sum()
        
        #add whatever data frames to the list you want to be included in the x^2 calculation
        return x_sum([full0,full1,full2], current_column) 
        
        
        
    def chi_plot_values(frame):
        #Uses a dataframe as the list of all of the independant variables
        #and appends the chi squared value to the frame 
      
        def chisquaredfinal(row):
            x2 = general_fit(row['H0s'], row['Oms'], row['Ols'], row['Oks'])
            
            #comment out the following line to get rid of the terminal spam       
            print row.name, x2        
            return x2 
            
        frame['x^2'] = frame.apply(lambda row: chisquaredfinal(row), axis = 1)
        
        return frame
    
    #making a dataframe of a bunch of values for Om and Ol, generating the Ok so that 
    #the three of the sum to 1 (see the lambda function(woo I know how lambda functions work!))
    def IVFrame(itt):
        Om = list(np.linspace(0,1,itt))
        Ol = list(np.linspace(0,1,itt))
        
        Oms=[]
        Ols=[]
        
        for f in Om:
            for x in Ol:
                Oms.append(f)
                Ols.append(x)
            
        Oks = list(map(lambda m, l: 1.-m-l, Oms, Ols))
        
        IVs = [("H0s",70.),("Oms",Oms),("Ols", Ols),('Oks',Oks)]
        
        return pd.DataFrame.from_items(IVs)
    
    
    #Make this higher if you want a higher defintion graph
    chi_itt = 30
    
    final_x2 = chi_plot_values(IVFrame(chi_itt))
    
    print 'it took', (dt.datetime.now() - startTime), 'to generate all chi'
    print 'Drawing plot...'
    
    # look here to see what I am doing:
    # https://stackoverflow.com/questions/9152958/matplotlib-3d-plot-2d-format-for-input-data
    
    exx = np.linspace(0,1,chi_itt)
    why = np.linspace(0,1,chi_itt)
    
    
    fig, chnt = plt.subplots()
    
    xx, xy = np.meshgrid(exx,why)
    xz = np.array(final_x2['x^2'])
    xz = xz.reshape((chi_itt,chi_itt))
    
    CS = chnt.contour(xx,xy,xz)
    
    
    chnt.set_title(r'Chi Squared vs $\Omega_m$ and $\Omega_\lambda$', fontsize = 15)
    
    #subtitle code doenst work for some reason
    #fig.set_subtitle('the line from top left to bottom right is curvature 0, because the three Omegas sum to 1', fontsize = 8.)
    
    chnt.set_xlabel(r"$\Omega_m$", fontsize = 14.5)
    chnt.set_ylabel(r"$\Omega_\lambda$", fontsize = 14.5)
    chnt.clabel(CS, inline = 1, fontsize = 10, colors = 'k')
    
    #for the k = 0 guideline:
    Ohm_m = np.linspace(0,1,chi_itt)
    Ohm_l = np.linspace(1,0,chi_itt)
    chnt.plot(Ohm_m, Ohm_l,'r--')
    
    
    print 'it took', (dt.datetime.now() - startTime), 'to do everything'
    
    plt.show()
    


########################################## SIMULATION #####################################################

#TO DO

######################################### MAKING GRAPHS ######################################################


#uncomment the function at the bottom to graph everything

def biggraph():
    #Labels:    
    title = 'Luminosity  Distance  vs  Z  For  Different  Cosmologies'
    axis_x =  "$Redshift$"
    axis_y1 = 'Luminosity Distance $(MPc)$'
    axis_y2 = r"$ \frac{\tau}{\sqrt{\lambda  F}}$"
    
    z = np.arange(0,1,.001)
    norm = lumdis(z, "s")
    bs0 = lumdis(z, "bs0")
    bs1 = lumdis(z, "bs1")
    bs2 = lumdis(z, "bs2")
    bs3 = lumdis(z, "test")
    
    fig, ax1 = plt.subplots()
    
    #from data 1
    Z = full0['z']
    D = full0['dis_mpc']
    D_e= full0['dis_Ft_err']
    D_ee= full0['dis_Ftab_err'] 
    
    #from data 2
    Z1 = full1['redshift']
    D1 = full1['dis_mpc']
    D_e1= [full1['dis_Ft_loerr'],full1['dis_Ft_uperr']]
    D_ee1= [full1['dis_Ftab_loerr'],full1['dis_Ftab_uperr']] 

    #from data 3
    Z2 = full2['z']
    D2 = full2['dis_mpc']
    D_e2= [full2['dis_Ft_loerr'],full2['dis_Ft_uperr']]
    D_ee2= [full2['dis_Ftab_loerr'],full2['dis_Ftab_uperr']] 
            
    ax1.plot(z, norm, 'k', label=r"$\Omega_m = .3089$  $\Omega_\lambda = .6911$  $\Omega_k = 0$  (Current Best Cosmolgy)", linewidth = 1)
    ax1.plot(z, bs0, 'r:', label=r"$\Omega_m = .1$  $\Omega_\lambda = .9 $  $\Omega_k = 0$", linewidth = 1)
    ax1.plot(z, bs1, 'y--', label=r"$\Omega_m = 0$   $\Omega_\lambda = 0$    $\Omega_k = 1$", linewidth = 1)
    ax1.plot(z, bs2, 'g-.', label=r"$\Omega_m = .5$  $\Omega_\lambda = 1$  $\Omega_k = -.5$", linewidth = 1)
    ax1.plot(z, bs3, 'b-.', label=r"$\Omega_m = 0$  $\Omega_\lambda = 1$  $\Omega_k = 0$", linewidth = 1)
    ax1.set_xlabel(axis_x, fontsize = 14.5)
    plt.title(title, fontsize = 15)
    ax1.set_ylabel(axis_y1, fontsize = 14.5)
    
    
    "test"
    #set 

    label0 = r'Redshift Independant Distance Bentz $without$ $\alpha$ and $\beta$ error' 
    label1 = r'Redshift Independant Distance Bentz $with$ $\alpha$ and $\beta$ error'

    ax1.errorbar(Z, D,xerr=0, yerr = D_e, fmt = 'k.', ecolor='r'  , elinewidth= 1.4, label=label0,zorder=10)
    ax1.errorbar(Z, D,xerr=0, yerr = D_ee, fmt = 'k.', ecolor='#FF9900'  , elinewidth= 1.3, label=label1,zorder=9)
    
    #set 2
    
    label2 = r'Redshift Independant Distance $without$ $\alpha$ and $\beta$ error'
    label3 = r'Redshift Independant Distance $with$ $\alpha$ and $\beta$ error'

    ax1.errorbar(Z1, D1,xerr=0, yerr = D_e1, fmt = 'b.', ecolor='r'  , elinewidth= 1.4, label=label2,zorder=10)
    ax1.errorbar(Z1, D1,xerr=0, yerr = D_ee1, fmt = 'b.', ecolor='#FF9900'  , elinewidth= 1.3, label=label3,zorder=9)

    #set 3
    
    label2 = r'Redshift Independant Distance Du $without$ $\alpha$ and $\beta$ error'
    label3 = r'Redshift Independant Distance Du $with$ $\alpha$ and $\beta$ error'

    ax1.errorbar(Z2, D2,xerr=0, yerr = D_e2, fmt = 'g.', ecolor='r'  , elinewidth= 1.4, label=label2,zorder=10)
    ax1.errorbar(Z2, D2,xerr=0, yerr = D_ee2, fmt = 'g.', ecolor='#FF9900'  , elinewidth= 1.3, label=label3,zorder=9)
    
    #ax1.axis([10**-3,4,4,1.1*10**5])
    
    
    legend = ax1.legend(loc="upper left")
    
    
    for label in legend.get_texts():
        label.set_fontsize('small')
    
    plt.show()


def main(state):
    if state == 0:
        fit()
    elif state ==1:
        biggraph()

main(0)