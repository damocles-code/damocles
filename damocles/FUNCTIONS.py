#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 16:45:22 2018

@author: maria
"""

###SMALL MISCELLANEOUS FUNCTIONS USED IN ALL MY CODES###

import numpy as np

import sys
#sys.path.append('/home/maria/PYTHONCODES/mazzyfuncts/modules')


import re
import os
#import matplotlib.pyplot as plt
import sys
import fileinput
import tkinter as tk
from scipy.ndimage import gaussian_filter1d
####read in data file and split each column into arrays!! 



def appendline(read,start,stop):    
    x = []    
    for line in read[start:stop]:    
        line = line.split()
        x.append(line)
    return x
 
def checkspacing(array):
    #returns True if spacing is equal between every point, used to test meshgrids have been made properly
    y = []
    for first, second in zip(array, array[1:]):
        x = abs(second - first)
        print(x)
        if x != 0:
            y.append(second - first)
            
    if len(set(y)) <= 1:
        print("GRID SPACED EVENLY")   
    else:
        print("ITS NOT EQUAL TRY AGAIN")
        

def find_nearest(array, values):
    #tell the function a value, finds closest number to it in an array 
    #and returns index
    indices = np.abs(np.subtract.outer(array, values)).argmin()
    return array[indices]


def convert_wav_to_vel(wavarray,obspeak,labpeak):
	#used for spectra to convert from wavelength to velocity units where 0km/s is set by the peak of the line
    final_vel = []
    
    zpeak = (float(obspeak) - float(labpeak))/float(obspeak)
    
    radvel = 299792 * zpeak
  

    for i in range(0,len(wavarray)):
        
        z = ((wavarray[i] - labpeak)/wavarray[i])
        j = (299792 * z) - radvel
        final_vel.append(j)
    return(final_vel)
  


    
def datafile_2_array(filename,isint=True,zipped=True):
    '''
    need to change to directory and specify datafile beforehand
    returns a list of lists where each list is a column in the datafile
    can then allocate each array after calling this function
    if isint=True then convert rows to ints, otherwise converts to float
    this is specifically for spectra (two/three rows), sometimes with a header
    '''
    
    opening = open(filename, 'r')
    reading = opening.readlines()
    
    obsdata = appendline(reading,0,len(reading))
    
    
    if zipped == False:
        return obsdata
           
    no_rows = (len(obsdata[1])) #just incase theres a header in first line
 
    ###ROWS NEED TO BE ALL SAME LENGTH FOR THIS TO WORK
    lists = [[] for _ in range(no_rows)]
    lists = list(zip(*obsdata))
    
    
    for i in range(0,no_rows):        
        if isint == False:
                    lists[i] = [float(j) for j in lists[i]]
            
        else:            
                    lists[i] = [int(float(j)) for j in lists[i]]
                    
    return lists                



def trim_wav_flux(lambd_points,flux_points,point1,point2):
        #given two points on the x axis of a 1d spectra, 
        #we trim both the x and y axis to be only between these two points
        x = []
        lambd_points = list(lambd_points)
        
        x1 = find_nearest(lambd_points,point1)
        x2 = find_nearest(lambd_points,point2)                
        ind1 = lambd_points.index(x1)
        ind2 = lambd_points.index(x2)
        new_lam = lambd_points[ind1:ind2]
        new_fl = flux_points[ind1:ind2]
        
        return new_lam, new_fl 


def snip_spect(x_axis,flux_axis,*args):
	#args are values you want to snip between, from left to right
	#you're creating a straight line of x = y betweeen the two points 
	#eg. to interpolate between narrow line regions, to snip cosmic rays etc
	#only works on a list
        
        for first, second in zip(args, args[1:]):
               
           if (args.index(first)) % 2 == 0:
                             
               x1 = find_nearest(x_axis,first)              		
               x2 = find_nearest(x_axis,second) 
               ind1 = x_axis.index(x1)       
               ind2 = x_axis.index(x2)
               y1 = flux_axis[ind1]
               y2 = flux_axis[ind2]
               
               arraysize = ind2 - ind1
               
               y_replacementpoints = np.linspace(y1,y2,arraysize)
               flux_axis[ind1:ind2] = y_replacementpoints
               
        return flux_axis
    
#'''    


def replace_str(value,place,linetally):
#USED in conjuction with the fileinput module, which then allows the changes we make to string to be printed back into the file        
        x = linetally.split()       
        x[place] = str(value)
        newx = ' '.join(x)
        string = newx + '\n'
        
        return string

#Replace values in files in input.in fortran file with values defined in this script



def chi_sq(obs_flux,mod_flux,obs_err,mod_err):
    #cutting off last value as damocles returns rebinned model flux array missing final because of how the unequal frequency bins for the model are rebinned
    obs_flux=obs_flux[0:-1]
    
    if len(obs_flux) == len(mod_flux):
        chi_sq=0
        for i in range(len(obs_flux)):         
            chi_sq += ((mod_flux[i] - obs_flux[i])**2/(obs_err**2 + mod_err[i]**2))
    #want to return reduced chi sq so you divide by the array length 
        return chi_sq/len(obs_flux)
    else:
        print(len(obs_flux),len(mod_flux))
        print("observed and modelled flux arrays need to be same length to perform chi sq calculation")
        


def convolve_spectra(res,velarray,fluxarray,Velocity=True):
    
            binwidth = abs(velarray[2] - velarray[1])
           
            sigma = res/(binwidth *2.3548)
            
            #make a gaussian using sigma
            mod_convolve=gaussian_filter1d(fluxarray, sigma)
            
            #mod_convolve = np.convolve(fluxarray,g,boundary = 'extend')
            return(mod_convolve)
             

def make_Grid(v_max,Rrat,rho_index,age_d,divno):
	#creating a gridfile of the supernova here
#where we have 4 1d arrays; a list of x,y,z and density points
#this allows us to plot what the model looks like 

        v_min = v_max * Rrat
        Rout = v_max * age_d * 8.64e-6 * 1e15
        Rin = v_min* age_d * 8.64e-6 * 1e15
      
        #grid divisions for a uniform grid
        grid_arr = np.linspace(-Rout,Rout,divno)

        #These values contain every point in a 40*40*40 grid from the limits of -Rout,Rout 
        Y,X,Z = np.meshgrid(grid_arr,grid_arr,grid_arr) 
        #turning these into 1d arrays
        Xf = X.flatten()
        Yf = Y.flatten()  
        Zf = Z.flatten()
        rads= [np.sqrt(Xf[i]**2 + Yf[i]**2 + Zf[i]**2) for i in range(len(Xf))]
        rho = np.zeros((len(Xf)))
        
        plotdens,plotx,ploty,plotz = [],[],[],[]

	    #looping through radius of every grid point.
	    #if rad is within Rout and Rin, set density using r^-(rho_index) law
        for j in range(len(rads)):
           if abs(rads[j]) <= Rout and abs(rads[j]) >= Rin:
              if Xf[j] > 0 and Yf[j] > 0 and Zf[j] > 0:
                 rho[j] = (abs(rads[j]))**(-rho_index)*1e20   #randomly rescaling the density to make it a reasonable number
              else:
	           
                rho[j] = (abs(rads[j]))**(-rho_index)*1e20
                plotdens.append(rho[j])
                plotx.append(Xf[j])
                ploty.append(Yf[j])
                plotz.append(Zf[j])

        return plotx,ploty,plotz,plotdens



def setax(axis,g_s):
       axis.view_init(elev=30, azim=50)
       axis.set_xlabel('X axis (cm)')
       axis.set_ylabel('Y axis (cm)')
       axis.set_zlabel('Z axis (cm)')
       axis.set_title("Model of gas distribution in a Supernova")
       axis.set_xlim([-1.5*g_s,1.5*g_s])
       axis.set_ylim([-1.5*g_s,1.5*g_s])
       axis.set_zlim([-1.5*g_s,1.5*g_s])

    
   
def initialise_grid_axis(v_max_init,Rrat_init,rho_index_init,age,axis,figur,figur_canv,frame):
        
        #plotting initial grid made by specified initial parameter values and setting up color-bar scale and axis 
        
        x,y,z,d = make_Grid(v_max_init,Rrat_init,rho_index_init,age,20)
        grid_lim= np.amax(x)
        setax(axis,grid_lim)

        l = axis.scatter(x,y,z,c=d,cmap="nipy_spectral")
        cbar = figur.colorbar(l)
        cbar.set_label('density', rotation=270,size=18,labelpad=20)
        cbar.ax.tick_params(labelsize=13)
        figur.tight_layout()
        figur_canv.get_tk_widget().pack(side=tk.TOP, fill='x', expand=1)
       

def update_gasgrid(axis,figur_canv,frame,age,params):
         print("IN GAS GRID",params)
         try: 
             figur_canv.get_tk_widget().pack_forget()
         except AttributeError: 
            pass 
    
         x,y,z,d = make_Grid(params[0],params[1],params[2],age,20) 
         l = axis.scatter(x,y,z,c=d,cmap="nipy_spectral")
               
         figur_canv.get_tk_widget().pack(side=tk.TOP, fill='x', expand=1)
         figur_canv.draw()