#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 14:58:36 2021

@author: Maria Niculescu-Duvaz
"""
import tkinter as tk
import tkinter.font as TkFont

import os
import damocleslib as model
import numpy as np
import matplotlib
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg,
    NavigationToolbar2Tk
)
import sys
import fileinput
from FUNCTIONS import *
from mpl_toolkits import mplot3d
import numpy as np
from matplotlib.figure import Figure


           
class DamoclesInput(object):
    
    "This class contains all functions and variables which modify input files in the DAMOCLES code from the GUI"
    
    
    def __init__(self):
        
        self.obsfile = InputWindow.start_vars["Data Filename"][1]
        self.z = InputWindow.start_vars["Host Redshift"][1]
        self.wavelength_peak_1 = InputWindow.start_vars["Lab. wavelength peak 1 (A)"][1]
        self.wavelength_peak_2 = InputWindow.start_vars["Lab. wavelength peak 2 (A)"][1]
        self.age_d = int(InputWindow.start_vars["SN age (days)"][1])
        self.snip_regs = InputWindow.start_vars["Snip region range (A)"][1]
        self.bg_lims = InputWindow.start_vars["Continuum region range (A)"][1]
        self.trim_lims = InputWindow.start_vars["Emission Line region range (A)"][1]
        
        self.bg_lims = tuple([float(i) for i in self.bg_lims.split(",")])
        self.trim_lims = tuple([float(i) for i in self.trim_lims.split(",")])
        
        if self.snip_regs == "" :
           self.snip_regs=() 
        else:
          self.snip_regs = tuple([float(i) for i in self.snip_regs.split(",")])
          
         
    
        self.buttonfont = TkFont.Font(family='bitstream charter', size=20)
        
        self.obswav_init,self.obsflux_init= datafile_2_array(self.obsfile,isint=False,zipped=True)
        self.obswav,self.obsflux = trim_wav_flux(self.obswav_init,self.obsflux_init,self.trim_lims[0],self.trim_lims[1])
        self.obsflux = snip_spect(self.obswav,self.obsflux,*self.snip_regs)
        
        self.obsvels = convert_wav_to_vel(self.obswav,(1+self.z)*(self.wavelength_peak_1*10.0),self.wavelength_peak_1*10)
    
  
        self.obs_err = self.get_obserr()
        self.write_obsfile_out()
        
       
        
    def get_obserr(self):
            #calculate observational uncert on input spectrum using background regions provided
                bg_vels,bg_flux = trim_wav_flux(self.obswav_init,self.obsflux_init,self.bg_lims[0],self.bg_lims[1])
                return np.std(bg_flux)

    
    def write_obsfile_out(self):
    #write observed line to line.out file, which is used to rebin the damocles model to the bin number of the observed line
        filey = open(path+"/input/line.in",'w')
        filey.write(str(len(self.obsflux))+ " " + str(self.obs_err) + "\n")
        for j in range(len(self.obsflux)):    
                        filey.write(str(self.obsvels[j]) + ' ' + str(self.obsflux[j]) + "\n")
        filey.close()


    
    def make_clump_button(self,frame):
        
        confirmation = tk.BooleanVar()
        clump_button = tk.Checkbutton(frame, text="Clump?",variable=confirmation,onvalue='True',offvalue='False',command=lambda: self.is_Clump(confirmation),bg='red',
                                      activebackground='red',font=self.buttonfont,borderwidth=5,indicatoron=0)
        clump_button.pack(fill='x')
            

           
    def is_Clump(self,conf): 
     #called when the Clump? button is pressed and changes the input file in damocles
      fi2 = fileinput.FileInput(files=(dust_file),inplace=True)
      
      if conf.get() == True:        
          for line in fi2:	             
            if  'fraction' in line:                
                 line=replace_str('1.0',0,line)
            sys.stdout.write(line)     
            
      else:
         for line in fi2:	
             if  'fraction' in line:
                 line=replace_str('0.0',0,line)         
             sys.stdout.write(line)
      fi2.close()
    
    
    def update_damocfile_input(self,params):
         #writing values we've updated via sliders (where params comes from) to damocles input files 
         fi2 = fileinput.FileInput(files=(dust_file,spec_file),inplace=True)
         for line in fi2:	
             if 'max dust velocity' in line:
                 line=replace_str(params[0],0,line)
             if 'Rin/Rout' in line:
                 line=replace_str(params[1],0,line)
             if 'rho~r^-q' in line:
                 line=replace_str(params[2],0,line)
             if 'Total dust mass' in line:
                 line=replace_str(params[3],0,line)
             if  'dustData' in line:
                 line=replace_str(params[4],3,line)
                 line=replace_str(params[4],4,line)  
             if  'amC' in line:
                  line=replace_str(params[5],2,line)
             if 'sil' in line:
                  line=replace_str(str(float(1-params[5])),2,line)
             sys.stdout.write(line)
         fi2.close()
         
    def initialise_damocfile_input(self):
        
        #values are specified by the user in the opening pane, and in this function are used to 
        #change variables in the damocles input files upon running the GUI
        #so this is only run once
        
        phot_no = InputWindow.start_vars['Photon packet number'][1] 
        is_doublet = InputWindow.start_vars["Doublet?"][1]
        doublet_ratio = InputWindow.start_vars["Doublet ratio"][1]
      
        fi = fileinput.FileInput(files=(input_file,dust_file,gas_file,spec_file),inplace=True)
        
        for line in fi:
            if 'day' in line:
                line=replace_str(self.age_d,0,line)
            if  'number of photons' in line:
                line=replace_str(phot_no,0,line)
            if  'total flux of line to scale model' in line:
                line=replace_str(str(np.amax(self.obsflux)*60),0,line)
            if 'doublet?' in line:
                line=replace_str(is_doublet,0,line)
            if 'Flux ratio' in line:
                line=replace_str(doublet_ratio,0,line)
            if "first doublet component" in line:
                line=replace_str(self.wavelength_peak_1,0,line)
            if "second doublet component" in line:
                line=replace_str(self.wavelength_peak_2,0,line)
       #unless otherwise specified via the interactive button, the default dust distribition is smooth
            if  'fraction in clumps' in line:                
                     line=replace_str('0.0',0,line)
            if  'using observed data' in line:                
                     line=replace_str('True',0,line)           
            sys.stdout.write(line)  

        fi.close()   
        
        



class Plotting_window(DamoclesInput):
    
     "This class contains functions that relate to and define the plotting window of the observed and modelled emission line overplots "
    
     def __init__(self,frame_a_pw,frame_b_pw,frame_c_pw):
        super(Plotting_window,self).__init__()
       
        self.frame_a_pw = frame_a_pw
        self.frame_b_pw = frame_b_pw
        self.frame_c_pw = frame_c_pw
        
        self.buttonfont = TkFont.Font(family='bitstream charter', size=16)
        
        self.fig = Figure(figsize=(7.5, 7.7), dpi=100)
        self.figure_canvas= FigureCanvasTkAgg(self.fig,self.frame_a_pw)
        self.toolbar = NavigationToolbar2Tk(self.figure_canvas, self.frame_a_pw)
        self.toolbar.update()
        
        self.res_kms = InputWindow.start_vars["Resolution (A)"][1]/(self.wavelength_peak_1*10) * 299792 
        
        
     def initialise_chitau_box(self):   
        tk.Label(self.frame_b_pw,text = 'Optical depth (\u03C4) ',font=self.buttonfont).pack(anchor='w',side=tk.LEFT)
        self.tau_text = tk.Text(self.frame_b_pw,height=2,width=20,font=self.buttonfont)
        self.tau_text.pack(anchor='w',side=tk.LEFT)
        tk.Label(self.frame_b_pw,text = '\u03A7^2: ',font=self.buttonfont).pack(side=tk.LEFT)
        self.chi_text = tk.Text(self.frame_b_pw,height=2,width=20,font=self.buttonfont)
        self.chi_text.pack(side=tk.LEFT)
       
        
         
     def initialise_plotwindow(self,frame):
         
          trim_lims_vel = convert_wav_to_vel(self.trim_lims,(1+self.z)*(self.wavelength_peak_1*10.0),self.wavelength_peak_1*10.0)
          
          self.ax = self.fig.add_subplot(111)
          self.ax.axes.set_xlabel("Velocity (km/s)",fontsize=20)
          self.ax.axes.set_ylabel("Flux ($ergs$  $cm^{-2}$  $s^{-1}$  $\AA^{-1}$)",fontsize=20)
          self.ax.tick_params(axis='both', which='major',labelsize=20)
          self.ax.set_xlim([trim_lims_vel[0],trim_lims_vel[1]])
          self.ax.plot(self.obsvels,self.obsflux) 
          self.figure_canvas.get_tk_widget().pack(fill='x', expand=1)

    
     def plot_model(self,frame):
          
          self.modvel,self.modflux,self.modflux_e = datafile_2_array(outfile,isint=False,zipped=True)
          self.modflux = convolve_spectra(self.res_kms, self.modvel,self.modflux)
          self.modflux= [i * np.amax(self.obsflux)/np.amax(self.modflux) for i in self.modflux]
         
          
          self.ax.plot(self.modvel,self.modflux) 
          self.figure_canvas.draw()
            
     def make_reset_button(self,frame):
         
        reset_button = tk.Button(frame, text="Reset",command=lambda: self.clear_pane(frame),bg='orange',
                                 activebackground='grey',font=self.buttonfont,borderwidth=5)
        reset_button.pack(fill='x',side=tk.BOTTOM,anchor='s')
     
     def make_model_scalebox(self):
            labelText=tk.StringVar()
            labelText.set("Scale model by:")
            labelDir= tk.Label(self.frame_c_pw, textvariable=labelText, height=3,font=self.buttonfont)
            labelDir.pack(fill='x',side=tk.LEFT,padx='20',pady='10')
         
            scale_var = tk.StringVar(value="1.0")
            
            def scalebox_command(var):
                sf = float(scale_var.get())
                modflux= [i * sf for i in self.modflux]
                chi = chi_sq(self.obsflux,modflux,self.obs_err,self.modflux_e) 
                chi = round(chi,2)
                self.chi_text.delete(1.0,6.0)
                self.chi_text.insert(tk.END,str(chi))
                
                self.ax.plot(self.modvel,modflux) 
                self.figure_canvas.draw()
            
            scale_var_entry = tk.Entry(self.frame_c_pw, textvariable=scale_var)
            scale_var_entry.bind("<Return>", scalebox_command)   
            scale_var_entry.pack(fill = 'x', side=tk.LEFT, padx='20',pady='10')
  

     def clear_pane(self,frame):
         self.fig.clear() #clear your figure
         self.initialise_plotwindow(frame)
         
     def update_tautextbox(self,frame):
         props = datafile_2_array(mod_prop_file,isint=False,zipped=False)
         for i in props:      
             if 'optical' in i and 'depth' in i and 'wavelength' in i and 'cell' not in i and '(absn)' not in i:
                 tau = i[-1]
         
         self.tau_text.delete(1.0,4.0)
         self.tau_text.insert(tk.END,str(tau))

         
     def update_chitextbox(self,frame):
         
         chi = chi_sq(self.obsflux,self.modflux,self.obs_err,self.modflux_e) 
         chi = round(chi,2)
         
         self.chi_text.delete(1.0,6.0)
         self.chi_text.insert(tk.END,str(chi))
  
    
        
class Slider(Plotting_window):
    def __init__(self,frame_a,frame_b,frame_c,frame_d,frame_e):
        
        self.frame_a=frame_a
        self.frame_b=frame_b
        self.frame_c=frame_c
        self.frame_d=frame_d
        self.frame_e=frame_e
        
        super(Slider,self).__init__(frame_c,frame_d,frame_e)

        
        self.sliderfont = TkFont.Font(family='bitstream charter', size=14)
        self.slider_input_values = {"v_slider": [v_max_init,(1000, 15000),"Vmax (km/s)",1], "r_slider": [Rrat_init,(0.01, 1),"Rin/Rout",0.0005],
                               "rho_slider": [rho_index_init,(-6, 6),"Density index (\u03B2)",0.01], "md_slider": [mdust_init,(-9, 0.2),"log(Dust mass (M\u2609))",0.001],
                               "gs_slider": [grain_size_init,(-2.3, 0.5),'log(Grain radius (\u03BCm))',0.001], "amc_frac_slider": [0.0,(0.0,1.0),'AmC Fraction',0.01]   }
    
        
        #making a widget for every slider 
        for i in self.slider_input_values:
            self.initialise_slider(self.slider_input_values[i])

        
        #initialising all other widgets here, which plot the 
        #as sliders are the only way information interacts with these widgets
        self.fig2 = Figure(figsize=(8.5, 5.3), dpi=100)
        self.ax2 = self.fig2.add_subplot(111, projection='3d')
        self.figure_canvas2= FigureCanvasTkAgg(self.fig2,frame_b)
        initialise_grid_axis(v_max_init,Rrat_init,rho_index_init,self.age_d,self.ax2,self.fig2,self.figure_canvas2,frame_b)
        
        self.initialise_plotwindow(frame_c)
        self.initialise_chitau_box()
        self.make_reset_button(frame_c)
        self.make_model_scalebox()
    
    def initialise_slider(self,slide_params):
        lab = tk.Label(self.frame_a,text=slide_params[2],font=self.sliderfont)
        slider = tk.Scale(self.frame_a,from_=slide_params[1][0],to=slide_params[1][1],orient='horizontal',resolution=slide_params[3],width=17,length=900,font=self.sliderfont)      
        slider.set(slide_params[0])
        slider.bind("<ButtonRelease-1>", self.slider_command)      
        lab.pack(fill='x', padx=1)      
        slider.pack()
        

    def slider_command(self,event):
       
        slider_vals = []
        for child_widget in self.frame_a.winfo_children():
	          if child_widget.winfo_class() == 'Scale':
    	            slider_vals.append(child_widget.get())
        
        #de-log the dust mass and grain size values
        slider_vals[3] = 10**(slider_vals[3])
        slider_vals[4] = 10**(slider_vals[4])
        
        #pass all slider values upon a slider change 
        update_gasgrid(self.ax2,self.figure_canvas2,self.frame_b,self.age_d,slider_vals)
        self.update_damocfile_input(slider_vals)
        #run damocles after sliders have been changed and input files have been changed with slider values
        model.run_damocles_gui_wrap() 
        self.plot_model(self.frame_c)
        self.update_tautextbox(self.frame_d)
        self.update_chitextbox(self.frame_d)


  
class InputWindow(tk.Tk):
    
    #This class sets the "main window" of the application which opens up the interactive damocles modelling environment
    #and into which important variables such as age, redshift of supernova, are set
    
    def __init__(self):
       
       super().__init__()
       
       self.wm_title("Damocles interactive")
       self.filename = tk.StringVar(value='iPTF14hls_2016-11-08_14-31-56_FTN_FLOYDS-N_iPTF_contsub.ascii')
       self.z_var = tk.DoubleVar(value=0.034)
       self.SN_name = tk.StringVar()
       self.Line_name = tk.StringVar()
       self.bg_lims = tk.StringVar(value="5809,6300")
       self.trim_lims = tk.StringVar(value="6650,7050")
       self.snip_reg = tk.StringVar()
       self.is_doublet= tk.StringVar(value="false")
       self.wavelength_peak_1= tk.DoubleVar(value=656.3) 
       self.wavelength_peak_2= tk.DoubleVar(value=0.0) 
       self.doublet_ratio = tk.DoubleVar(value=3.13)
       self.resolution= tk.DoubleVar(value=10.0)
       #age of supernova remnant at the time that observational data was taken, in days.
       self.age_d = tk.DoubleVar(value=778)
       ##no of grid cells in x,y,z direction used to make the spherical shell model in damocles  
       #no of photon packets w. which to run simulation. more packets = more time for each simulation to run and higher SNR model line profile
       self.phot_no = tk.IntVar(value=30000)
       
       

       self.start_vars = {"Data Filename": [self.filename,None],"Host Redshift": [self.z_var,None],"SN name" : [self.SN_name,None],
                          "Line name" : [self.Line_name,None], "Continuum region range (A)":[self.bg_lims,None],
                          "Emission Line region range (A)":[self.trim_lims,None],
                          "Snip region range (A)":[self.snip_reg,None], "Doublet?":[self.is_doublet,None],
                          "Lab. wavelength peak 1 (A)":[self.wavelength_peak_1,None],"Lab. wavelength peak 2 (A)":[self.wavelength_peak_2,None],
                          "Doublet ratio":[self.doublet_ratio,None],"Resolution (A)":[self.resolution,None],
                          "SN age (days)":[self.age_d,None],
                          "Photon packet number":[self.phot_no,None]}
       
       #creating widget for entry boxes interacted with by the user
       count = 0
       for i in list(self.start_vars.keys()):
           self.make_label_entry(i,self.start_vars.get(i)[0],count)
           count += 1
        
      
      # Button to be clicked which opens up modelling app when fields are complete
       tk.Button(self,
                 text='Open modelling app',
                 command=self.open_window,pady=5).grid(columnspan=2)
       
    def make_label_entry(self,labelname,variablename,label_no):
       
       labelText=tk.StringVar()
       labelText.set(labelname + ": ")
       labelDir= tk.Label(self, textvariable=labelText)
       labelDir.grid(column=0,row=label_no,stick=tk.W)
      
       z_entry = tk.Entry(self, textvariable=variablename)
       z_entry.grid(column=1,row=label_no)
     
    
    def open_window(self,event=None):
        #upon pressing button to start app after entry fields have been filled,
        #loop through the dictionary and collect values
            
        for i in list(self.start_vars.keys()):                
                self.start_vars.get(i)[1] = self.start_vars.get(i)[0].get()
      
      
        window = App(self)
        window.grab_set()
      
    
    
class App(tk.Toplevel,Slider):
    def __init__(self,parent):
      
        #this class is the window in which all widgets for the damocles modelling environment are placed
        
        super().__init__(parent)
        

        self.geometry("2000x1500")
        self['bg'] = 'blue'
 
        
        frame_2 = tk.Frame(self)
        frame_1 = tk.Frame(self)
        frame_3 = tk.Frame(self)
        frame_4 = tk.Frame(self)
        frame_5 = tk.Frame(self)
        
        frame_2.pack()
        frame_1.pack()
        frame_3.pack()
        frame_4.pack()
        frame_5.pack()
 
        
        
        self.DamoclesInput = DamoclesInput()
        self.Plotting_window = Plotting_window(frame_3,frame_4,frame_5)
        
        self.DamoclesInput.initialise_damocfile_input()
        self.DamoclesInput.make_clump_button(frame_1)
         
        #create sliders for parameters that are changed by user
        #plotting window things are re-initialised here
        Slider(frame_1,frame_2,frame_3,frame_4,frame_5)

        frame_1.place(x=950,y=515)
        frame_2.place(x=980,y=0) 
        frame_3.place(x=100,y=10)
        frame_4.place(x=105,y=865)
        frame_5.place(x=235,y=920)
        
    


path = os.path.dirname(os.path.realpath(__file__))

input_file = "input/input.in"
gas_file = "input/gas.in"
spec_file = "input/species.in"
dust_file = "input/dust.in"

outfile = path + "/output/integrated_line_profile_binned.out"
mod_prop_file = path + "/output/model_properties.out" 


v_max_init = 4130  #maximum velocity at edge of shell
Rrat_init = 0.31 #radius of inner shell boundary divided by ratio of outer shell. Determines how 'thick' your gas shell is
rho_index_init=1.330   #density is proportional to radius^(-rho_index)
mdust_init=-7.67       #mass of dust in solar masses
grain_size_init=-1.235    #size of dust grain in microns


###########################################################################################
######################## running code  #########################
###########################################################################################


if __name__ == '__main__':
  
  InputWindow = InputWindow()
  InputWindow.mainloop()
