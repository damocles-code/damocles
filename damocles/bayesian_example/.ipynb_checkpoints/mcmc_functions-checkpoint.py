import numpy as np
import damocleslib as model
#import matplotlib
#matplotlib.use('agg')
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve, convolve_fft
import sys

#--------------function definitions ---------------
def read_filenames(n_lines,input_folder):

      #open multiline file to read in data file names
      with open('multiline.in') as f:
            f.next
            filenames = [x.split(" ") for x in f.readlines()[1:-3]]
            for i_line in range(n_lines):
                  for i in range(len(filenames[i_line])):
                        filenames[i_line][i] = input_folder[i_line].strip("'")+filenames[i_line][i].strip("'")
      f.close()
      return filenames

def read_profile(multi_line,cont,n_lines,filenames):
      if (multi_line == 1):
            profile = []
            
            for i_line in range(n_lines):
                  data_file = filenames[i_line][1]
                  prof = np.array(np.genfromtxt(data_file,skip_header=1,usecols=(0,1)))
#                 prof[:,1] = np.maximum(prof[:,1]-cont,0)
                  prof[:,1]=prof[:,1]-cont[i_line]
                  profile.append(prof)
            return (profile, n_lines)
      else: 
            profile = np.genfromtxt('input/line.in',skip_header=1)
#           profile[:,1] = np.maximum(profile[:,1]-cont,0)
            profile[:,1] = profile[:,1]-cont
            return profile
            
def read_err(multi_line,n_lines,filenames):
      if (multi_line == 1):
            n_bins = []
            err = []

      #open data file and read first line to store number of bins and array
            for i_line in range(n_lines):
                  data_file = filenames[i_line][1]
                  n_bins.append(np.genfromtxt(data_file, max_rows = 1)[0])
                  err.append(np.genfromtxt(data_file, max_rows = 1)[1])
      else:
            (n_bins, err) = np.genfromtxt('input/line.in', max_rows = 1)
            n_bins = int(n_bins)
      return (n_bins, err)

def read_line_exclusions(profile,multi_line,n_lines,filenames):

      if (multi_line == 1):
      #for every observed line profile, 
      #work out which regions to exclude from chi2calculation
            mask = [[] for i in range(n_lines)]
            for i_line in range(n_lines):
                  data_file = filenames[i_line][2]

                  #read in the limits of the exclusion zones
                  #iterate over each zone to construct mask
                  line_lims = np.genfromtxt(data_file, skip_header = 1)
                  if np.size(line_lims) == 0:
                        #i.e. no exclusion zones
                        mask_single = np.zeros(np.shape(profile[i_line])[0],dtype=bool)
                  elif np.size(line_lims) == 2:
                        mask_single = (profile[i_line][:,0]<line_lims[1]) & (profile[i_line][:,0]>line_lims[0])
                  else:
                        for i in range(np.shape(line_lims)[0]):
                              c = (profile[i_line][:,0]<line_lims[i,1]) & (profile[i_line][:,0]>line_lims[i,0])                     
                              
                              if (i==0):
                                    mask_single = c 
                              mask_single = mask_single | c
                              print i_line,i,mask_single
                  mask[i_line]=(~mask_single).astype(float)
                  print 'mask',i_line,mask
      else:
            #single line case
            line_lims = np.genfromtxt('input/line_exclusions.in', skip_header = 1)
            for i in range(0,np.shape(line_lims)[0]):
                  if np.size(line_lims) == 2:
                        c = (profile[:,0]<line_lims[1]) & (profile[:,0]>line_lims[0])
                        if (i==0):
                              mask = c
                  else:
                        c = (profile[:,0]<line_lims[i,1]) & (profile[:,0]>line_lims[i,0])
                        if (i==0):
                              mask = c
                              mask = mask | c
            mask = (~mask).astype(float)
      return mask

def convolve_model(model,sd):
   g = Gaussian1DKernel(stddev=sd)
   mod_convolve = convolve(model,g,boundary = 'extend')
   return mod_convolve
      
def lnlike(theta,n_bins,profile,err,sigma,multi_line,include_arr,flags,n_lines):
   if (multi_line == 1):
      
         #calculate total number of rows required to store all model line profiles
         n = int(sum(n_bins))

         #run model
         theta_pad = np.zeros((21,6))
         theta_pad[flags==1] = theta
         mod  = model.run_damocles_wrap(theta_pad,  flags, n, multi_line)

         #split model into individual modelled lines
         n_indices = np.cumsum(np.array(n_bins).astype(int))
         mod_multi = np.vsplit(mod,n_indices)
         mod_multi = mod_multi[:-1]


         #initialise arrays to store 
         mod_convolved = [[] for i in range(n_lines)]
         mod_chi = [[] for i in range(n_lines)]

         #calculate chi^2 for each line
         for i_line in range(n_lines):

               #convolve the model
               mod_convolved[i_line] = convolve_model((mod_multi[i_line][:,0]),sigma[i_line])

               #scale the model and uncertainties to the observed flux of that line
               obs_flux = np.trapz(profile[i_line][include_arr[i_line]==1,1], \
                                         profile[i_line][include_arr[i_line]==1,0])

               mod_multi[i_line] = (mod_multi[i_line]*obs_flux/np.trapz( \
                           mod_convolved[i_line][:],profile[i_line][:,0]))
               mod_convolved[i_line] = (mod_convolved[i_line]*obs_flux/np.trapz( \
                           mod_convolved[i_line][:],profile[i_line][:,0]))

               #calculate chi^2
               mod_chi[i_line] = np.nansum((include_arr[i_line]/ \
                                                  sum(include_arr[i_line]))*( \
                           (mod_convolved[i_line] - profile[i_line][:,1])**2/ \
                           ((mod_multi[i_line][:,1]**2 + err[i_line]**2))))

               print 'line no',i_line,'chi model',mod_chi
      
               #test plot
              # plt.figure()
              # plt.plot(profile[i_line][:,0],profile[i_line][:,1])
              # plt.plot(profile[i_line][:,0],mod_convolved[i_line][:])
#               plt.errorbar(profile[i_line][:,0],mod_convolved[i_line][:], yerr=(mod_multi[i_line][:,1]**2 + err[i_line]**2)**0.5)

         #total chi for all lines
         chi2 = sum(mod_chi)

         #plt.show()

   else:
         obs_flux = np.trapz(profile[include_arr[:]==1,1],profile[include_arr[:]==1,0])
         theta_pad = np.zeros((21,6))
         theta_pad[flags==1]=theta
         mod  = model.run_damocles_wrap(theta_pad,  flags, n_bins, multi_line)
#     mod[:,0] = convolve_model(mod[:,0],sigma)
         mod = (mod*obs_flux/np.trapz(mod[:,0],profile[:,0]))
         #      plt.plot(mod[1:-1,0])
         #      plt.plot(profile[1:-1,1])
         #      plt.ylim(0,15000)
         chi2 = np.nansum(include_arr[1:-1]*((mod[1:-1,0]-profile[1:-1,1])**2)/(mod[1:-1,1]**2+(err)**2))
   print 'chi', -chi2
   plt.show()
   return -chi2

def lnprior(theta,lower_bounds,upper_bounds): 
   inrange = True
   for i in range(len(theta)):
      if (lower_bounds[i] < theta[i] < upper_bounds[i]): 
         inrange = inrange & True
      else:
         inrange = False
         return -np.inf
         break
   if inrange == True:
      return 0.0      
   return -np.inf

def lnprob(theta,n_bins,profile,err,sigma,multi_line,lower_bounds,upper_bounds,include_arr,flags,n_lines):
      lp = lnprior(theta,lower_bounds,upper_bounds)
      if not np.isfinite(lp):
            return -np.inf
      return lp + lnlike(theta,n_bins,profile,err,sigma,multi_line,include_arr,flags,n_lines)
#---------------------------------------------------
