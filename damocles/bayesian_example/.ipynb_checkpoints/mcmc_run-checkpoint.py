import os
import damocleslib as model
import emcee
import numpy as np
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve, convolve_fft
#import matplotlib.pyplot as plt
import cPickle as pickle
from mcmc_functions import *
import sys
from io import StringIO

#----------------run multiple lines?--------------
multi_line = 1 #set to 1 for true, 0 for false
if (multi_line == 1):
  n_lines = np.fromfile('multiline.in',dtype=int,count=1,sep=' ')
  #read in the names of the input folders from multiline file
  f = open('multiline.in','r')
  lines = f.readlines()[1:n_lines[0]+1]
  input_folder = [[] for i in range(n_lines)] 
  for i in range(n_lines):
    input_folder[i] = lines[i].split(' ')[0]
  f.close()
else:
  n_lines=1
#-------------------------------------------------

#----------------fix continuum---------------------
cont = [2500,300,0.0]
#0.21e-13
#-------------------------------------------------

#--------------initialise--------------------------
#print the current time and date
print os.system('date')

#open file to store chain
f = open("chain.dat", "w")
f.close()
p = open('percent_complete.dat','w')
p.close()
#----------------------------------------------------

#--------read in variable flags and priors---------
flags=np.zeros((21,6))
prior_min=[[] for i in range(n_lines)]
prior_max=[[] for i in range(n_lines)]
ball_mu=[[] for i in range(n_lines)]
ball_sd=[[] for i in range(n_lines)]

print 'input folder', input_folder[0]

if (multi_line == 1):
  for i_line in range(n_lines):
    input = np.genfromtxt(input_folder[i_line].strip("'") + 'bayesian_input.in',skip_header=9,usecols=(0,1,2,3,4))
    flags[:,i_line] = input[:,0]

  #initialise priors and initial starting ball
    prior_min[i_line] = []*int(sum(flags[:,i_line]))
    prior_max[i_line] = []*int(sum(flags[:,i_line])) 
    ball_mu[i_line] = []*int(sum(flags[:,i_line]))
    ball_sd[i_line] = []*int(sum(flags[:,i_line]))

#assign priors
    prior_min[i_line][:] = input[flags[:,i_line]==1,1]
    prior_max[i_line][:] = input[flags[:,i_line]==1,2]
    ball_mu[i_line][:] = input[flags[:,i_line]==1,3]
    ball_sd[i_line][:] = input[flags[:,i_line]==1,4]

print 'input',input
print 'prior',prior_min

if (multi_line == 1):
  n_dims=sum(flags)
else:
  ndim = 6
#----------------------------------------------------

#-------------read in data and set up arrays---------
#read in all relevant files names from multiline file
if (multi_line == 1):
  filenames = read_filenames(n_lines,input_folder)
else:
  filenames = []

#read in observed line profile for chi2 calculation
(n_bins, err) = read_err(multi_line,n_lines,filenames)
if (multi_line == 1):
  (profile, n_lines) = read_profile(multi_line,cont,n_lines,filenames)
else:
  profile = read_profile(multi_line,cont,n_lines,filenames)

#read in and set up include/exclude array
include_arr = read_line_exclusions(profile,multi_line,n_lines,filenames)

#calculate sigma to use in convolution to gaussian
resolution_vel = 30
if (multi_line == 1):
  sigma=np.zeros(n_lines)
  for i_line in range(n_lines):
    av_bin_width = np.mean([profile[i_line][i,0]-profile[i_line][i-1,0] for i in range(1,np.shape(profile[i_line])[0])])
    sigma[i_line]=resolution_vel/(av_bin_width*2.3548)
else:
  av_bin_width = np.mean([profile[n,0]-profile[n-1,0] for n in range(1,np.shape(profile)[0])])
  sigma = resolution_vel/(av_bin_width*2.3548)
#-----------------------------------------------------

print 'prior', prior_min[0]
print 'prior', prior_min[1]
print 'prior', prior_min[2]

#initialise walkers in a ball around best estimate
if (multi_line==1):
  ndim = int(sum(n_dims))
lower_bounds=prior_min[0]
upper_bounds=prior_max[0]
ball_mean = ball_mu[0]
ball_sigma = ball_sd[0]
for i_line in range(1,n_lines):
  lower_bounds = lower_bounds + prior_min[i_line]
  upper_bounds = upper_bounds + prior_max[i_line]
  ball_mean = ball_mean + ball_mu[i_line]
  ball_sigma = ball_sigma + ball_sd[i_line]
print (ball_mean)
print (ball_sd)
print (lower_bounds)
print (upper_bounds)


nwalkers = 500
pos=emcee.utils.sample_ball(ball_mean, ball_sigma, size=nwalkers)
for i in range(0,nwalkers):
    while np.isinf(lnprior(pos[i,:],lower_bounds,upper_bounds)):
        pos[i,:]=emcee.utils.sample_ball(ball_mean, ball_sigma, size=1)
    print pos[i,:],i

#set up sampler
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args = (n_bins, profile, err, sigma,multi_line,lower_bounds,upper_bounds,include_arr,flags,n_lines))

q=open('chi_tracker.dat','w')
q.close()

#run sampler
nsteps = 50000
chi_min = -10000000000
count=0
for pos, lnprob, rstate in sampler.sample(pos, iterations=nsteps):
   count=count+1
   print 'step no',count
   with open("state.pkl", "w") as filestate:
              pickle.dump([pos, lnprob, rstate], filestate, -1)
   filestate.close()
   p = open('percent_complete.dat','a')
   p.write('{}'.format(count))
   position = pos
   indx = np.argmax(lnprob)
   if (lnprob[indx]>chi_min):
      chi_min = lnprob[indx]
      chi_min_pos = position[indx,:]
   q = open('chi_tracker.dat','a')
   q.write('{:f}'.format(chi_min))
   q.write('{:s} \n'.format(chi_min_pos))
   q.close()

   #print 'best fit chi', chi_min
   #print 'best fit chi pos', chi_min_pos

   position = position.astype(float)
   f = open("chain.dat", "a")
   i = np.arange(np.shape(position)[0])
   position_w = np.concatenate((i.reshape(np.shape(position)[0],1),position),axis=1)
   np.savetxt(f,position_w)
#    for k in range(position.shape[0]):
#        print k,position[k],type(position[k])
#        np.savetxt(f,position[k])
   f.close()

   print("Mean acceptance fraction: {0:.3f} ".format(np.mean(sampler.acceptance_fraction)))
   p.write("Mean acceptance fraction: {0:.3f} ".format(np.mean(sampler.acceptance_fraction)))

   try:
      print("autocor time: {0:.3f }".format(np.mean(sampler.get_autocorr_time(low=5,c=1))))
      p.write("autocor time: {0:.3f }".format(np.mean(sampler.get_autocorr_time(low=5,c=1))))
   except:
      print 'cannot determine autocorr time (c=1)'

   try:
      print("autocor time: {0:.3f }".format(np.mean(sampler.get_autocorr_time(low=5,c=10))))
      p.write("autocor time: {0:.3f }".format(np.mean(sampler.get_autocorr_time(low=5,c=10))))
   except:
      print 'cannot determine autocorr time (c=10)'

   try:
      print('acor time: {0:.3f }'.format(np.mean(sampler.acor)))
      p.write('acor time: {0:.3f }'.format(np.mean(sampler.acor)))
   except:
      print 'cannot determine acor time'

   p.write('\n')
   p.close()
