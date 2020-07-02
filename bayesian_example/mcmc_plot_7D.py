import numpy as np
import matplotlib.pyplot as plt
import corner
import emcee
from matplotlib.ticker import FormatStrFormatter

#set up parameters
ndim = 7
titles = [r"$v_{\rm max}$","$\beta$",
 r"$v_{\rm min}$",
		"$\beta_{v}$",
r"$\log \ M_{\rm d}$","ff","$\log \ a$"
]
titles_units = [r"$v_{\rm max}$"
				"\n"
				r"$(10^3 {\rm km s}^{-1})$", 
				r"$\beta$",
r"$v_{\rm min}$"
                                "\n"
                                r"$(10^3{\rm km s}^{-1})$",
                                r"$\beta_v$",
				r"$\log \ (M_{\rm d}/M_{\odot})$","ff",r"$\log \ (a/{\rm \mu m})$",#r"$F_{\rm 6300}/F_{\rm 6363}$"#,"$f$",r"$\beta_{\rm \, clump}$"
]
nwalkers = 500
best_fit_vals =[2.4936849,   4.73954609,  0.45378683, -0.47526938, -3.19960678, -0.67382181]
#[ 2.26632542,  7.39682673,  0.45148303, -0.93285957, -2.24742658,  0.11934331]

#read in chain and lop off the ] at the end of each line
data = np.genfromtxt('chain.dat',usecols=(1,2,3,4,5,6,7))

#create chain 
nsteps = int(np.floor(np.shape(data)[0]/nwalkers))
chain = np.empty([nsteps,nwalkers,ndim])

#this is stupid. fix with reshape.
for i in range(nsteps):
    for j in range(nwalkers):
        for k in range(ndim):
            p=nwalkers*i+j
            chain[i,j,k]=data[p,k]

#plot 1 chain for each dimension
for k in range(0,ndim):
   plt.figure()
   for j in range(0,nwalkers):
       plt.plot(chain[:,j,k])
       plt.title(titles[k])

plt.show()

nburnin =50
samples_lin = chain[nburnin:,:, :].reshape((-1, ndim))

print('chain shape', np.shape(chain))

#plot pdfs
samples = chain[nburnin:,:, :].reshape((-1, ndim))
#samples[:,[0,1]] = samples[:,[0,1]]*10

#thin the samples for plotting
samples_thinned = samples[:nwalkers,:]
for i in range(nwalkers,np.shape(samples)[0],nwalkers):
    if (i % (nwalkers*30) == 0):
        samples_thinned = np.concatenate((samples_thinned[:,:],samples[i:i+nwalkers,:]),axis=0)
print('size samples_thinned', np.shape(samples_thinned))

fig1, axes = plt.subplots(ndim, ndim, figsize=(16.5, 16.5))
axes = np.array(fig1.axes).reshape((ndim,ndim))
corner.corner(samples,quantiles=(0.16,0.5,0.84),fig=fig1,show_titles = True,label_kwargs={"fontsize": 20},title_kwargs={"fontsize": 20, "y": 1.06},fontsize = 20)

#for i in range(ndim):
#	print i*ndim
#	fig1.axes[i*ndim].set_ylabel(titles_units[i],fontsize = 20)

#	if i>0:
#		fig1.axes[i*ndim].tick_params(axis='y',labelleft = True)
#	#	fig1.axes[i*6].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#	fig1.axes[ndim*(ndim-1)+i].set_xlabel(titles_units[i],fontsize = 20)
#	fig1.axes[ndim*(ndim-1)+i].tick_params(axis='x',labelbottom = True)
#	fig1.axes[ndim*(ndim-1)+i].xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
#	for j in range(i,ndim-1):
#		fig1.axes[i+(j+1)*5].tick_params(labelsize = 16)
#		
#
#for i in range(ndim):
#        print i*ndim
#        fig1.axes[i*ndim].set_ylabel(titles_units[i])
#        for j in range(i,ndim-1):
#                fig1.axes[i+(j+1)*ndim].tick_params(labelsize = 16)
#                fig1.axes[i+(j+1)*ndim].plot(best_fit_vals[i],best_fit_vals[j+1],marker='o',markerfacecolor='magenta',markeredgecolor='none',linestyle='none')
#                fig1.axes[i+(j+1)*ndim].set_xlabel(titles_units[i])


#fig1.axes[0].set_ylabel('')	
#fig1.axes[0].tick_params(labelsize = 16)
#fig1.axes[(ndim-1)*(ndim-1)].tick_params(labelsize = 16,)

#plt.savefig('05ip_d4075_Ha_Hb_HeI_7D.pdf')

#plot the autocorrelation for each parameter
for i in range(ndim):
   plt.figure()
   g = emcee.autocorr.function_1d(np.mean(chain[:, :, i], axis=1))
   plt.plot(g)
   plt.title(titles[i])

plt.show()
