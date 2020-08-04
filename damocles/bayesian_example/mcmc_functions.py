import numpy as np
import damocleslib as model
import matplotlib.pyplot as plt
from astropy.convolution import Gaussian1DKernel
from astropy.convolution import convolve

def convolve_model(model, sd):
    '''

    Parameters
    ----------
    model : array
        modelled line profile fluxes
    sd : float
        standard deviation of gaussian for convolution

    Returns
    -------
    convolved modelled fluxes

    '''

    g = Gaussian1DKernel(stddev=sd)
    mod_convolve = convolve(model, g, boundary="extend")

    return mod_convolve


def ln_like_damocles(theta,
            lines,
            line_summary,
            line_profiles,
            flags,
            check_with_plots
            ):
    '''

    Parameters
    ----------
    theta : array
        the complete parameter list for all lines
    lines : list
        a list of emission line names to be modelled
    line_summary : pandas dataframe
        metadata of the lines
    line_profiles : dictionary of masked arrays
        the observed line profiles with excluded regions masked
    flags : boolean array
        combine with theta to indicate which variables to model
    check_with_plots : boolean
        if true, show plots to check set-up

    Returns
    -------
    chi-squared value of modelled lines (equally weighted)

    '''

    n_lines = len(lines)

    # calculate total number of rows required to store all model line profiles
    n_bins = np.array(line_summary.n_bins)
    n = line_summary.n_bins.sum()

    # run model
    theta_pad = np.zeros((21, 6))
    theta_pad[flags == 1] = theta
    mod = model.run_damocles_wrap(theta_pad, flags, n, n_lines)

    # split model into individual modelled lines
    n_indices = np.cumsum(np.array(n_bins).astype(int))
    mod_multi = np.vsplit(mod, n_indices)
    mod_multi = mod_multi[:-1]

    # initialise arrays to store
    mod_convolved = [[] for i in range(n_lines)]
    mod_chi = [[] for i in range(n_lines)]
    mod_err = [[] for i in range(n_lines)]

    # calculate chi^2 for each line
    for i_line, line in enumerate(lines):

        # calculate flux of observed line
        x_obs = line_profiles[line][:,0].compressed()
        y_obs = line_profiles[line][:,1].compressed()
        y_obs[y_obs<0] = 0
        obs_flux = np.trapz(x_obs, y_obs)

        # convolve the model
        mod_convolved[i_line] = convolve_model(
            (mod_multi[i_line][:, 0]), line_summary.convolution_sigma[i_line]
        )

        # calculate flux of modelled line
        x_mod = line_profiles[line].data[:,0]
        mod_flux = np.trapz(x_mod, mod_convolved[i_line])
        print('min x mod', min(x_mod))
        print('max x mod', max(x_mod))

        # scale the model and uncertainties to the observed flux of that line
        mod_err[i_line] = (mod_multi[i_line][:, 1] * obs_flux / mod_flux)
        mod_convolved[i_line] = (mod_convolved[i_line] * obs_flux / mod_flux)

        # calculate chi^2
        line_mask = np.ma.getmask(line_profiles[line][:,1])
        include_arr = np.logical_not(line_mask)
        mod_chi[i_line] = np.nansum(
            (include_arr / (sum(include_arr)-np.sum(flags)))
            * (
                (mod_convolved[i_line] - np.array(line_profiles[line][:,1])) ** 2
                / ((mod_err[i_line] ** 2 + line_summary.err[line] ** 2))
            )
        )

        print("line no", i_line, "chi model", mod_chi)

        # test plots
        if (check_with_plots):
            plt.figure()
            plt.plot(x_obs, y_obs)
            plt.plot(x_mod, mod_convolved[i_line])
            plt.title(line)
            plt.show()

    # total chi^2 for all lines
    chi2 = sum(mod_chi)
    print("chi", chi2)

    return -chi2


def ln_uniform_prior(theta, lower_bounds, upper_bounds):
    '''


    Parameters
    ----------
    theta : array
        the complete parameter list for all lines
    lower_bounds : array
        lower bounds of uniforms priors
    upper_bounds : array
        upper bounds of uniform priors

    Returns
    -------
    TYPE
        uniform prior return -np.inf if outside bounds

    '''
    inrange = True
    for i in range(len(theta)):
        if lower_bounds[i] < theta[i] < upper_bounds[i]:
            inrange = inrange & True
        else:
            inrange = False
            return -np.inf
            break
    if inrange == True:
        return 0.0
    return -np.inf


def ln_prob(theta,
            lines,
            line_summary,
            line_profiles,
            lower_bounds,
            upper_bounds,
            flags,
            check_with_plots
            ):
    '''

    Parameters
    ----------
    theta : array
        the complete parameter list for all lines
    lines : list
        a list of emission line names to be modelled
    line_summary : pandas dataframe
        metadata of the lines
    line_profiles : dictionary of masked arrays
        the observed line profiles with excluded regions masked
    lower_bounds : array
        lower bounds of uniforms priors
    upper_bounds : array
        upper bounds of uniform priors
    flags : boolean array
        combine with theta to indicate which variables to model
    check_with_plots : boolean
        if true, show plots to check set-up


    Returns
    -------
         value of log posterior as sum of log likelihood and log prior

    '''
    ln_prior = ln_uniform_prior(theta, lower_bounds, upper_bounds)
    
    if np.isinf(ln_prior):
        return -np.inf
    return ln_prior + ln_like_damocles(theta,lines,line_summary,line_profiles,flags,check_with_plots)
