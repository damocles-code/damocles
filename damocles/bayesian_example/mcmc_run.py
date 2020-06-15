import os
import emcee
import numpy as np
import pickle
import mcmc_functions as mcmc
import mcmc_set_up

# input for user to vary
resolution = 30 #resolution of the observed data in km/s
n_walkers = 500
check_with_plots = True

# read in all input folders and continuum parameters for each line
# set up all variable parameters with priors and starting distributions
lines, line_summary, init_params, line_profiles = mcmc_set_up.mcmc_set_up(resolution)


# open files to track progress and record output
f = open("chain.dat", "w")
p = open("percent_complete.dat", "w")

# print the current time and date
print(os.system("date"))


# define useful quantities
n_dims = sum(init_params[init_params['variable']==1].variable)
n_lines = len(lines)
assert n_lines > 0, "You have not specified any lines."
assert n_dims > 0, "You must specify at least 1 parameter to vary."


# create mask to pass to fortran
# array is 21 (10 dust params, 11 gas params)
flags = np.zeros((21,6))
for i_line, line in enumerate(lines):
    flags[10:,i_line] = init_params[init_params.line == line].variable.astype('int')
# dust params specified in first column
flags[:10,0] = init_params[init_params.line == 'dust'].variable.astype('int')


# initialise walkers in a ball around best estimate
lower_bounds = np.array(init_params.loc[init_params.variable == 1].prior_min)
upper_bounds = np.array(init_params.loc[init_params.variable == 1].prior_max)
ball_mu = np.array(init_params.loc[init_params.variable == 1].ball_mu)
ball_sd = np.array(init_params.loc[init_params.variable == 1].ball_sd)
pos = emcee.utils.sample_ball(ball_mu, ball_sd, size=n_walkers)

# check initial positions inside prior bounds
for i in range(0, n_walkers):
    while np.isinf(mcmc.ln_uniform_prior(pos[i, :],
                                    lower_bounds,
                                    upper_bounds)):
        pos[i, :] = emcee.utils.sample_ball(ball_mu, ball_sd, size=1)


# set up sampler
sampler = emcee.EnsembleSampler(n_walkers,
                                n_dims,
                                mcmc.ln_prob,
                                args=(lines,
                                      line_summary,
                                      line_profiles,
                                      lower_bounds,
                                      upper_bounds,
                                      flags,
                                      check_with_plots)
                                )

q = open("chi_tracker.dat", "w")
q.close()

# run sampler
nsteps = 50000
chi_min = -10000000000
count = 0
for pos, ln_prob, rstate in sampler.sample(pos, iterations=nsteps):
    count = count + 1
    print("step no", count)
    with open("state.pkl", "w") as filestate:
        pickle.dump([pos, lnprob, rstate], filestate, -1)
    filestate.close()
    p = open("percent_complete.dat", "a")
    p.write("{}".format(count))
    position = pos
    indx = np.argmax(lnprob)
    if lnprob[indx] > chi_min:
        chi_min = lnprob[indx]
        chi_min_pos = position[indx, :]
    q = open("chi_tracker.dat", "a")
    q.write("{:f}".format(chi_min))
    q.write("{:s} \n".format(chi_min_pos))
    q.close()

    # print 'best fit chi', chi_min
    # print 'best fit chi pos', chi_min_pos

    position = position.astype(float)
    f = open("chain.dat", "a")
    i = np.arange(np.shape(position)[0])
    position_w = np.concatenate((i.reshape(np.shape(position)[0], 1), position), axis=1)
    np.savetxt(f, position_w)
    #    for k in range(position.shape[0]):
    #        print k,position[k],type(position[k])
    #        np.savetxt(f,position[k])
    f.close()

    print(
        "Mean acceptance fraction: {0:.3f} ".format(
            np.mean(sampler.acceptance_fraction)
        )
    )
    p.write(
        "Mean acceptance fraction: {0:.3f} ".format(
            np.mean(sampler.acceptance_fraction)
        )
    )

    try:
        print(
            "autocor time: {0:.3f }".format(
                np.mean(sampler.get_autocorr_time(low=5, c=1))
            )
        )
        p.write(
            "autocor time: {0:.3f }".format(
                np.mean(sampler.get_autocorr_time(low=5, c=1))
            )
        )
    except:
        print("cannot determine autocorr time (c=1)")

    try:
        print(
            "autocor time: {0:.3f }".format(
                np.mean(sampler.get_autocorr_time(low=5, c=10))
            )
        )
        p.write(
            "autocor time: {0:.3f }".format(
                np.mean(sampler.get_autocorr_time(low=5, c=10))
            )
        )
    except:
        print("cannot determine autocorr time (c=10)")

    try:
        print("acor time: {0:.3f }".format(np.mean(sampler.acor)))
        p.write("acor time: {0:.3f }".format(np.mean(sampler.acor)))
    except:
        print("cannot determine acor time")

    p.write("\n")
    p.close()
