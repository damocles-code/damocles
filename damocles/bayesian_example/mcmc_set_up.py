"""
Set up necessary lists to specify variable parameters and features
"""

import pandas as pd
import numpy as np


def list_param_names():
    """ create list of all possible parameters
    that could be varied for each line """
    global param_list_gas, param_list_dust
    param_list_gas = [
        "gas_max_vel",
        "gas_Rout",
        "gas_radius_ratio",
        "vel_law",
        "gas_density_index",
        "gas_emissivity_index",
        "gas_min_vel",
        "vel_law_r_independent",
        "doublet_flux_ratio",
        "gas_clump_filling_factor",
        "gas_clump_num_density"
    ]

    param_list_dust = [
        "dust_mass",
        "dust_fraction_in_clumps",
        "clump_filling_factor",
        "clump_num_density",
        "dust_max_vel",
        "dust_Rout",
        "dust_radius_ratio",
        "dust_vel_law",
        "dust_density_index",
        "dust_grain_radius"
    ]

def set_up_lines():
    """ determine how many lines are being modelled
    read in names of input directories and continuums for each line """
    global lines, line_summary
    line_list = open("line_list.in", "r").readlines()[4:]
    line_summary = pd.DataFrame([line.split() for line in line_list])
    line_summary.columns = ["line_name", "input_folder", "fixgas", "continuum"]
    lines = list(line_summary["line_name"])
    line_summary.set_index("line_name", inplace=True)
    assert len(lines) > 0, "List of lines to model is empty!"


def read_inputs():
    """ read in priors and starting distributions for all parameters """
    column_names = [
        "variable",
        "prior_min",
        "prior_max",
        "ball_mu",
        "ball_sd",
        "line",
        "parameter_name",
    ]
    global init_params
    init_params = pd.DataFrame(columns=column_names)
    lines_and_dust = ['dust']
    lines_and_dust.extend(lines)

    for line in lines_and_dust:
        if line == 'dust':
            input_file = "bayesian_input_dust.in"
            param_list = param_list_dust
        else:
            input_file = line_summary.input_folder[line] + "/bayesian_input.in"
            param_list = param_list_gas

        line_params = pd.read_csv(input_file, comment="!", header=None, sep="\s+")
        line_params.columns = column_names[:5]
        line_params["line"] = line
        line_params["parameter_name"] = param_list
        init_params = pd.concat([init_params, line_params])
      #re-arranging using a "sort" so that indices so that in multiline case gas geometry parameters are read out to theta variable in the correct order
          
    if len(lines) > 1:
      
        init_params=init_params[init_params.line =='dust'].append(init_params[init_params.line != 'dust'].sort_index())


def read_line_errs():
    """ read in number of bins and line errors from line.in files """
    for line in lines:
        file_name = line_summary.input_folder[line] + "/line.in"
        f = open(file_name, "r")
        n_bins, err = f.readline().split()
        line_summary.loc[line, "n_bins"] = n_bins
        line_summary.loc[line, "err"] = err
    line_summary.loc[:, ["continuum","err","n_bins"]] = line_summary.loc[
        :, ["continuum","err","n_bins"]
    ].astype("int")


def update_line_masks(profile, line):
    """ update line profile mask according to line exclusion regions in file """
    file_name = line_summary.input_folder[line] + '/line_exclusions.in'
    file_lines = open(file_name, 'r').readlines();
    exclude = [line.split() for line in file_lines[1:]];
    for region in exclude:
        profile[:,0] = np.ma.masked_inside(profile[:,0], float(region[0]), float(region[1]));
        profile = np.ma.mask_rowcols(profile, axis=0);
    return profile


def read_line_profiles():
    """ read in all line profiles and subtract continuum """
    global line_profiles
    line_profiles = {}
    for line in lines:
        file_name = line_summary.input_folder[line] + "/line.in"
        profile = np.genfromtxt(file_name, skip_header=1);

        # create profile as masked array with default false
        # i.e. no excluded values (True = exclude)
        profile = np.ma.array(profile, mask=False);
        profile[:, 1] = profile[:, 1] - line_summary.continuum[line];
        profile = update_line_masks(profile, line);


        # create dictionary of masked line profiles
        line_profiles.update({line: profile});


def calculate_sigma(resolution):
    """ calculate sigma to use in convolving models to gaussian """
    for line in lines:
        profile = line_profiles[line].data
        av_bin_width = (np.abs(profile.data[0,0])+np.abs(profile.data[-1,0]))/len(profile)
        sigma = resolution / (av_bin_width * 2.3548)
        line_summary.loc[line,'convolution_sigma'] = sigma


def mcmc_set_up(resolution):
    """ Set up the environment variables including
    environment variables, parameter lists and masked line profile arrays

    Args:
        resolution (float): resolution of the observed data in km/s

    Returns:
        lines (list): list of names of lines to be modelled
        line_summary (dataframe): contains metadata on lines to be modelled
        init_params (dataframe): properties of the parameters to be varied
        line profiles (dataframe): the observed line profiles themselves
    """
    list_param_names()
    set_up_lines()
    read_inputs()
    read_line_errs()

    # print out variable params and starting distributions
    print(
        f"You are modelling {len(lines)} lines: {lines} \n"
        f"with associated input folders {list(line_summary.input_folder)}."
    )
    print()
    print(f"You are varying the following parameters:")
    print(init_params.loc[init_params.variable == 1, "prior_min":])

    read_line_profiles();
    calculate_sigma(resolution)


    return lines, line_summary, init_params, line_profiles


if __name__ == "__main__":
    mcmc_set_up(30)

