import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import shapiro
from scipy.stats import gaussian_kde
from inspect import signature

def single_exp(x, a, b, c):
    return a * np.exp(-b * x) + c

def double_exp(x, a, b, c, d, e):
    return a * np.exp(-b * x) + c * np.exp(-d * x) + e

def get_r_squared(y, y_pred):
    """
    Return R-squared as float defined as:

    r_sq = sum of squared model variation from mean / sum of squared data varation from mean

    Note that r-squared is not a valid qualify of fit indicator for non-linear
    regression models. Instead use std. error of the estimate (regression) implemented here
    as get_estimate_std_error
    """
    if (len(y) == len(y_pred)):
        try:
            # Should always work unless y and y_pred aren't np.arrays or
            # pd.DataFrame columns (ie pd.Series)
            y_bar = y.mean()
            r_sq = (np.power(y_pred - y_bar, 2).sum()) / np.power(y - y_bar, 2).sum()
        except:
            print('Failed to calculate r_sq, y and y_pred were not arrays')
            r_sq = False
    else:
        print("Lengths of y and y_bar not the same")
        r_sq = False

    return r_sq

def get_estimate_std_err(y, y_pred):
    """
    Return standard error of the estimate as float defined as:

    std_err = sqrt(sum of squared residuals / n - 2)
    """
    if (len(y) == len(y_pred)):
        try:
            # Should always work unless y and y_pred aren't np.arrays or
            # pd.DataFrame columns (ie pd.Series)
            n = len(y)
            std_err = np.sqrt(np.power(y_pred - y, 2).sum() / (n - 2))
        except:
            print('Failed to calculate std_err, y and y_pred were not arrays')
            std_err = False
    else:
        print("Lengths of y and y_bar not the same")
        std_err = False

    return std_err

def get_shapiro_p(residuals):
    """
    Return shapiro p value of residuals with null hypothesis that
    residuals are normally distributed
    """
    try:
        # Should always work unless y and y_pred aren't np.arrays or
        # pd.DataFrame columns (ie pd.Series)
        shapiro_p = shapiro(residuals)[1]
    except:
        print('Failed to calculate shapiro p, perhaps y and y_pred were not arrays')
        shapiro_p = np.NaN

    return shapiro_p

def exp_fit(cell_df, start_frame, end_frame, fit_func=single_exp, col_name='yfp_norm'):
    """
    Fit a single or double exponential function to the data in cell_df.col_name.
    X axis is cell_df.hours[0:end_frame - startframe]
    """
    param_options = list('abcdefghijklmnop')
    # Check if these explanatory variables are annotated in the cell_df
    expl_var_names = ['cell_index',
                      'age_at_chase',
                      'rls',
                      'div_duration',
                      'dist_from_sen',
                      'late_daughter_shape']
    expl_vars_values = []
    for expl_var_name in expl_var_names:
        try:
            expl_var_value = cell_df.loc[0, expl_var_name]
        except:
            expl_var_value = np.NaN

        expl_vars_values.append(expl_var_value)

    params_dict = dict(zip(expl_var_names, expl_vars_values))

    # If the cell trace isn't as long as end_frame-start_frame, just fit
    # until the end of the trace
    if cell_df.index.max() < end_frame:
        adj_end_frame = cell_df.index.max() + 1
    else:
        adj_end_frame = end_frame + 1

    # Background substract the cell_df column to be fit    
    y_raw = cell_df[col_name][start_frame: adj_end_frame]
    background = y_raw.min()
    y_raw = y_raw - background
    y_raw.index = range(adj_end_frame - start_frame)
    y_norm = y_raw / y_raw[0]
    y_norm.index = range(adj_end_frame-start_frame) # do this so y_norm won't have an index and
    # residuals can be properly calculated from y_norm - y_output_norm
    x = cell_df['hours'][0: adj_end_frame-start_frame]
    
    # Define a list of paramaters to be added to the final params_dict
    names_to_add = ['x_input', 'y_input_raw', 'y_input_norm',
                    'y_pred_norm', 'residual', 'r_sq',
                    'est_std_err', 'shapiro_p']
    n_model_params = len(signature(fit_func).parameters) - 1
    for i in range(n_model_params):
        names_to_add.append(param_options[i])
    # fit the t0 normalized data and add explanatory vars, fit params,
    # and error measurements to params_dict
    try:
        # Fit the data using the fit_func passed this function
        popt, pcov = curve_fit(fit_func, x, y_norm)
        y_pred_norm = fit_func(x, *popt)
        # Define r_sq and standard error of the estimate (est_std_err)
        residuals = y_norm - y_pred_norm
        r_sq = get_r_squared(y=y_norm, y_pred=y_pred_norm)
        est_std_err = get_estimate_std_err(y=y_norm, y_pred=y_pred_norm)
        shapiro_p = get_shapiro_p(residuals)
        # Add fit output and quality measures etc. to params_dict
        values = [x, y_raw, y_norm,
                  y_pred_norm, residuals, r_sq,
                  est_std_err, shapiro_p]
        # Tack on parameters to end of params_dict values lists 
        for param in popt:
            values.append(param)
        # Make sure names_to_add and values are the same length,
        # and add them to the params_dict created above
        if len(names_to_add) == len(values):
            for i in range(len(values)):
                params_dict[names_to_add[i]] = values[i]
        else:
            print(f'names_to_add and values not the same length')
    except:
        print(f"Could not fit cell with cell_index={params_dict['cell_index']}")
        # If the cell can't be fit, then just put an np.NaN in params_dict for
        # every fit-dependent field
        for i in range(len(values)):
            params_dict[names_to_add[i]] = np.NaN

    return params_dict

def get_all_fits_df(dfs_list, start_frame, window_size, fit_func=single_exp, col_name='yfp_norm'):
    if window_size == 'max':
        window_size = np.array([len(df) for df in dfs_list]).max()
    else:
        pass
    end_frame = start_frame + window_size

    fit_params_dicts = []
    fit_params_dfs = []
    for i in range(len(dfs_list)):
        try:
            fit_params_dict = exp_fit(dfs_list[i], start_frame, end_frame, fit_func=fit_func, col_name=col_name)
            fit_params_dfs.append(pd.DataFrame(fit_params_dict))
        except:
            print('fit failed for cell', i)
            fit_params_dict = False
            
        fit_params_dicts.append(fit_params_dict)    
        
         
    all_fits_df = pd.concat(fit_params_dfs, ignore_index=True)
    return all_fits_df


def scan_start_frames(cell_df, col_name='yfp_norm', fit_func=single_exp, window_size=30):
    """ 
    Fit all possible windows of size window_size using the function.
    Return a dictionary called scan_fit_dict with (start_frame, end_frame) 
    of each window as keys and the params_dict output from fit_exp() using
    curve_fit() as the value. 
    """

    fit_results = []
    fit_windows = []

    for timepoint in cell_df.Slice:

        start_frame = timepoint
        end_frame = timepoint + window_size
        # Check if end frame is out of index, if not 
        if start_frame in cell_df.Slice and end_frame in cell_df.Slice:
            #print(f'Found window at ({start_frame}, {end_frame})')
            fit_windows.append(tuple((start_frame-1, end_frame-1)))
            try:
                # should only fail if data cannot be fit with exponential
                params_dict = exp_fit(cell_df,
                                      start_frame-1,
                                      end_frame,
                                      fit_func,
                                      'yfp_norm')
                fit_results.append(pd.DataFrame(params_dict))
            except:
                print(f'fit failed at frame {start_frame}')
                fit_results.append(False)

        else: # This window is out of range of the data
            pass
    # Get all the fits aggregated into a data frame with different fit
    # params vs. start_frame for each window
    scan_fit_dict = dict(zip(fit_windows, fit_results))
    start_frames = [key_tuple[0] for key_tuple in scan_fit_dict.keys()]
    frame_index = 0
    df_index = 0
    fit_result_dfs_list = []
    for fit_result in fit_results:
        try:
            fit_result_row_df = fit_result.groupby(by='cell_index').median().reset_index()
            fit_result_row_df.loc[:, 'start_frame'] = start_frames[frame_index]
            fit_result_dfs_list.append(fit_result_row_df)
            df_index += 1

        except:
            pass
        frame_index += 1
    scan_fit_df = pd.concat(fit_result_dfs_list, ignore_index=True)

    return scan_fit_df