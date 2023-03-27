import os

import numpy as np
import pandas as pd

from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.stats import shapiro
from scipy.stats import gaussian_kde
from scipy.stats import mannwhitneyu
from scipy.stats import linregress
from inspect import signature

from lmfit.model import ModelResult
from lmfit import Model
from lmfit import Parameters

from lifelines import KaplanMeierFitter
from lifelines import WeibullFitter

from byc import constants

import matplotlib.pyplot as plt

def single_exp(x, a, b, c):
    return a * np.exp(-b * x) + c

def double_exp(x, a, b, c, d, e):
    return a * np.exp(-b * x) + c * np.exp(-d * x) + e

def line_exp(x, m, intercept, a, k, c):
    y = (m*x + intercept) + (a*np.exp(-k*x) + c)
    return y

def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def line(x, m, b):
    return m*x + b

def logistic(x, L, k, x_center, offset):
    y = L/(1 + np.exp(-1*k*(x - x_center))) + offset
    return y

def sigmoid(x, L ,x0, k, b):
    """
    Standard usage in fitting:

    p0 = [max(y), np.median(x),1,min(y)] # this is an mandatory initial guess
    popt, pcov = curve_fit(fitting_tools.sigmoid, x, y,p0, method='dogbox')
    """
    y = L / (1 + np.exp(-k*(x-x0))) + b
    return (y)

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
            r_sq = 1 - (np.power(y - y_pred, 2).sum()) / np.power(y - y_bar, 2).sum()
        except:
            print('Failed to calculate r_sq, y and y_pred were not arrays')
            r_sq = False
    else:
        print("Lengths of y and y_pred not the same")
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

def get_stderr_by_x(
    x_smooth,
    y_pred,
    sub_fits_df,
    xvar='dist_from_sen',
    yvar='b',
):
    """
    For each unique x value in <x_smooth>, calculate the standard
    error of y_predicted from y data in <sub_fits_df>. If there is
    only one (or zero) y data point at x, set that standard error
    to np.nan

    Return a list of standard errors the same length as 
    <y_pred> and <x_smooth>
    """
    stderrs = []
    sorted_subdf = sub_fits_df.sort_values(by=xvar, ascending=False)
    for i, xval in enumerate(list(x_smooth)):
        ypredval = y_pred[i]
        ydata = sorted_subdf.loc[sorted_subdf[xvar]==xval, yvar]
        if len(ydata) > 0:
            stderr = np.std(ydata)/np.sqrt(len(ydata))
            stderr = np.mean(np.array([np.abs(ypredval - ydataval) for ydataval in ydata]))
            if stderr == 0:
                stderr = np.nan
        else:
            stderr = np.nan
        stderrs.append(stderr)

    return stderrs

def get_r_sq_with_multi_y_per_x(
    x_smooth,
    y_pred,
    sub_fits_df,
    xvar='dist_from_sen',
    yvar='b',
):
    """
    """
    sorted_subdf = sub_fits_df.sort_values(by=xvar, ascending=False)
    residuals_across_all_x = []
    y_pred_values = []
    y_data_values = []
    for i, xval in enumerate(list(x_smooth)):
        yval_pred = y_pred[i]
        ydata = sorted_subdf.loc[sorted_subdf[xvar]==xval, yvar]
        if len(ydata) > 0:
            for value in ydata:
                residual =  yval_pred - value
                residuals_across_all_x.append(residual)
                y_data_values.append(value)
                y_pred_values.append(yval_pred)
        else:
            pass
    y_pred_values = np.array(y_pred_values)
    y_data_values = np.array(y_data_values)
    r_sq = get_r_squared(y_data_values, y_pred_values)
    
    return r_sq

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

def get_standard_exponential_bounds():
    mins = np.array([
        0.99,
        0,
        0
    ])

    maxes = np.array([
        1.01,
        np.inf,
        np.inf
    ])

    bounds = (mins, maxes)

    return bounds

def exp_fit(cell_df, start_frame, fit_func=single_exp,
            col_name='yfp_norm', background_subtract=True,
            transfer_all_cell_df_information=True, bounds=None, **kwargs):
    """
    Fit a single or double exponential function to the data in cell_df.col_name.
    X axis is cell_df.hours[0:end_frame - startframe]
    """
    param_options = list('abcdefghijklmnop')
    # Choose how many total frames to use for fit
    window_width = kwargs.get('window_width', 30)
    background = kwargs.get('background', None)
    # Set limit on how far how many timepoints to include
    # in the background setting (which uses minimum y value
    # in trace as bg)
    background_x_limit = kwargs.get('background_limit', 30)
    x_var = kwargs.get('x_var', 'hours')
    # Check if these explanatory variables are annotated in the cell_df
    if not transfer_all_cell_df_information:
        expl_var_names = kwargs.get('expl_var_names', None)
        if expl_var_names == None:
            expl_var_names = ['cell_index',
                            'age_at_chase',
                            'rls',
                            'div_duration',
                            'dist_from_sen',
                            'late_daughter_shape',
                            'first_bud_frame']
        expl_vars_values = []
        for expl_var_name in expl_var_names:
            try:
                expl_var_value = cell_df[expl_var_name].iloc[0]
            except:
                expl_var_value = np.nan

            expl_vars_values.append(expl_var_value)
            params_dict = dict(zip(expl_var_names, expl_vars_values))
    else:
        params_dict = {}
        expl_var_names = list(cell_df.columns)
        for expl_var_name in expl_var_names:
            params_dict[expl_var_name] = cell_df.loc[cell_df.index[0], expl_var_name]

    # Background substract the cell_df column to be fit
    end_frame = start_frame + window_width  
    y_raw = cell_df.loc[start_frame:end_frame, col_name]
    print(f'Length of y_raw = {len(y_raw)}')
    if background_subtract:
        if background == None:
            background = cell_df[col_name][start_frame:start_frame + background_x_limit].min()
            print(f'Using background val: {background}')
        y_raw = y_raw - background
    else:
        pass
    y_raw.index = range(0, len(y_raw))
    y_norm = y_raw / y_raw.iloc[0]
    print(f'Length of y_raw = {len(y_raw)}. Length of y_norm = {len(y_norm)}')
    x = cell_df.loc[0: len(y_norm)-1, x_var]
    print(f'Fitting with x length {len(x)} and y length {len(y_norm)}')
    
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
        if bounds:
            popt, pcov = curve_fit(fit_func, x, y_norm, bounds=bounds)
        else:
            popt, pcov = curve_fit(fit_func, x, y_norm)
        y_pred_norm = fit_func(x, *popt)
        # Define r_sq and standard error of the estimate (est_std_err)
        print(f'Calculating resids')
        print(f'Length of y_input {len(y_norm)}. Length of y_pred {len(y_pred_norm)}')
        residuals = np.array(y_norm) - np.array(y_pred_norm)
        print(f'calculating R_squared')
        r_sq = get_r_squared(y=y_norm, y_pred=y_pred_norm)
        print(f'Calculating estimated standard error')
        est_std_err = get_estimate_std_err(y=y_norm, y_pred=y_pred_norm)
        print(f'Calculating shapiro_p')
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
    except Exception as e:
        print(f"Could not fit cell with cell_index={params_dict['cell_index']}")
        print(f'Error: {e}')
        # If the cell can't be fit, then just put an np.NaN in params_dict for
        # every fit-dependent field
        for i in range(len(names_to_add)):
            params_dict[names_to_add[i]] = np.NaN

        params_dict['y_input_norm'] = y_norm
        params_dict['y_raw'] = y_raw

    return params_dict

def get_all_fits_df(dfs_list, universal_start_frame, window_size,
                    fit_func=single_exp, col_name='yfp_norm',
                    background_subtract=True, expl_vars=None,
                    bounds=None, manual_bg=None, **kwargs):
    
    background = manual_bg
    if window_size == 'max':
        window_size = np.array([len(df) for df in dfs_list]).max()
    else:
        pass

    fit_params_dicts = []
    fit_params_dfs = []
    for i in range(len(dfs_list)):
        cell_index = dfs_list[i].cell_index.iloc[0]
        print(f'Fitting cell with index {cell_index}')
        # If the user passes None as start_frame, determine start
        # frame from chase_frame colum in cell df
        if universal_start_frame is None:
            start_frame = int(dfs_list[i]['chase_frame'].unique()[0])
            print('Got start frame from trace dataframe')
        else:
            print('Using universal start frame')
        print(f'Using start frame {start_frame} for fit')
        try:
            fit_params_dict = exp_fit(dfs_list[i], start_frame,
                                      fit_func=fit_func, col_name=col_name,
                                      background=background,
                                      background_subtract=background_subtract,
                                      expl_var_names=expl_vars, window_width=window_size,
                                      bounds=bounds)
            fit_params_dfs.append(pd.DataFrame(fit_params_dict))
        except Exception as e:
            print('fit failed for cell ', i)
            print(f'Error: {e}')
            if col_name not in dfs_list[i].columns:
                print(f"Warning, no column with name {col_name} in cell df {i}")
            fit_params_dict = False
            
        fit_params_dicts.append(fit_params_dict)    
        
         
    all_fits_df = pd.concat(fit_params_dfs, ignore_index=True)
    return all_fits_df


def scan_start_frames(cell_df, col_name='yfp_norm', fit_func=single_exp, window_size=20, **kwargs):
    """ 
    Fit all possible windows of size window_size using the function.
    Return a dictionary called scan_fit_dict with (start_frame, end_frame) 
    of each window as keys and the params_dict output from fit_exp() using
    curve_fit() as the value. 
    """
    xvar = kwargs.get('xvar', 'Slice')
    xmax = kwargs.get('xmax', 25)
    xmin = kwargs.get('xmin', 0)
    cell_df.index = range(len(cell_df))
    fit_results = []
    fit_windows = []

    for timepoint in cell_df.loc[xmin:xmax, xvar]:
        print(f'Fitting cell {cell_df.cell_index.iloc[0]} frame {timepoint}')
        start_frame = timepoint
        end_frame = cell_df[xvar].max()
        # end_frame = timepoint + window_size
        # Check if end frame is out of index, if not 
        if start_frame in list(cell_df[xvar]) and end_frame in list(cell_df[xvar]):
            start_frame = int(start_frame)
            print(f'Start frame = {start_frame}')
            end_frame = int(end_frame)
            #print(f'Found window at ({start_frame}, {end_frame})')
            fit_windows.append(tuple((start_frame, end_frame)))
            try:
                # should only fail if data cannot be fit with exponential
                params_dict = exp_fit(cell_df,
                                      start_frame,
                                      fit_func,
                                      col_name=col_name,
                                      x_var=xvar)
                print(f'Attempting to make params_df')
                params_df = pd.DataFrame(params_dict)
                ss_resids = np.sum(np.power(params_df.residual, 2))
                params_df.loc[:, 'ss_residuals'] = ss_resids
                fit_results.append(params_df)
            except Exception as e:
                print(f'fit failed at frame {start_frame}')
                print(e)
                fit_results.append(None)

    # Get all the fits aggregated into a data frame with different fit
    # params vs. start_frame for each window
    scan_fit_dict = dict(zip(fit_windows, fit_results))
    start_frames = [key_tuple[0] for key_tuple in scan_fit_dict.keys()]
    frame_index = 0
    df_index = 0
    fit_result_dfs_list = []
    print(f'Found {len(fit_results)} fit_results')
    for fit_result in fit_results:
        try:
            fit_result_row_df = fit_result.groupby(by='cell_index').median().reset_index()
            fit_result_row_df.loc[:, 'start_frame'] = start_frames[frame_index]
            fit_result_dfs_list.append(fit_result_row_df)
            df_index += 1

        except Exception as E:
            print(E)
            fit_result_dfs_list.append(pd.DataFrame(None))
        frame_index += 1
    scan_fit_df = pd.concat(fit_result_dfs_list, ignore_index=True)

    return scan_fit_df


# Here starts a gen 2 set of fitting functions inspired by how buggy
# fitting_tools.fit_exp() and fitting_tools.get_all_fits_df() etc. are

def params_dict_to_df(params_dict, **kwargs):
    """
    Take a fit parameters dictionary output from 
    fitting_tools.fit_df() and convert it to two 
    pd.DataFrames, one with smoothed x and y prediction
    (longer dataframe) and one without (shorter dataframe)

    Return the dataframe
    """
    pred_keys = ['x_smooth', 'y_pred_smooth']
    other_keys = [key for key in params_dict.keys() if key not in pred_keys]

    pred_df = pd.DataFrame()
    idx = np.arange(0, len(params_dict['x_input']), 1)
    other_df = pd.DataFrame(index=idx)

    for key in pred_keys:
        pred_df.loc[:, key] = params_dict[key]
        
    for key in other_keys:
        # print(f'Adding {key}')
        other_df.loc[:, key] = params_dict[key]

    return other_df, pred_df

def fit_df(df, **kwargs):
    """
    Fit a function (single_exp by default) *xvar* and *yvar*,
    which are column names in *df*, and
    return a dictionary of paramaters found in the fit

    Parameters dict includes many different length
    arrays, from single values to arrays

    If fitting fails for the *df*, (df, pd.DataFrame(None))
    will be returned so that in concatenating multiple samples,
    those for whom fitting failed will have np.nan in place of
    fit parameters and predicted values etc.

    Typically this function is used for fitting flow cytometry data
    where chase start is always t0, rather than byc data where chase
    start can change expt to expt. Need to add option to change
    start frame so I can use it for everything
    """
    fit_func = kwargs.get('fit_func', single_exp)
    xvar = kwargs.get('xvar', 'hours')
    yvar = kwargs.get('yvar', 'yfp_norm')
    # Reset index just in case a multi-indexed dataframe
    # was passed
    if isinstance(df.index, pd.MultiIndex):
        df.reset_index(inplace=True)
    if (xvar in list(df.columns)) and (yvar in list(df.columns)):
        pass
    else:
        print(f'xvar {xvar} and and yvar {yvar} not found in df columns\n{df.columns}')
        return df, pd.DataFrame(None)
    x_input = df[xvar]
    y_input = df[yvar]
    # If fitting fails at this step, return None and 
    # print out the exception. Usually occur because 
    # of bad data or the curve_fit failed to converge
    try:
        popt, pcov = curve_fit(fit_func, x_input, y_input)
    except Exception as e:
        if e == RuntimeError:
            print(f'Runtime isue fitting data using {fit_func}\n')
        else:
            print('Error encountered:\n', e)
        return df, pd.DataFrame(None)
    # Increase density of x points for smoothed
    # y predicted
    x_smooth_min = 0
    x_smooth_max = np.max(x_input)*1.2
    x_smooth_delta = x_smooth_max/50
    x_smooth = np.arange(x_smooth_min,
                         x_smooth_max,
                         x_smooth_delta)
    
    y_pred = fit_func(x_input, *popt)
    y_pred_smooth = fit_func(x_smooth, *popt)
    # Define r_sq and standard error of the estimate (est_std_err)
    residuals = y_input - y_pred
    r_sq = get_r_squared(y=y_input, y_pred=y_pred)
    est_std_err = get_estimate_std_err(y=y_input, y_pred=y_pred)
    shapiro_p = get_shapiro_p(residuals)
    
    dict_keys = ['xvar',
                 'yvar',
                 'x_input',
                 'y_input',
                 'x_smooth',
                 'y_pred_smooth',
                 'y_pred',
                 'residuals',
                 'r_sq',
                 'est_std_err',
                 'shapiro_p']
    
    # Have to set the scope of eval() to local variables
    # only since by default it looks only for global variables
    func_locals = locals()
    dict_vals = [eval(key, func_locals) for key in dict_keys]
    params_dict = dict(zip(dict_keys, dict_vals))
    # Expand popt (params found by curve_fit(x, y)) into
    # single values insteady of the original tuple
    popt_names = list('abcdefghijklmnop')
    for i in range(0, len(popt)):
        name = popt_names[i]
        params_dict[name] = popt[i]
    
    short_df, smooth_df = params_dict_to_df(params_dict)
    # Add back in all the original data in the sampledf
    # to output alogside the fits. Will include a bunch
    # of labels like genotype, chase_method, substrate
    # etc.
    for col in df.columns:
        if col not in list(short_df.columns):
            short_df.loc[:, col] = df.loc[:, col]
    # Now add all the labels to the smooth df as well
    for col in short_df.columns:
        if col not in list(smooth_df.columns):
            smooth_df.loc[:, col] = short_df.loc[:, col].unique()[0]
    return short_df, smooth_df

def fit_experiment(alldf, **kwargs):

    cols = ['sample_id',
            'expt_date']
    cols = kwargs.get('cols', cols)

    sample_fit_dfs = []
    sample_smooth_fit_dfs = []

    indices = alldf.set_index(cols).index.unique()
    for idx in indices:
        print(f'Fitting sample {idx}')
        sampledf = alldf.set_index(cols).loc[idx, :]
        paramsdf, smoothdf = fit_df(sampledf)
        sample_fit_dfs.append(paramsdf)
        sample_smooth_fit_dfs.append(smoothdf)

    allfitsdf = pd.concat(sample_fit_dfs, ignore_index=True)
    allfitsdf_smooth = pd.concat(sample_smooth_fit_dfs, ignore_index=True)

    return allfitsdf, allfitsdf_smooth

def annotate_mdf_censors(mdf):
    """
    To prepare to fit survival curves to RLS data in the 
    master index dataframe (<mdf>), annotate whether death
    and beginning of life were observed
    """
    mdf.loc[mdf.end_event_type=='death', 'rls_observed'] = True
    mdf.loc[mdf.end_event_type=='sen', 'rls_observed'] = False
    mdf.loc[mdf.end_event_type=='escape', 'rls_observed'] = False
    mdf.loc[mdf.observed_since_start==False, 'rls_observed'] = False

    return mdf

def survival_fit(table, **kwargs):
    wbf = WeibullFitter()
    kmf = KaplanMeierFitter()
    T_col = kwargs.get('T_col', 'rls')
    E_col = kwargs.get('E_col', 'rls_observed')
    n_cells = len(table)
    T = table[T_col]
    E = table[E_col]

    kmf.fit(T, event_observed=E)
    wbf.fit(T, event_observed=E)

    x_kmf = kmf.survival_function_.index
    y_kmf = kmf.survival_function_.values

    x_wbf = np.arange(0, T.max()+5, 1)
    y_wbf = [wbf.survival_function_at_times(x) for x in x_wbf]
    # above gives a series for each value so need to cut
    # it down to floats
    y_wbf = [val.values[0] for val in y_wbf]
    # Find halflife of the kaplan meier and
    # weibull curves
    kmf_halflife = np.nan
    wbf_halflife = np.nan
    kmf_halflife = kmf.median_survival_time_
    wbf_halflife = wbf.median_survival_time_

    keys = ['y_kmf',
            'x_kmf',
            'x_wbf',
            'y_wbf',
            'kmf', 
            'wbf',
            'kmf_halflife',
            'wbf_halflife',
            'n_cells',
            'kmf_fit_object',
            'wbf_fit_object']

    data = [y_kmf, 
            x_kmf, 
            x_wbf, 
            y_wbf, 
            kmf, 
            wbf,
            kmf_halflife,
            wbf_halflife,
            n_cells,
            kmf,
            wbf]
    
    fit_dict = dict(zip(keys, data))
    
    return fit_dict

def univariate_spline(df, **kwargs):
    """
    Use scipy.interpolate.UnivariateSpline to create a spline
    estimating the relationship between the <xvar> and <yvar>
    columns of the <df>. By default, fit to the median <yvar>
    at each <xvar> and weight each <xvar> by sample number

    Return xpred, ypred, and the spline object created by UnivariateSpline
    """
    xvar = kwargs.get('xvar', 'dist_from_sen')
    yvar = kwargs.get('yvar', 'b')
    use_sample_weights = kwargs.get('use_sample_weights', True)
    xmax = kwargs.get('xmax', np.nan)

    # Get rid of outlier measurements because they
    # may throw off curve results
    if np.isnan(xmax):
        pass
    else:
        df = df[df[xvar] <= xmax]
    
    medtable = df.pivot_table(index=[xvar], aggfunc=np.median).reset_index()
    counts = df.pivot_table(index=[xvar], aggfunc='count')[yvar]
    counts_norm = counts/np.sum(counts)
    medtable.loc[:, 'weight'] = np.array(counts_norm)

    if use_sample_weights:
        weights = medtable.weight
    else:
        weights = None
    spline = UnivariateSpline(medtable[xvar], medtable[yvar], w=weights)
    xpred = np.arange(0, df.dist_from_sen.max()+1, 1)
    ypred = spline(xpred)

    return xpred, ypred, spline

def df_from_ModelResult(result: ModelResult):
    """
    Extract values, std errors, and inital guesses and ranges
    from the parameters in <result> and return them in as a 
    pd.DataFrame
    """
    fit_params = result.params

    fitdict = {}
    for key in fit_params.keys():
        param = fit_params.get(key)
        stderr = param.stderr
        if stderr:
            std_err_fraction = stderr/param
        else:
            std_err_fraction = None
        fitdict[key] = param.value
        fitdict[f"{key}_stderr"]= stderr
        fitdict[f"{key}_stderr_fraction"] = std_err_fraction

    fitdf = pd.DataFrame(fitdict, index=[0])

    return fitdf

def fit_logistic_to_fits_df(
    sub_fits_df: pd.DataFrame,
    yvar='b',
    xvar='dist_from_sen',
    name=None,
    fitting_func=logistic,
    plot_results=False,
    return_result=False
    ):
    if name is None:
        name = sub_fits_df.strain_name.iloc[0]
    sorted_subdf = sub_fits_df.sort_values(by=xvar, ascending=False)
    # In case you want to fit to median of data per x instead of all data
    table = sorted_subdf.pivot_table(index=xvar, aggfunc=np.median).reset_index()
    df = sorted_subdf

    logistic_model= Model(fitting_func)

    params = Parameters()
    bounds_dict = {
        'L': (-2, 5),
        'k': (-10, 10),
        'x_center': (0, 12),
        'offset': (-2, 5)
    }

    guesses_dict = {
        'L': 2,
        'k': 2,
        'x_center': 3,
        'offset': 1
    }

    for key, val in bounds_dict.items():
        params.add(key, value=guesses_dict[key], min=np.min(val), max=np.max(val))

    result = logistic_model.fit(df[yvar], params, x=df[xvar])

    fitsdf = df_from_ModelResult(result)
    # Calculate young cells decay rate, taken to be the maximum value
    # of decay rate moving away from senescence
    fitsdf.loc[:, 'young_decay_rate'] = fitsdf.loc[0, 'L'] + fitsdf.loc[0, 'offset']
    fitsdf.loc[:, 'young_decay_rate_stderr_fraction'] = fitsdf.loc[0, 'L_stderr_fraction'] + fitsdf.loc[0, 'offset_stderr_fraction']
    fitsdf.loc[:, 'young_decay_rate_stderr'] = fitsdf.young_decay_rate_stderr_fraction*fitsdf.young_decay_rate

    # Offset is the value of the curve at 0 generations from senescence
    fitsdf.loc[:, 'old_decay_rate'] = fitsdf.loc[0, 'young_decay_rate'] - fitsdf.loc[0, 'offset']
    fitsdf.loc[:, 'old_decay_rate_stderr_fraction'] = fitsdf.loc[0, 'offset_stderr_fraction']
    fitsdf.loc[:, 'old_decay_rate_stderr'] = fitsdf.loc[0, 'offset_stderr']
    fitsdf.loc[:, 'strain_name'] = name
    # Create a dataframe that includes standard error at each x and 
    # as well as all the fit params created above
    params = (fitsdf.L.iloc[0], fitsdf.k.iloc[0], fitsdf.x_center.iloc[0], fitsdf.offset.iloc[0])
    x_smooth = range(int(df[xvar].max())+1)
    y_pred = logistic(x_smooth, *params)
    stderrs = get_stderr_by_x(x_smooth, y_pred, sub_fits_df)
    smoothdf = pd.DataFrame(
        {
            'x_input_smooth': x_smooth,
            'y_pred': y_pred,
            'stderr': stderrs
        }
    )
    smoothdf.fillna(method='ffill', inplace=True)
    for column in fitsdf.columns:
        smoothdf.loc[:, column] = fitsdf.loc[0, column]
    # Plot results of the fit
    if plot_results:
        ax = result.plot_fit()
        fig = plt.figure()
        fig.axes.append(ax)
        ax.set_xlim(25, -1)
        ax.set_ylim(0, 4)
    # Save the fit params and error tables
    filename = f'{name}_logistic-fit_b-vs-dist-from-sen.csv'
    filepath = os.path.join(constants.byc_data_dir, f'meta/{filename}')
    filepath = os.path.abspath(filepath)
    smoothdffilepath = filepath.replace('.csv', '_smoothdf.csv')
    figfilepath = filepath.replace('.csv', '.png')
    fitsdf.to_csv(filepath, index=False)
    smoothdf.to_csv(smoothdffilepath, index=False)
    print(f'Saved logistic fit results at\n{filepath}')
    if plot_results:
        fig.savefig(figfilepath)
        print(f'Saved figure at\n{figfilepath}')

    if return_result:
        return fitsdf, smoothdf, result
    else:
        return fitsdf

def fit_line_to_fits_df(
    sub_fits_df: pd.DataFrame,
    yvar='b',
    xvar='dist_from_sen',
    name=None,
    fitting_func=line,
    return_result=False,
    plot_result=False
    ):
    if name is None:
        name = sub_fits_df.strain_name.iloc[0]
    sorted_subdf = sub_fits_df.sort_values(by='dist_from_sen', ascending=False)
    # table = sorted_subdf.pivot_table(index='dist_from_sen', aggfunc=np.median).reset_index()
    # Use df = table if you want to fit to the medians per x value instead of all the data
    # df = table
    df = sorted_subdf

    linear_model = Model(fitting_func)
    result = linear_model.fit(df[yvar], m=0.1, b=1, x=df[xvar])

    fitsdf = df_from_ModelResult(result)
    fitsdf.loc[:, 'strain_name'] = name
    params = (fitsdf['m'].iloc[0], fitsdf['b'].iloc[0])
    x_smooth = range(int(df[xvar].max())+1)
    y_pred = line(x_smooth, *params)
    stderrs = get_stderr_by_x(x_smooth, y_pred, sub_fits_df)
    smoothdf = pd.DataFrame(
        {
            'x_input_smooth': x_smooth,
            'y_pred': y_pred,
            'stderr': stderrs
        }
    )
    smoothdf.fillna(method='ffill', inplace=True)
    for column in fitsdf.columns:
        smoothdf.loc[:, column] = fitsdf.loc[0, column]

    # Plot results of the fit
    if plot_result:
        ax = result.plot_fit()
        fig = plt.figure()
        fig.axes.append(ax)
        ax.set_xlim(25, -1)
        ax.set_ylim(0, 4)
    # Write the fit params and their standard errors as well as 
    filename = f'{name}_linear-fit_b-vs-dist-from-sen.csv'
    filepath = os.path.join(constants.byc_data_dir, f'meta/{filename}')
    filepath = os.path.abspath(filepath)
    figfilepath = filepath.replace('.csv', '.png')
    fitsdf.to_csv(filepath, index=False)
    print(f'Saved logistic fit results at\n{filepath}')
    if plot_result:
        fig.savefig(figfilepath)
        print(f'Saved figure at\n{figfilepath}')

    if return_result:
        return fitsdf, smoothdf, result
    else:
        return fitsdf