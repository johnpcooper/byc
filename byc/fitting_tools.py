import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import UnivariateSpline
from scipy.stats import shapiro
from scipy.stats import gaussian_kde
from scipy.stats import mannwhitneyu
from scipy.stats import linregress
from inspect import signature

from lifelines import KaplanMeierFitter
from lifelines import WeibullFitter

def single_exp(x, a, b, c):
    return a * np.exp(-b * x) + c

def double_exp(x, a, b, c, d, e):
    return a * np.exp(-b * x) + c * np.exp(-d * x) + e

def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def line(x, m, b):
    return m*x + b

def logistic(x, L, k, x_center, offset):
    y = L/(1 + np.exp(-1*k*(x - x_center))) + offset
    return y

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

def exp_fit(cell_df, start_frame, fit_func=single_exp,
            col_name='yfp_norm', background_subtract=True, **kwargs):
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

def get_all_fits_df(dfs_list, start_frame, window_size,
                    fit_func=single_exp, col_name='yfp_norm',
                    background_subtract=True, expl_vars=None,
                     **kwargs):
    
    background = kwargs.get('background', None)
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
        if start_frame is None:
            start_frame = int(dfs_list[i]['chase_frame'].unique()[0])
        try:
            fit_params_dict = exp_fit(dfs_list[i], start_frame,
                                      fit_func=fit_func, col_name=col_name,
                                      background=background,
                                      background_subtract=background_subtract,
                                      expl_var_names=expl_vars)
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
    and retention time in device were observed
    """
    for i in mdf.index:
        end_event_type = mdf.loc[i, 'end_event_type']
        if end_event_type == 'death':
            mdf.loc[i, 'rls_observed'] = True
            mdf.loc[i, 'retention_observed'] = False
        elif end_event_type == 'escape':
            mdf.loc[i, 'rls_observed'] = False
            mdf.loc[i, 'retention_observed'] = True
        elif end_event_type == 'sen':
            mdf.loc[i, 'rls_observed'] = False
            mdf.loc[i, 'retention_observed'] = False
    
    return mdf

def survival_fit(table, **kwargs):
    wbf = WeibullFitter()
    kmf = KaplanMeierFitter()
    T_col = kwargs.get('T_col', 'rls')
    E_col = kwargs.get('E_col', 'observed')
    T = table[T_col]
    E = table[E_col]

    kmf.fit(T, event_observed=E)
    wbf.fit(T, event_observed=E)

    x_kmf = kmf.survival_function_.index
    y_kmf = kmf.survival_function_.values

    x_wbf = np.arange(0, T.max()+5, 1)
    y_wbf = [wbf.survival_function_at_times(x) for x in x_wbf]
    # Find halflife of the kaplan meier and
    # weibull curves
    kmf_halflife = np.nan
    wbf_halflife = np.nan
    kmf_diffs = np.abs(0.5 - np.array(y_kmf))
    # wbf_diffs = np.abs(0.5 - np.array(y_wbf))
    kmf_idx = np.argmin(np.array(kmf_diffs))
    # wbf_idx = np.argmin(np.array(wbf_diffs))
    # wbf_halflife = x_wbf[int(wbf_idx)]
    kmf_halflife = x_kmf[int(kmf_idx)]


    keys = ['y_kmf',
            'x_kmf',
            'x_wbf',
            'y_wbf',
            'kmf', 
            'wbf',
            'kmf_halflife',
            'wbf_halflife']

    data = [y_kmf, 
            x_kmf, 
            x_wbf, 
            y_wbf, 
            kmf, 
            wbf,
            kmf_halflife,
            wbf_halflife]
    
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