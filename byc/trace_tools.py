import numpy as np
import pandas as pd
from scipy.signal import find_peaks, medfilt
from pykdtree.kdtree import KDTree

def local_normalize(df, column_name, kernel_size, name_with_kernel=False):
    assert np.mod(kernel_size, 2) == 1, 'kernel_size must be odd and > 1'    
    
    offset = (kernel_size - 1)/2
    
    new_col_name = f'{column_name}_local_mean_norm'
    
    for i in range(len(df.loc[:, column_name])):
        
        lower = i - offset
        upper = i + offset
        if lower < 0:
            lower = 0
        elif upper >=len(df.loc[:, column_name]):
            upper = i
            
        local_mean = np.mean(df.loc[lower:upper, column_name])
        df.loc[i, new_col_name] = [df.loc[i, column_name]] / local_mean
        
    return df

def mean_filter(df, column_name, kernel_size, name_with_kernel=False):
    assert np.mod(kernel_size, 2) == 1, 'kernel_size must be odd and > 1'    
    
    offset = (kernel_size - 1)/2
    
    new_col_name = f'{column_name}_meanfilt'
    if name_with_kernel:
        new_col_name = f'{new_col_name}_{kernel_size}'
    else:
        pass
    
    for i in range(len(df.loc[:, column_name])):
        
        lower = i - offset
        upper = i + offset
        if lower < 0:
            upper = upper + abs(lower)
            lower = 0

        elif upper >=len(df.loc[:, column_name]):
            upper = i
            
        mean_of_kernel = np.mean(df.loc[lower:upper, column_name])
        df.loc[i, new_col_name] = mean_of_kernel
        
    # return df

def median_filter(df, column_name, kernel_size, name_with_kernel=False):
    assert np.mod(kernel_size, 2) == 1, 'kernel_size must be odd'    
    
    offset = (kernel_size - 1)/2
    
    new_col_name = f'{column_name}_medfilt'
    if name_with_kernel:
        new_col_name = f'{new_col_name}_{kernel_size}'
    else:
        pass
    
    for i in range(len(df.loc[:, column_name])):
        
        lower = i - offset
        upper = i + offset
        if lower < 0:
            upper = upper + abs(lower)
            lower = 0
            
        elif upper >=len(df.loc[:, column_name]):
            upper = i
            
        mean_of_kernel = np.median(df.loc[lower:upper, column_name])
        df.loc[i, new_col_name] = mean_of_kernel
        
    # return df

def add_filtered_columnns(cell_df, column_name, local_norm_kernel_size=63, med_filter_kernel_size=3):
    """
    Return the cell_df datafarme with added column_name_local_mean_norm and
    column_name_local_mean_norm_meanfilt columns added
    """
    filtered_df = local_normalize(cell_df,
                                  column_name=column_name,
                                  kernel_size=local_norm_kernel_size)
    filtered_df = median_filter(filtered_df,
                              column_name=f'{column_name}_local_mean_norm',
                              kernel_size=med_filter_kernel_size)
    return filtered_df

def make_bud_neighbor_df(manual_bud_indices, auto_bud_indices, **kwargs):
    
    """
    auto_bud_indices should be made using find_peaks from scipy.signal
    as follows:
    
    peak_indices, peak_dict = find_peaks(new_df.dsred_mean_local_mean_norm_meanfilt,
                                         distance=min_cycle_frames)
    """
    
    collection_interval = kwargs.get('collection_interval', 10)
    death_cutoff_hr = (np.max(manual_bud_indices)*collection_interval)/60

    df1 = pd.DataFrame({'manual_timepoint': manual_bud_indices[:]})
    df2 = pd.DataFrame({'auto_timepoint': auto_bud_indices[auto_bud_indices <= death_cutoff_hr*6]})

    comparison_df = pd.concat([df1, df2], ignore_index=True, axis=1)
    # Find nearest manual frame for each auto frame
    X = np.array(comparison_df.loc[:, 0])
    Y = np.array(comparison_df.loc[:, 1])    
    tree = KDTree(X)
    nearest_manual_dists, neighbor_indices = tree.query(Y)
    
    nearest_manual_frame = []
    for i in neighbor_indices:
        try:
            nearest_manual_frame.append(X[i])
        except:
            nearest_manual_frame.append(np.NaN)
    # Find nearest auto frame for each manual frame
    X = np.array(comparison_df.loc[:, 0]) # X is manual
    Y = np.array(comparison_df.loc[:, 1]) # Y is auto
    tree = KDTree(Y)
    nearest_auto_dists, neighbor_indices = tree.query(X)
    
    nearest_auto_frame = []
    for i in neighbor_indices:
        try:
            nearest_auto_frame.append(Y[i])
        except:
            nearest_auto_frame.append(np.NaN)

    neighbor_df = pd.DataFrame({'manual_bud_frame': X,
                                'auto_bud_frame': Y,
                                'nearest_manual_frame': nearest_manual_frame,
                                'nearest_auto_frame': nearest_auto_frame,
                                'nearest_auto_dist': nearest_auto_dists,
                                'nearest_manual_dist': nearest_manual_dists})
    
    return neighbor_df


