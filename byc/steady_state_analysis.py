import os, re, shutil, skimage

import numpy as np
import pandas as pd
from skimage import io

import tkinter as tk
import tkinter.filedialog as tkdia

from functools import reduce

import nd2

from byc import constants


def select_files(prompt):
    # return a list of files that user selects in the dialogue box
    root = tk.Tk()
    files = tkdia.askopenfilenames(parent=root, title=prompt)
    files_list = root.tk.splitlist(files)
    root.destroy()
    return sorted(files_list)

def select_file():
    # return a list of files that user selects in the dialogue box
    root = tk.Tk()
    file_path = tkdia.askopenfilename(parent=root, title='Choose file')
    root.destroy()
    return file_path

def make_dfs(prompt):
    # return a list of dataframes created by reading the list of filenames
    # that you pass the function
    fns = select_files(prompt)
    dfs = []
    for i in range(0, len(fns)):
        dfs.append(pd.read_csv(fns[i]))
    return dfs, fns

def set_steady_state_dfs_list(master_df, max_n_fovs):    
    """
    Return a list of Dataframes, one for each distinct condition in the dataset
    e.g. plasmid, clone, genetic background
    """
    if max_n_fovs == None:
        max_n_fovs = 20

    dfs_list = []
    for dataset_index in range(0, len(master_df)):
    
        info = master_df.loc[dataset_index, :]
        channel_names = info.fluor_channel_names.split()
        fovs = list(np.arange(0, max_n_fovs+1))

        fov_dfs_list = []
        for fov in fovs:
            channel_dfs = []
            # Info is a series so what would
            # be 'columns' is keys()
            for channel_name in channel_names:
                if ('measdirname' in info.keys() and 'path' in info.keys()):
                    measdirname = info.measdirname
                    if '.tif' in measdirname:
                        measdirname = measdirname.replace('.tif', '')
                    possible_suffs = [f'_{num}' for num in [1, 2, 3, 4, 5]]
                    print(f'Reading from measdirname:\n{measdirname}')
                    # micromanager typically adds a _1 string to the end 
                    # of your base filename in multi-dimensional acquisition
                    # even if it's the first file with that filename. If it
                    # is the second file with that same basefilename micromanager
                    # adds a 2, 3...
                    if measdirname[-2:] in possible_suffs:
                        basefilename = measdirname[0:-2]
                    else:
                        basefilename = measdirname
                    filename = f'{basefilename}_{str(fov).zfill(3)}_{channel_name}.csv'
                    filepath = os.path.join(info.path, filename)
                    print(f'Looking for data at {filepath}')
                    if os.path.exists(filepath):
                        print(f'Found file at path above')
                        channel_df = pd.read_csv(filepath)
                        channel_df = channel_df.rename(columns={'Mean': str(f'{channel_name}_mean'), ' ': 'cell_index'})
                        channel_df = channel_df.rename(columns={'RawIntDen': str(f'{channel_name}_int'), ' ': 'cell_index'})
                        channel_dfs.append(channel_df)
                    else:
                        if len(channel_dfs)>2:
                            pass
                        else:
                            print(f'No file found')
                else:
                    print(f'No measdirname column in master index')
            if len(channel_dfs) > 0:
                fov_merged_df = reduce(lambda x, y: pd.merge(x, y, on='cell_index'), channel_dfs)
                fov_dfs_list.append(fov_merged_df)
        print(f'Found {len(fov_dfs_list)} .csvs')
        if len(fov_dfs_list) > 0:
            final_df = pd.concat(fov_dfs_list, ignore_index=True, sort=False)
        else:
            print('No data found for this condition, returning empty dataframe')
            final_df = pd.DataFrame(None)
            
        # add identifying information to final dataset:
        for i in range(0, len(master_df.columns)):
            column_name = list(master_df.columns)[i]
            value = info[i]
            if not final_df.empty:
                final_df.loc[:, f'{column_name}'] = value
        
        dfs_list.append(final_df)
        
    return dfs_list


def make_expt_df(master_index_path, bg_channel='yfp', filter_cells=False, **kwargs):
    """
    Find a {exptdate}_master_index.csv file at master_index_df_path,
    read in all steady state imaging measurement .csvs found 
    recorded in that master_index as dataframes. Then take some ratios,
    clean, etc. that data and return a concatenated dataframe of all
    those measurement csvs
    """
    master_df = pd.read_csv(master_index_path)
    all_dfs_list = set_steady_state_dfs_list(master_df, max_n_fovs=10)
    all_data_df = pd.concat(all_dfs_list, sort=False, ignore_index=True)
    
    # Set a no-plasmid slice to define background values later
    no_plasmid = all_data_df[all_data_df.plasmid == 'no-plasmid']
    
    # Get the fluor channel names manually recorded in master index
    fluor_channels = master_df.fluor_channel_names.iloc[0].split(' ')
    
    # Set background autofluorescence values for each channel and normalize
    for channel in fluor_channels:
        
        # Normalize signal to median of integrated fluorescence within cell 
        # for no plasmid, BY4741 cells
        channel_bg = no_plasmid[f'{channel}_mean'].median()
        all_data_df.loc[:, f'{channel}_norm'] = all_data_df[f'{channel}_mean']/channel_bg
        
        # Normalize signal to median of integrated fluorescence within cell 
        # for no plasmid, BY4741 cells
        channel_bg = no_plasmid[f'{channel}_int'].median()
        all_data_df.loc[:, f'{channel}_int_norm'] = all_data_df[f'{channel}_int']/channel_bg
        
    # Set coloumns containing ratios of each channel normalized
    # and raw to each other channel
    for channel in fluor_channels:
        
        other_channels = [name for name in fluor_channels if name != channel]
        
        for channel2 in other_channels:
            
            # Ratios of normalized mean within cell channel signals
            ratios = all_data_df[f'{channel}_norm'] / all_data_df[f'{channel2}_norm']
            all_data_df.loc[:, f'{channel}_{channel2}'] = ratios
            # Ratios of raw mean within cell channel signals
            rawratios = all_data_df[f'{channel}_mean'] / all_data_df[f'{channel2}_mean']
            all_data_df.loc[:, f'{channel}_{channel2}_raw'] = rawratios
            # Ratios of normalized integrated within cell channel signals
            ratios = all_data_df[f'{channel}_int_norm'] / all_data_df[f'{channel2}_int_norm']
            all_data_df.loc[:, f'{channel}_{channel2}_int'] = ratios
            # Ratios of raw integrated within cell channel signals
            rawratios = all_data_df[f'{channel}_int'] / all_data_df[f'{channel2}_int']
            all_data_df.loc[:, f'{channel}_{channel2}_int_raw'] = rawratios

    # Filter out background expr cells
    if filter_cells:

        std = all_data_df.loc[all_data_df.plasmid=='no-plasmid', f'{bg_channel}_mean'].std()
        med = all_data_df.loc[all_data_df.plasmid=='no-plasmid', f'{bg_channel}_mean'].median()
        thresh = med + 2*std
        print(f'No-plasmid {bg_channel} median + 2*stdev = {thresh}')
        all_data_df = all_data_df[all_data_df[f'{bg_channel}_mean'] > thresh]
        
    else:
        pass
        
    return all_data_df

def set_flow_cyto_dfs_list(master_df):
    
    """ Return a list of Dataframes, one for each distinct condition in the dataset
        e.g. plasmid, clone, genetic background """
    
    n_datasets_found = 0
    dfs_list = []
    for dataset_index in range(0, len(master_df)):
    
        info = master_df.loc[dataset_index, :]
        channel_names = info.fluor_channel_names.split()
        dataset_id = f'{info.plasmid}_{info.genotype}_{info.clone}'
        
        try:
            filepath = f'{info.path}\\{info.file_name}'
            df = pd.read_csv(filepath)
            print(f"Found .csv for {dataset_id} at {filepath}\n")
            n_datasets_found += 1
            
        except:                
            print(f"No file found at {filename} for dataset {dataset_id}\n")

        # add identifying information to final dataset:
        for i in range(0, len(master_df.columns)):
            column_name = list(master_df.columns)[i]
            value = info[i]
            df.loc[:, f'{column_name}'] = value
        
        dfs_list.append(df)
        
    if len(dfs_list) != len(master_df):
        print("WARNING, did not find a .csv for all rows in master_df")
    else:
        print("Found .csv for all rows in master_df")
        
    return dfs_list

def set_proportional_weights_by_plasmid(df, colname='plasmid'):
    
    """ Return the dataframe passed to this function with a new column called 'weight'. 
        The weight for a row (cell) is 1 - the number of cells with that cell's unique
        plasmid / total number of cells in the data set. 
        
        This allows evenly selecting from each plasmid (or potentially other) group 
        when using df.sample(n=some_number, weights=df.weight). 
        
        WARNING: currently this function applies weight by df.index.levels[0]"""
    
    df.loc[:, 'weight'] = 0
    df = df.set_index([colname, 'cell_index'])

    for plasmid_level in df.index.levels[0]:
        # print(plasmid_level)

        n_cells = len(df.loc[plasmid_level, :])
        # print(f"Number of cells in {plasmid_level} group = {n_cells}")
        proportion = n_cells / len(df)
        # print(f"Fraction of all cell in {plasmid_level} = {proportion}")
        weight = 1 / proportion
        # print(f"weight={weight}")k
        df.loc[plasmid_level, 'weight'] = weight
        
    return df.reset_index()

def set_proportional_weights(df, by=['plasmid', 'genotype', 'clone']):
    
    """ Set the index of the dataframe using the list 
        provided in by. Count the number of cells in 
        each unique group according to that index and set
        a weight value = 1 / proportion of total cells
        in that unique group"""
    
    df.loc[:, 'weight'] = 0
    df = df.set_index(by)

    for group in df.index.unique():
        # print(group)
        n_cells = len(df.loc[group, :])
        # print(f"Number of cells in {group} group = {n_cells}")
        if n_cells != 0:
            proportion = n_cells / len(df)
            # print(f"Fraction of all cells in {group} = {proportion}")
            weight = 1 / proportion
            # print(f"weight={weight}")
            df.loc[group, 'weight'] = weight
        else:
            # print(f'No cells found in {group}')
            df.loc[group, 'weight'] = 0
        
    return df.reset_index()


def crop_image(image):
    """
    Return the central third of the image
    """
    print('Cropping image')
    h = np.round(image.shape[0]/3, 0)
    h = int(h)
    w = np.round(image.shape[1]/3, 0)
    w = int(w)
    crop = image[h:h+300, w:w+300]
    return crop

def read_and_write_rep_images(tif_paths, roi_paths, strain_df, channels, **kwargs):
    """"""
    write_dir = constants.steady_state_rep_image_dir
    contrast = kwargs.get('contrast', False)
    crop = kwargs.get('crop', False)
    # Extract no-pasmid fluorescence level for experiment in 
    # which strain was imaged
    df = strain_df
    bgs = {}
    channel_images = {}
    for channel in channels:
        fov_channel_paths = [path for path in tif_paths if channel in path]
        # Only keep the first FOV for each channel
        tif_path = fov_channel_paths[0]
        image = skimage.io.imread(tif_path)
        tif_write_path = os.path.join(write_dir, os.path.basename(tif_path))
        # If image is brightfield, don't record background
        # intensity value or contrast the image
        if channel != 'bf':
            # Derive background mean fluorescence based on
            # channel_norm = channel_mean/median(channel_mean of no-plasmid cells)
            bg = df.loc[:, f'{channel}_mean']/df.loc[:, f'{channel}_norm']
            bg = bg.unique()[0]
            bgs[channel] = bg
            if contrast:
                image = image/bg
                # Get image back into uint numbers using 100
                # as a unit background fluorescence
                image = np.round(image*100, 0)            
                image = np.array(image, dtype=np.uint8)

        if crop:
            image = crop_image(image)
            skimage.io.imsave(tif_write_path.replace('.tif', '_crop.tif'), image)
        else:
            shutil.copyfile(tif_path, tif_write_path)
        channel_images[channel] = image

        print(f'Saved {channel} at {tif_write_path}')
    return channel_images

def collect_rep_images(plasmid, background, alldf, return_images_rois=False, crop=False, **kwargs):
    contrast = kwargs.get('contrast', False)
    print(f'Finding representative image channels for {plasmid} {background}')
    strainsdf = pd.read_excel(constants.strains_db_path)
    args = (plasmid, background)
    ssdata_dir = constants.steady_state_data_dir
    strainsdf.set_index(['Plasmid', 'Background'], inplace=True)
    strain_expt_dir = strainsdf.loc[args, 'steady_state_image_data'].iloc[0]
    if type(strain_expt_dir) != str:
        print(f'Found non-string expt dir: {strain_expt_dir}')
        return
    elif '|' in strain_expt_dir:
        print(f'Path does not exist {strain_expt_dir}')
        return
    print(f'Looking in {ssdata_dir} for {strain_expt_dir}')
    strain_data_dir = os.path.join(ssdata_dir, strain_expt_dir, 'data')
    
    fns = os.listdir(strain_data_dir)
    clone_number = int(strainsdf.loc[args, 'Clone'])
    strain_pattern = f'(.*)({plasmid})(.*)({background})(.*)(clone{clone_number})(.*)'
    tifs = []
    rois = []
    for fn in fns:
        tif_match = re.search(f'{strain_pattern}(.tif)', fn)
        roi_match = re.search(f'{strain_pattern}(.roi)', fn)
        if tif_match:
            tifs.append(fn)
        if roi_match:
            rois.append(fn)
    tif_paths = [os.path.join(strain_data_dir, t) for t in tifs]
    roi_paths = [os.path.join(strain_data_dir, r) for r in rois]
    # Check if there are any images treated with estradiol
    # found for the strain. If so, only use those ones
    # because they will have induced signal
    estradiol_bools = ['200nM-Estradiol' in path for path in tif_paths]
    estradiol_roi_bools = ['200nM-Estradiol' in path for path in roi_paths]
    if True in estradiol_bools:
        print('Found estradiol treatment')
        tif_paths = np.array(tif_paths)[estradiol_bools]
        roi_paths = np.array(roi_paths)[estradiol_roi_bools]
    # Find channel names in tif files
    channels = []
    for tif in tifs:
        lb = tif.rindex('_')
        ub = tif.rindex('.')
        channels.append(tif[lb+1:ub])
    channels = np.unique(channels)
    # Get imaging measurements for the strain of interest which
    # will be used to find background levels etc.
    strain_df = alldf.set_index(['plasmid', 'genotype']).loc[args, :]
    read_write_args = [tif_paths,
                       roi_paths,
                       strain_df,
                       channels]
    
    channel_images = read_and_write_rep_images(*read_write_args, contrast=contrast, crop=crop)
    if return_images_rois:
        return channel_images, roi_paths
    
def process_timecourse_exptdf(exptdf, mdf):
    channelnames = mdf.fluor_channel_names.loc[0].split()
    bgdict = {'yfp': 190,
              'rfp': 220}
    filter_cells = True
    for channel in channelnames:
        final_tp = exptdf.minutes.max()
        bg = exptdf.loc[exptdf.minutes==final_tp, f'{channel}_mean'].median()
        print(f'Found t180 background for channel {channel} of {bg}')
#         bg = bgdict[channel]
        print(f'Using {channel} background of {bg}')
        exptdf[f'{channel}_norm'] = exptdf[f'{channel}_mean']/bg
        # Define a universal background value and subtract
        # it from the t0 normalized channel column
        bgsub_col = f'{channel}_mean_bgsub'
        exptdf.loc[:, bgsub_col] = exptdf[f'{channel}_mean'] - bg

    if filter_cells:
        exptdf = exptdf[exptdf.rfp_mean_bgsub>0]
    sample_inds = ['plasmid', 'genotype', 'culture_condition']
    exptdf = exptdf.set_index(sample_inds)

    for channel in channelnames:
        t0norm_col = f'{channel}_t0norm'
        channel_col = f'{channel}_mean_bgsub'
        exptdf.loc[:, t0norm_col] = np.nan
        # For each unique sample in the experiment dataframe,
        # normalize fluorescence to the median fluorescence
        # at time 0
        for ind in exptdf.index.unique():
            t0_slice = exptdf[exptdf.minutes==0]
            t0_vals = t0_slice.loc[ind, channel_col]
            t0_fluor_med = t0_vals.median()
            exptdf.loc[ind, t0norm_col] = exptdf.loc[ind, channel_col]/t0_fluor_med

    exptdf.reset_index(inplace=True)
    
    return exptdf

def filter_small_rois(dataframe, **kwargs):
    """
    
    """
    df = dataframe
    xvar = kwargs.get('xvar', 'Area_x')
    yvar = kwargs.get('yvar', 'Feret_x')
    threshold = kwargs.get('threshold', (130, 16))

    keep_mask = (df.Area_x>=threshold[0]) & (df.Feret_x>=threshold[1])
    df = df[keep_mask]

    return df

class Constants():

    def __init__(self):
        # See analysis of 20230216 imaging of JPC121 in YFP and 
        # RFP channels for derivation of this value
        self.yfp_to_rfp_compensation_factor = 0.04019830128347812


def convert_nd2_files(exptname):
    """
    Takes a set of old .nd2 files from the <exptname> directory
    in constants.steady_state_data_dir and saves them as individual
    channel tifs using nd2 metadata
    """
    exptdir = os.path.join(constants.steady_state_data_dir, exptname)
    datasavedir = os.path.join(exptdir, 'data')
    if not os.path.exists(datasavedir):
        os.mkdir(datasavedir)
        print(f'Created directory at\n{datasavedir}')
    filenames = [fn for fn in os.listdir(exptdir) if fn[-4:] == '.nd2']
    metadatas = []
    img_arrs = []
    for filename in filenames:
        filepath = os.path.join(exptdir, filename)
        img_arr = nd2.imread(filepath)
        img_arrs.append(img_arr)
        with nd2.ND2File(filepath) as ndfile:
            metadata = ndfile.metadata
            metadatas.append(metadata)
        # Extract names of channels and their positions in the img_arr
        # then save the image with its channel name as a tif with 
        # skimage.io. These are the names used in updated steady
        # state image analysis
        for channel_idx, Channel in enumerate(metadata.channels):
            channelname = Channel.channel.name
            if channelname == '':
                channelname = 'bf'
            elif channelname == 'YFP':
                channelname = 'yfp'
            elif channelname == 'dsRed':
                channelname = 'rfp'
            else:
                print(f'Using original channel name {channelname}')
            channelfilename = filename.replace('.nd2', f'_{channelname}.tif')
            writepath = os.path.join(datasavedir, channelfilename)
            io.imsave(writepath, img_arr[channel_idx])
            print(f'Saved file at\n{writepath}')