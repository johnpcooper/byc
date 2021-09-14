import numpy as np
import pandas as pd
# Plotting tools
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
# Stats analysis tools
from scipy.stats import shapiro
from scipy.stats import gaussian_kde
# Stuff from byc
import byc.files as fm

plt_params_dict = {'font.sans-serif': 'Arial'}

def set_styles(plt, matplotlib):
    """
    Set fonts and default style
    """

    try:
        plt.style.use('default')        
        for param, value in plt_params_dict.items():
            matplotlib.rcParams[param] = value
            params = {'mathtext.default': 'regular' }          
            plt.rcParams.update(params)
    except:
        print("""
            Before running set_styles(), you must:

            import matplotlib.pyplot as plt
            import matplotlib
            """)

def make_yticks_0cent(y, **kwargs):
    
    decimals = kwargs.get('decimals', 1)
    pad = kwargs.get('pad', 0.15)
    n_ticks = kwargs.get('n_ticks', 5)

    min_max = np.abs((y.min(), y.max()))
    magnitude = min_max.max()
    y_ul = pad*magnitude + magnitude
    y_ul = np.round(y_ul, decimals)
    y_ll = -y_ul
    
    y_ticks = np.linspace(y_ll, y_ul, n_ticks)
    return y_ticks

def make_ticks(x, **kwargs):
    # can b an x or y. Generate a 5 tick axis starting
    # at 0 by default and ending at x_max + 0.15*x rounded
    # to nearest whole number by default
    
    start_at_zero = kwargs.get('start_at_zero', True) 
    n_ticks = kwargs.get('n_ticks', 5)
    decimals = kwargs.get('decimals', 0)
    pad = kwargs.get('pad', 0.15)
    x_max = np.max(x)
    
    if start_at_zero:
        x_ll = 0
        x_ul = pad*x_max + x_max
        x_ul = np.round(x_ul, decimals)
    else:
        x_range = np.max(x) - np.min(x)
        x_ll = np.min(x) - pad*x_range
        x_ll = np.round(x_ll, decimals)
        x_ul = np.max(x) + pad*x_range
        x_ul = np.round(x_ul)
        
    xticks = np.linspace(x_ll, x_ul, n_ticks)
    
    return xticks

def plot_fits_and_residuals(all_fits_df, dfs_list, expt_name, **kwargs):
    """
    Make and save a figure with 
    1) Scatter plot of decay traces with their y_pred line
       vs time after chase for each cell in all_fits_df
    2) Scatter plot of decay trace residuals vs time after
       chase for each cell in all_fits_df
    3) Kernel density estimate of residuals for each cell in
       all_fits_df as well as one (in black) of the resduals
       across all cells

    Return nothing. This function is ridiculously bloated at the moment
    but quite useful. 
    """

    colors = cm.rainbow(np.linspace(0, 1, len(dfs_list)))
    fig = plt.figure(figsize=(5, 5), tight_layout=True)
    fig.set_dpi(300)

    filename = f'{expt_name}_fits_and_residuals'
    fileformat = '.png'
    
    # Set parameters for decay traces plot
    xlabel_traces = kwargs.get('xlabel_traces', 'Time after Chase (Hrs.)')
    ylabel_traces = kwargs.get('ylabel_traces', 'YFP(t)/YFP(0)')
    ylim_traces = kwargs.get('ylim_traces', (0, 1.2))
    xticks_traces = kwargs.get('xticks_traces', make_ticks(all_fits_df.x_input, decimals=0))
    yticks_traces = kwargs.get('y_ticks_traces', make_ticks((0, 1), decimals=1, n_ticks=7))
    xlim_traces = kwargs.get('xlim_traces', (xticks_traces.min(), xticks_traces.max()))    
    # Set parameters for decay fit residuals plot
    xlabel_resids = kwargs.get('xlabel_resids', xlabel_traces)
    ylabel_resids = kwargs.get('ylabel_resids', 'Residuals')
    xlim_resids = kwargs.get('xlim_resids', xlim_traces)
    xticks_resids = xticks_traces
    yticks_resids = kwargs.get('yticks_resids', make_yticks_0cent(all_fits_df.residual))
    ylim_resids = kwargs.get('ylim_resids', (yticks_resids.min(), yticks_resids.max()))
    
    # Set parameters for decay fit residuals kernel density estimate
    # plot    
    xlabel_kde = kwargs.get('xlabel_kde', ylabel_resids)
    ylabel_kde = kwargs.get('ylabel_kde', 'Density')
    xlim_kde = kwargs.get('xlim_kde', ylim_resids)
    ylim_kde = kwargs.get('ylim_kde', None)
    xticks_kde = yticks_resids
    # yticks_kde will get set below during 
    # density calcuation
    
    # Set parameters used across all plots
    hidden_spines = kwargs.get('hidden_spine', ['top', 'right'])
    labelfontsize = kwargs.get('labelfontsize', 12)
    linewidth = kwargs.get('linewidth', 1)
    linealpha = kwargs.get('linealpha', 1)
    scatteralpha = kwargs.get('scatteralpha', 0.8)
    scattersize = kwargs.get('scattersize', 5)
   
    # Make the residuals scatter plot
    ax = fig.add_subplot(222)
    for cell_index in all_fits_df.cell_index.unique()[:]:
        cell_df = all_fits_df.loc[all_fits_df.cell_index == cell_index, :]
        ax.scatter(cell_df.x_input, cell_df.residual,
                   s=scattersize, alpha=scatteralpha,
                   facecolor='white', edgecolor=colors[cell_index])

    ax.axhline(0, linewidth=linewidth, alpha=linealpha, color='black')
    for spine in [ax.spines[hidden_spine] for hidden_spine in hidden_spines]:
        spine.set_visible(False)
    try:
        ax.set_xticks(xticks_resids)
    except:
        pass
    try:
        ax.set_yticks(yticks_resids)
    except:
        pass
    try:
        ax.set_ylim(ylim_resids)
    except:
        pass
    if xlabel_resids:
        ax.set_xlabel(xlabel_resids, fontsize=labelfontsize)
    if ylabel_resids:
        ax.set_ylabel(ylabel_resids, fontsize=labelfontsize)    

    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')

    # Scatter plot of traces and line plot of fitted decays
    ax2 = fig.add_subplot(221)

    for cell_index in all_fits_df.cell_index.unique()[:]:
        cell_df = all_fits_df.loc[all_fits_df.cell_index == cell_index, :]
        ax2.plot(cell_df.x_input, cell_df.y_pred_norm/cell_df.y_pred_norm.max(),
                linewidth=linewidth, alpha=linealpha, color=colors[cell_index])
        ax2.scatter(cell_df.x_input, cell_df.y_input_norm/cell_df.y_input_norm.max(),
                   s=scattersize, alpha=scatteralpha, color=colors[cell_index])

    if ylim_traces:
        ax2.set_ylim(ylim_traces)
    if xlim_traces:
        ax2.set_xlim(xlim_traces)        
    try:
        ax2.set_xticks(xticks_traces)
    except:
        pass
    try:
        ax2.set_yticks(yticks_traces)
    except:
        pass        
    if xlabel_traces:
        ax2.set_xlabel(xlabel_traces, fontsize=labelfontsize)
    if ylabel_traces:
        ax2.set_ylabel(ylabel_traces, fontsize=labelfontsize)    

    ax2.set_aspect(1.0/ax2.get_data_ratio(), adjustable='box')
    for spine in [ax2.spines[hidden_spine] for hidden_spine in hidden_spines]:
        spine.set_visible(False)

    # Smoothed hist of residuals for each cell (KDE plot)
    ax3 = fig.add_subplot(223)
    
    densities = []
    for cell_index in all_fits_df.cell_index.unique()[:]:
        cell_df = all_fits_df.loc[all_fits_df.cell_index == cell_index, :]
        density = gaussian_kde(cell_df.residual)
        xs = np.linspace(all_fits_df.residual.min(),all_fits_df.residual.max(),200)
        ax3.plot(xs,density(xs), color=colors[cell_index],
                 alpha=linealpha, linewidth=linewidth)
        densities.append(density(xs))

    # Also plot total residuals
    density = gaussian_kde(all_fits_df.residual)
    xs = np.linspace(all_fits_df.residual.min(),all_fits_df.residual.max(),200)
    ax3.plot(xs,density(xs), color='black',
             alpha=linealpha*2, linewidth=linewidth)
    densities.append(density(xs))
    
    # Figure out whech density outuput array has the highest y value and
    # set the yticks of the plot using that density array
    max_dens = np.array([np.max(arr) for arr in densities])
    longest_range_den = densities[max_dens.argmax()]
    yticks_kde = make_ticks(longest_range_den)

    if ylim_kde:
        ax3.set_ylim(ylim)
    if xlim_kde:
        ax3.set_xlim(xlim_kde)
    if xlabel_kde:
        ax3.set_xlabel(xlabel_kde, fontsize=labelfontsize)
    if ylabel_kde:
        ax3.set_ylabel(ylabel_kde, fontsize=labelfontsize)
    try:
        ax3.set_yticks(yticks_kde)
    except:
        pass
    try:
        ax3.set_xticks(xticks_kde)
    except:
        pass

    ax3.set_aspect(1.0/ax3.get_data_ratio(), adjustable='box')
    for spine in [ax3.spines[hidden_spine] for hidden_spine in hidden_spines]:
        spine.set_visible(False)

    if filename and fileformat:
        fig.savefig(f'{filename}{fileformat}', transparent=True)
        print(f'Saved plot at {filename}{fileformat}')

def plot_fluor_trace(cell_df, column_name, subplot_position, **kwargs):
    """
    Create and save a line plot of the column indicated vs. time
    for the whole cell_df 
    """
    # Check for kwargs and set defaults
    x = cell_df.hours
    y = cell_df.loc[:, column_name]
    cell_index = cell_df.cell_index.unique()[0]
    expt_name = kwargs.get('expt_name', cell_df.Label.unique()[0][0:12])
    filetype = kwargs.get('filetype', 'png')
    figsize = kwargs.get('figsize', ())
    filename = f'{expt_name}_cell{cell_index}_{column_name}_trace'
    # Create figure on which to plot traces or get it from kwargs
    fig=kwargs.get('fig', plt.figure(figsize=(4, 1.5), tight_layout=True))
    fig.set_dpi(300)
    linewidth = kwargs.get('linewidth', 0.8)
    tracecolor = kwargs.get('tracecolor', 'red')
    hidden_spines = kwargs.get('hidden_spines', ['top', 'right'])
    xlim = kwargs.get('xlim', (0, np.round(np.max(x), 0)))
    xlabel = kwargs.get('xlabel', 'Time (hrs.)')
    ylim = kwargs.get('ylim', (np.round(np.min(y) - np.max(y)*.1, 0),
                               np.round(np.max(y) + np.max(y)*.1, 0))
                      )
    ylabel = kwargs.get('ylabel', column_name)
    ylabelrotation = kwargs.get('ylabelrotation', 90)
    # Make the plot
    ax = fig.add_subplot(subplot_position)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_yticks(np.linspace(ylim[0], ylim[1], 3))
    ax.set_ylabel(ylabel, rotation=ylabelrotation)
    ax.set_xlabel(xlabel)

    ax.plot(x, y,
            linewidth=linewidth, color=tracecolor, alpha=0.8)

    for spine in [ax.spines[name] for name in hidden_spines]:
        spine.set_visible(False)

    return ax

def plot_cell_peak_detection(cell_df, peak_indices, **kwargs):
    """
    Perform peak detection using byc.trace_tools and plot the result.
    See byc.trace_tools for details on how filtering and detection is 
    achieved
    """    
    cell_index = cell_df.cell_index.unique()[0]
    expt_name = kwargs.get('expt_name', cell_df.Label.unique()[0][0:12])
    filetype = kwargs.get('filetype', 'png')
    filename = f'{expt_name}_cell{cell_index}_manual_vs_auto_bud'
    try:
        manual_bud_indices = kwargs.get('manual_bud_indices')
    except:
        manual_bud_indices = fm.read_roi_position_indices(fm.select_file("Choose manual bud rois .zip"))
    collection_interval = kwargs.get('collection_interval', 10)
    death_cutoff_hr = (np.max(manual_bud_indices)*collection_interval)/60
    linewidth = kwargs.get('linewidth', 0.8)
    hidden_spines = kwargs.get('hidden_spines', ['top', 'right'])
    xlim = kwargs.get('xlim', (0, 70))
    ylim = kwargs.get('ylim', (0.9, 1.1))
    ylim_raw = kwargs.get('ylim_raw', (1000, 5000))
    
    # Create figure on which to plot traces
    fig = plt.figure(figsize=(4, 3), tight_layout=True)
    fig.set_dpi(300)

    # Plot the automaticall discovered bud frames
    ax = fig.add_subplot(212)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_ylabel('Filtered')
    ax.set_xlabel('Time (hrs.)')

    ax.plot(cell_df.hours, cell_df.dsred_mean_local_mean_norm_medfilt,
            linewidth=linewidth, color='black')

    for peak_index in peak_indices:
        peak_hr = (peak_index*collection_interval)/60
        if peak_hr <= death_cutoff_hr:
            ax.axvline(peak_hr,
                       linewidth=linewidth, color='black',
                       linestyle='--')
        else:
            pass
    for spine in [ax.spines[name] for name in hidden_spines]:
        spine.set_visible(False)

    # Plot the manually discovered bud frames
    ax2 = fig.add_subplot(211)
    ax2.set_xlim(xlim)
    ax2.set_ylim(ylim_raw)
    ax2.set_yticks(np.linspace(ylim_raw[0], ylim_raw[1], 3))
    ax2.set_ylabel('Raw DsRed')

    ax2.plot(cell_df.hours, cell_df.dsred_mean,
            linewidth=linewidth, color='red', alpha=0.7)
    for frame in manual_bud_indices[:]:
        frame_hr = (frame*collection_interval)/60
        ax2.axvline(frame_hr,
                   linewidth=linewidth, color='black',
                   linestyle='--')
    for spine in [ax2.spines[name] for name in hidden_spines]:
        spine.set_visible(False)

    fig.savefig(f'{filename}.{filetype}')

def plot_auto_manual_corr(neighbor_df, cell_index, expt_name):
    """
    Really need to clean up defaults etc. on this one but so far
    it works. neighbor_df should be made with:

    neighbor_df = trace_tools.make_bud_neighbor_df(manual_bud_indices, auto_bud_indices=peak_indices)

    """

    collection_interval = 10
    n_auto_buds = len(neighbor_df.auto_bud_frame[~neighbor_df.auto_bud_frame.isnull()])
    n_manual_buds = len(neighbor_df.manual_bud_frame[~neighbor_df.manual_bud_frame.isnull()])
    ylim = (0, 70)
    xlim = ylim
    s=f'Cell {cell_index}\n# of Auto Buds: {n_auto_buds}\n# of Manual Buds: {n_manual_buds}'
    xy=(xlim[1]*0.33, xlim[1]*0.8)

    hidden_spines = ['top', 'right']
    filename = f'{expt_name}_cell{cell_index}_manual_vs_auto_bud_correlation'
    filetype = 'png'

    x = (neighbor_df.auto_bud_frame*collection_interval)/60
    y = (neighbor_df.nearest_manual_frame*collection_interval)/60

    xlabel = 'Auto Bud Hr.'
    ylabel = 'Nearest Manual Bud Hr.'

    fig = plt.figure(figsize=(2.5, 2.5), tight_layout=True)
    fig.set_dpi(300)

    ax = fig.add_subplot()
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.set_xticks(np.linspace(xlim[0], xlim[1], 8))
    ax.set_yticks(np.linspace(xlim[0], xlim[1], 8))
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)

    ax.plot(np.linspace(0, np.max(y), np.max(y)),
            linewidth=1, color='blue', linestyle='--',
            alpha=0.7)

    ax.scatter(x, y,
               s=8, color='white', edgecolor='black', linewidths=0.8)

    for spine in [ax.spines[name] for name in hidden_spines]:
        spine.set_visible(False)

    ax.set_aspect(1.0/ax.get_data_ratio(), adjustable='box')
    ax.annotate(s=s,
                xy=xy,
                fontsize=8)

    fig.savefig(f'{filename}.{filetype}')
 



