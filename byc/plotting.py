import os
from sys import maxunicode

from time import time

import numpy as np
import pandas as pd

# image tools
import cv2
import skimage
from skimage import io
# Plotting tools
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
# Stats analysis tools
from scipy.stats import shapiro
from scipy.stats import gaussian_kde
# Stuff from byc
from byc import files, constants, segmentation

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

def savefig_with_timestamp(fig, tight_layout=True):
    if tight_layout:
        plt.tight_layout()
    abspath = os.path.abspath(os.path.join('.', f'{time()}.png'))
    fig.savefig(abspath)
    print(f'Figure saved at \n{abspath}')

def fig_ax(figsize=(2.5, 2.5), subplotpos=(111)):
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(subplotpos)
    fig.set_dpi(300)
    remove_spines(ax)
    return fig, ax

class annoying_strings():
    """
    A place to store strings that I can never remember how to format
    """
    def __init__(self):

        self.k_inverse_hrs = 'k (hr$^{-1}$)'
        self.r_sq = 'R$^2$'
        self.mu = u"\u03bc"

def figure(**kwargs):
    figsize = kwargs.get('figsize', (2.5, 2.5))
    height_scale = kwargs.get('height_scale', 1)
    width_scale = kwargs.get('width_scale', 1)
    dpi = kwargs.get('dpi', 250)
    
    figsize = (figsize[0]*width_scale, figsize[1]*height_scale)
    fig = plt.figure(figsize=figsize)
    fig.set_dpi(dpi)

    return fig

def figure_ax(**kwargs):
    figsize = kwargs.get('figsize', (2.5, 2.5))
    height_scale = kwargs.get('height_scale', 1)
    width_scale = kwargs.get('width_scale', 1)
    dpi = kwargs.get('dpi', 250)
    n_subplots = kwargs.get('n_subplots', 1)
    
    figsize = (figsize[0]*width_scale, figsize[1]*height_scale)
    fig = plt.figure(figsize=figsize)
    fig.set_dpi(dpi)
    
    axes = []
    for subplot_index in range(n_subplots):
        ax = fig.add_subplot(1, n_subplots, subplot_index+1)
        axes.append(ax)

    remove_spines(ax)
    if n_subplots == 1:
        return fig, ax
    else:
        print(f'Returning fig, {n_subplots} axes as list')
        return fig, axes

def legend_outside():
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2)

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

def format_ticks(ax, **kwargs):
    """
    take <ax>, turn on minor ticks and set all ticks
    to inner direction. Can pass <xminorspace> and
    <yminorspace> as keyword args. Otherwise defaults
    to a minor tick every 1/5 of space between major 
    ticks.

    Return nothing
    """
    add_minor_x = kwargs.get('add_minor_x', True)
    add_minor_y = kwargs.get('add_minor_y', True)
    tickdirection = kwargs.get('tickdirection', 'in')      
    # Default to minor tick every 1/5 of space between
    # major ticks
    try:
        xmajorspace = np.diff(ax.get_xticks())[0]
    except:
        xmajorspace = np.nan
    ymajorspace = np.diff(ax.get_yticks())[0]
    xminorspace = np.abs(xmajorspace/2)
    yminorspace = np.abs(ymajorspace/2)

    yminorspace = kwargs.get('yminorspace', yminorspace)
    xminorspace = kwargs.get('xminorspace', xminorspace)

    if add_minor_x:
        ax.xaxis.set_minor_locator(MultipleLocator(xminorspace))
    if add_minor_y:
        ax.yaxis.set_minor_locator(MultipleLocator(yminorspace))

    for ticktype in ['minor', 'major']:
        ax.tick_params(axis="y", which=ticktype, direction=tickdirection)
        ax.tick_params(axis="x", which=ticktype, direction=tickdirection)

def transparent_boxes(ax, **kwargs):
    """
    Return nothing. Set box plot fill to white, edgecolor
    to black by default

    Sourced from https://stackoverflow.com/questions/43434020/black-and-white-boxplots-in-seaborn
    """
    facecolor = kwargs.get('facecolor', 'white')
    edgecolor = kwargs.get('edgecolor', 'black')
    edgewidth = kwargs.get('edgewidth', 1)
    for i,box in enumerate(ax.artists):
        box.set_edgecolor(edgecolor)
        box.set_facecolor(facecolor)
        box.set_linewidth(edgewidth)
            # iterate over whiskers and median lines
        for j in range(6*i,6*(i+1)):
            try:
                ax.lines[j].set_color(edgecolor)
                ax.lines[j].set_linewidth(edgewidth)
            except Exception as E:
                # print(f"Could not change color of line with j {j}")
                pass

def make_x_axis(**kwargs):
    """
    Plot and save an x axis with no other ticks
    or spines. Can pass path at which to save.

    Typically used to create an axis for a color map
    """
    xmax = kwargs.get('xmax', 40)
    xmin = kwargs.get('xmin', 0)
    stepsize = kwargs.get('stepsize', 5)
    xlabel = kwargs.get('xlabel', 'Generations from\nSenescence')
    figsize = kwargs.get('figsize', (2.5, 2.5))
    fontsize = kwargs.get('fontsize', 12)
    savepath = kwargs.get('savepath', f'legend_xlim_{xmin}-{xmax}.svg')
    # Create the figure and plot
    fig = plt.figure(figsize=figsize)
    fig.set_dpi(300)
    ax = fig.add_subplot(111)
    ax.set_xlim(xmin, xmax)
    ax.set_xticks(np.arange(xmin, xmax+1, stepsize))
    # Remove spines
    for spine in [ax.spines[key] for key in ['top', 'right', 'left']]:
        spine.set_visible(False)
    ax.set_xlabel(xlabel, fontsize=fontsize)
    # No y ticks!
    ax.set_yticks([])
    fig.savefig(savepath)

def remove_spines(ax):
    """
    Set top and right spines on <ax> to invisible
    """
    for spine in [ax.spines[name] for name in ['top', 'right']]:
        spine.set_visible(False)


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
        manual_bud_indices = files.read_roi_position_indices(files.select_file("Choose manual bud rois .zip"))
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
 
def plot_survival(survival_fit_dict,
                  ax,
                  color,
                  linestyle,
                  linewidth,
                  alpha,
                  size):
    
    x_kmf = survival_fit_dict['x_kmf']
    y_kmf = survival_fit_dict['y_kmf']
    x_wbf = survival_fit_dict['x_wbf']
    y_wbf = [df.iloc[0]for df in survival_fit_dict['y_wbf']]
    
    ax.scatter(x_kmf, y_kmf, color='white', edgecolor=color,
               alpha=alpha, s=size)
#     ax.plot(x_kmf, y_kmf, color=color,
#                alpha=alpha, linewidth=linewidth)
    ax.plot(x_wbf, y_wbf, color=color,
               alpha=alpha, linestyle=linestyle, linewidth=linewidth)


def plot_radial_avg_intensity_peaks(cell_index, crop_rois_df, stacksdict,
                                    **kwargs):

    savefig = kwargs.get('savefig', True)
    dpi = kwargs.get('dpi', 100)
    ylim = kwargs.get('ylim', (1500, 2500))
    
    celldf = crop_rois_df[crop_rois_df.cell_index==cell_index]
    cellstack = stacksdict[str(int(cell_index))]
    save_paths = []
    for frame_idx in celldf.frame_rel:
        cellframedf = celldf[celldf.frame_rel==frame_idx]
        peak_radius = cellframedf.peak_avg_intensity_radius.iloc[0]
        frame_idx = int(frame_idx)
        mpl = matplotlib
        max_distance = 30
        # Adjust plot parameters
        mpl.rcParams['font.family'] = 'Arial'
        mpl.rcParams['font.size'] = 16
        mpl.rcParams['axes.linewidth'] = 2
        mpl.rcParams['axes.spines.top'] = False
        mpl.rcParams['axes.spines.right'] = False
        mpl.rcParams['xtick.major.size'] = 7
        mpl.rcParams['xtick.major.width'] = 2
        mpl.rcParams['ytick.major.size'] = 7
        mpl.rcParams['ytick.major.width'] = 2
        # Create figure and add subplot
        fig = plt.figure(figsize=(5, 5))
        fig.set_dpi(dpi)
        ax = fig.add_subplot(211)
        # Plot data
        rad = cellframedf.circle_slice_radius.iloc[0].split('|')
        rad = np.array([float(r) for r in rad])
        intensity = cellframedf.intensity_avg.iloc[0].split('|')
        intensity = np.array([float(inte) for inte in intensity])
        ax.plot(rad, intensity, linewidth=2, color='red')
        # Edit axis labels
        ax.set_xlabel('Radial Distance (pixels)', labelpad=10)
        ax.set_ylabel('Average Intensity (AU)', labelpad=10)
        ax.set_xlim(0, max_distance)
        ax.set_ylim(ylim)
        ax.set_title(f'Cell {cell_index}, frame {frame_idx}')
        intensity_slice = intensity[0:max_distance+1]
        intensity_slice = intensity_slice - intensity_slice.min()
        # peaks[0] is the index at which the peak is found. 
        # So index + 1 = actual distance
        ax.axvline(x=peak_radius+1, color='black', linestyle='--')
        ax = fig.add_subplot(212)

        ax.imshow(cellstack[frame_idx])
        x_center_rel = celldf.x_center_rel.iloc[frame_idx]
        y_center_rel = celldf.y_center_rel.iloc[frame_idx]
        circ = mpl.patches.Circle(xy=(x_center_rel, y_center_rel), radius=peak_radius,
                                  color='white', fill=False)
        ax.add_patch(circ)
        plt.tight_layout()

        if savefig:
            reldir = celldf.compartment_reldir.iloc[0]
            exptname = f'{celldf.date.iloc[0]}_{celldf.expt_type.iloc[0]}'
            save_dir = os.path.join(constants.byc_data_dir, reldir)
            save_dir = f'{save_dir}_auto'
            save_dir = os.path.join(save_dir, 'plots')
            if not os.path.exists(save_dir):
                os.mkdir(save_dir)
            filename = f'{exptname}_cell{str(cell_index).zfill(3)}_frame{str(frame_idx).zfill(3)}.tif'
            save_path = os.path.join(save_dir, filename)
            fig.savefig(save_path)
            print(f'Saved figure at\n{save_path}')
            save_paths.append(save_path)
            # Close the frame figure so it stops taking up memory
            plt.close()
    # Read the individual frames back in and concatanate
    # them into a stack and save
    images = [skimage.io.imread(p) for p in save_paths]
    stack = skimage.io.concatenate_images(images)
    rdx = save_path.rindex(f'cell{str(cell_index).zfill(3)}')
    stack_save_path = f'{save_path[0:rdx+7]}_rad_avg_stack.tif'
    skimage.io.imsave(stack_save_path, stack)
    print(f'Saved full stack at {stack_save_path}')

# for cell_index in mdf.cell_index.unique():
def plot_cell_radial_segmentation(allframesdf, stacksdict, **kwargs):

    manual_savedir = kwargs.get('manual_savedir', None)
    tracemax = allframesdf.intensity.max() 
    tracemin = allframesdf.intensity.min()
    r = tracemax - tracemin
    buff = 0.1*r
    ylim = (tracemin - buff, tracemax + buff)
    cell_index = allframesdf.cell_index.unique()[0]
    savefig = True
    dpi = 100
    xlim = (0, 20)
    celldf = allframesdf
    if type(list(stacksdict.keys())[0]) == str:
        cellstack = stacksdict[str(int(cell_index))]
    else:
        cellstack = stacksdict[cell_index]

    # Adjust plot parameters
    mpl = matplotlib
    mpl.rcParams['font.family'] = 'Arial'
    mpl.rcParams['font.size'] = 16
    mpl.rcParams['axes.linewidth'] = 2
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['xtick.major.size'] = 7
    mpl.rcParams['xtick.major.width'] = 2
    mpl.rcParams['ytick.major.size'] = 7
    mpl.rcParams['ytick.major.width'] = 2

    save_paths = []
    for frame_idx in celldf.frame_rel.unique():
        cellframedf = celldf[celldf.frame_rel==frame_idx]
        frame_idx = int(frame_idx)
        max_distance = 30
        palette = sns.color_palette("icefire",
                                    n_colors=len(cellframedf.theta.unique()))
        # Create figure and add subplot
        fig = plt.figure(figsize=(5, 8))
        fig.set_dpi(dpi)
        ax = fig.add_subplot(311)
        ax.set_title(f'Cell {cell_index}, frame {frame_idx}')
        # Plot radial intensity traces by theta. Sometimes this doesn't work because
        # of mismatch between palette length and actual number of theta angles.
        # Maybe beacuse of NaNs?
        try:
            sns.lineplot(x='radial_distance', y='intensity', hue='theta', data=cellframedf,
                         hue_order=cellframedf.theta.unique(), palette=palette, ax=ax)
        except Exception as e:
            print(f'Could not plot radial intensity curves for cell {cell_index}')
            print(f'Exception:\n{e}')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_ylabel('Intensity', fontsize=16)
        ax.set_xlabel('Radial Distance from Center (px)', fontsize=16)
        # ax.annotate(f'Frame {frame_idx}', )
        plt.legend(bbox_to_anchor=(1.05, 1), loc=2)
        ax.legend_.remove()
        # Plot cell crop frame with scatter plot of intensity peaks
        ax = fig.add_subplot(312)
        cellcrop = cellstack[frame_idx]
        ax.imshow(cellcrop)
        for i, t in enumerate(list(cellframedf.theta.unique())):
            c = palette[i]
            x = cellframedf[cellframedf.theta==t].peak_X.iloc[0]
            y = cellframedf[cellframedf.theta==t].peak_Y.iloc[0]
            ax.scatter(x, y, color=c)
        x_center_rel = cellframedf.x_center_rel.iloc[0]
        y_center_rel = cellframedf.y_center_rel.iloc[0]
        ax.scatter(x_center_rel, y_center_rel, color='white')

        # Plot cell outline mask
        frame_mask = segmentation.get_frame_cell_mask(allframesdf, cellstack, frame_idx)
        ax = fig.add_subplot(313)
        ax.imshow(frame_mask)
        plt.tight_layout()
        if savefig:
            if manual_savedir == None:
                reldir = celldf.compartment_reldir.iloc[0]
                save_dir = os.path.join(constants.byc_data_dir, reldir)
            else:
                save_dir = manual_savedir
                
            exptname = f'{celldf.date.iloc[0]}_{celldf.expt_type.iloc[0]}'
            save_dir = os.path.join(save_dir, 'plots')
            if not os.path.exists(save_dir):
                os.mkdir(save_dir)
            filename = f'{exptname}_cell{str(cell_index).zfill(3)}_frame{str(frame_idx).zfill(3)}_rad_theta_peaks.tif'
            save_path = os.path.join(save_dir, filename)
            fig.savefig(save_path)
            print(f'Saved figure at\n{save_path}')
            save_paths.append(save_path)
            # Close the frame figure so it stops taking up memory
            plt.close()
    # Read the individual frames back in and concatanate
    # them into a stack and save
    images = [skimage.io.imread(p) for p in save_paths]
    stack = skimage.io.concatenate_images(images)
    rdx = save_path.rindex(f'cell{str(cell_index).zfill(3)}')
    stack_save_path = f'{save_path[0:rdx+7]}_rad_theta_bin_stack.tif'
    skimage.io.imsave(stack_save_path, stack)
    for path in save_paths:
        os.remove(path)
    print(f'Saved full stack at {stack_save_path}')


def plot_cell_chase_stack(tracedf, cellstacksdict, cellkey, manual_contrast=False, **kwargs):
    """
    Make a tif stack plot of YFP vs. time for the data in <tracedf>
    with inset axes showing raw data for the bf and yfp channels of
    the cell contained in <cellstacksdict>

    <cellkey> should be compartment_name-cell0 or some other index
    used to uniquely identify individual trace measurements in a dataset
    """
    
    time_delta_hrs = kwargs.get('time_delta_hrs', 10/60)
    tracedf.sort_values(by='x_input', inplace=True)

    channel_stacks = []
    missing_channels = []
    for key in cellstacksdict.keys():
        channel_stack = cellstacksdict[key]
        channel_stacks.append(channel_stack)
        if type(channel_stack) == np.ndarray:
            channel = key
        else:
            missing_channels.append(key)
            print(f'{cellkey} missing {key} channel')
    stack = cellstacksdict[channel]
    frame = stack[0]
    height = frame.shape[0]
    width = frame.shape[1]
    width_scale = width/width
    height_scale = height/width
    mins = [np.min(frame) for frame in list(stack)]
    meds = [np.median(frame) for frame in list(stack)]
    maxes = [np.max(frame) for frame in list(stack)]

    basename = f'{cellkey}_chase_{channel}'
    savedir = f'{constants.byc_data_dir}meta\\plots'
    if not os.path.exists(savedir):
        os.mkdir(savedir)
    x, y = tracedf.x_input, tracedf.y_input_norm
    y_pred = tracedf.y_pred_norm
    b = tracedf.b.iloc[0]
    dist_from_sen = tracedf.dist_from_sen.iloc[0]
    chase_frame = int(tracedf.chase_frame.iloc[0])

    xlim = (0, 3)
    ylim = (0, 1.2)

    xlabel = 'Hours'
    ylabel = 'YFP/YFP(t=0)'

    xticks = np.arange(np.min(xlim), np.max(xlim) + 0.1, 0.5)

    vmin = np.median(meds[-10:])
    vmax = np.max(maxes)

    if not manual_contrast:
        vmin = None
        vmax = None

    linekwargs = {
        'color': 'black'
    }

    scatterkwargs = {
        's': 15,
        'facecolor': 'white',
        'edgecolor': 'black'
    }

    imkwargs = {
        'cmap': 'gray',
        'vmin': vmin,
        'vmax': vmax
    }
    savepaths = []
    for i, frame in enumerate(list(stack)):

        fig, ax = figure_ax()
        fig.set_dpi(150)
        ax.set_xlim(xlim)
        ax.set_xticks(xticks)
        ax.set_ylim(ylim)
        ax.set_ylabel(ylabel)
        ax.set_xlabel(xlabel)
        remove_spines(ax)
        format_ticks(ax)
        
        unit_dim = 0.4
        axins1 = ax.inset_axes([0.5, 0.8, unit_dim*width_scale, unit_dim*height_scale])
        axins2 = ax.inset_axes([0.5, 0.6, unit_dim*width_scale, unit_dim*height_scale])
        for a in [axins1, axins2]:
            a.get_xaxis().set_visible(False)
            a.get_yaxis().set_visible(False)
        

        bxy = (np.max(xlim)*0.3, np.max(ylim)*0.4)
        dxy = (np.max(xlim)*0.3, np.max(ylim)*0.3)
        cxy = (np.max(xlim)*0.3, np.max(ylim)*0.2)

        ax.annotate(f'k={np.round(b, 3)}', xy=bxy)
        ax.annotate(f'Gen. from sen.={int(dist_from_sen)}', xy=dxy)
        ax.annotate(f'Chase frame={chase_frame}', xy=cxy)

        ax.scatter(x, y, **scatterkwargs)
        ax.plot(x, y_pred, **linekwargs)

        axins1.imshow(cellstacksdict['bf'][i])
        axins2.imshow(cellstacksdict['yfp'][i], **imkwargs)

        chase_hr = chase_frame*time_delta_hrs
        hr_rel_chase = i*time_delta_hrs - chase_hr
        axins1.set_title(f'{str(np.round(hr_rel_chase, 2)).zfill(4)} hrs. post chase', fontsize=7)

        plt.tight_layout()

        savepath = os.path.join(savedir, f'{basename}_frame{str(i).zfill(3)}.tif')
        savepaths.append(savepath)
        fig.savefig(savepath)
        plt.close(fig)

    allframes = [io.imread(path) for path in savepaths]
    plotstack = io.concatenate_images(allframes)
    plotstacksavepath = os.path.join(savedir, f'{basename}.tif')
    io.imsave(plotstacksavepath, plotstack)
    for path in savepaths:
        os.remove(path)
    print(f'Saved plot stack at\n{plotstacksavepath}')


def save_skimage_stack_as_mp4(filepaths, savepath, **kwargs):
    """
    Read in the image files at <filepaths> as cv2 image objects
    then encode as mp4 video and save at <savepath>

    Return nothing, print save location
    """
    # Save as .mp4
    cv2_frames = [cv2.imread(filepath) for filepath in filepaths]
    height, width, layers = cv2_frames[0].shape
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    fps = 5
    video = cv2.VideoWriter(savepath, fourcc, fps,
                            (width, height))
    for img in cv2_frames:
        video.write(img)
    print(f'mp4 saved at\n{savepath}')
    # Release file from memory
    video.release()

def save_segmentation_visualization(allframesdf, channel, cellstacksdict, draw_outline=True):
    """
    Using the segmentation data found in <allframesdf>, annotate cell outline
    ROI and write timestamp on each frame of the <channel> stack in 
    <cellstacksdict>. Save output as mp4 and tif stack

    Returns nothing
    """
    outlinedf = segmentation.save_outline_rois_df(allframesdf, write_df=False, return_outline_df=True)
    frames_table = allframesdf.set_index('frame_rel')
    cell_index = frames_table.cell_index.iloc[0]
    cellname = f'{frames_table.compartment_name.iloc[0]}_cell{str(cell_index).zfill(3)}_roi_{channel}_stack.tif'
    savepath = os.path.join(frames_table.compartment_dir.iloc[0], f'{cellname}')
    if not draw_outline:        
        savepath = savepath.replace('.tif', '_no-outline.tif')
    
    t0_frame = 0
    linewidth = 0.4
    fontsize = 8
    filenames = []
    
    for frame_index in frames_table.index.unique():
        frame_index = int(frame_index)
        frame_img = cellstacksdict[cell_index][frame_index]
        interval_hours = 10/60
        hours = (frame_index - t0_frame)*interval_hours
        fig, ax = figure_ax(
            dpi=300,
            height_scale=2*frame_img.shape[0]/300,
            width_scale=2*frame_img.shape[1]/300)
        ax.imshow(frame_img, cmap=plt.get_cmap('gray'))
        xs = list(outlinedf.loc[outlinedf.frame==frame_index, 'x'])
        ys = list(outlinedf.loc[outlinedf.frame==frame_index, 'y'])
        #repeat the first point to create a 'closed loop'
        xs.append(xs[0])
        ys.append(ys[0])

        if draw_outline:
            ax.plot(xs, ys, color='white', linewidth=linewidth)
        # Timestamp the frame
        ax.annotate(f'{np.round(hours, 2)} hrs', xy=(5, 35), color='white', fontsize=fontsize)
        # Get rid of the axis since we only want to see the image
        plt.axis('off')
        filename = f'frame{str(frame_index).zfill(3)}.tif'
        # Save the frame. It will be read in and deleted later after 
        # saving stack as tif and mp4
        fig.savefig(filename, bbox_inches=0)
        filenames.append(filename)
        plt.close(fig)
        plt.close('all')


    frames = []
    for filename in filenames:
        img = io.imread(filename)
        frames.append(img)
    filepaths = [os.path.join(os.getcwd(), fname) for fname in filenames]
    save_skimage_stack_as_mp4(filepaths, savepath.replace('.tif', '.mp4'))
    stack = io.concatenate_images(frames)
    print(f'Saved stack at \n{savepath}')
    # Save original tif stack
    io.imsave(savepath, stack)
    # Delete tifs for individual frames now that the stack has been saved
    for filename in filenames:
        os.remove(filename)

def filename_from_kwargs(kwargs_dict, ext='.png'):
    """
    Using the dict in <kwargs> (which are kwargs passed to a 
    seaborn plotting method e.g. sns.lineplot(**kwargs)),
    create a filename for the plot and return it
    """
    excluded_kws = ['data', 'estimator', 'hue_order', 'palette', 'ax', 'hue_order', 'order']
    filename = '_'.join([f'{key}={val}' for key, val in kwargs_dict.items() if key not in excluded_kws])
    filename = f'{filename}{ext}'
    return filename

def save_figure(fig, kwargs_dict, ext='.png'):
    """
    Save the figure to byc.constants.source_dir.plots_dir
    and print the saved location
    """
    filename = filename_from_kwargs(kwargs_dict, ext=ext)
    savepath = os.path.abspath(os.path.join(constants.plots_dir, filename))
    fig.savefig(savepath)
    print(f'Saved figure at\n{savepath}')

def get_dist_from_sen_color_map(max_dist_from_sen=20):
    colors = sns.color_palette('viridis_r', max_dist_from_sen+1)
    return colors

def get_pre_post_sep_palette(
        max_dist_from_sen=20,
        index_for_post_sep=2,
        index_for_pre_sep=None):
    if not index_for_pre_sep:
        index_for_pre_sep = max_dist_from_sen - 3
    colors = sns.color_palette('viridis_r', max_dist_from_sen+1)
    colors_list = [c for c in colors]
    pre_post_SEP_palette = [colors_list[index_for_pre_sep], colors_list[index_for_post_sep]]

    return pre_post_SEP_palette

def plot_dist_from_sen_palette_key(
        max_dist_from_sen=20,
        major_tick_space=5,
        minor_tick_space=1
    ):

    colors = sns.color_palette('viridis_r', max_dist_from_sen+1)
    colors_list = [c for c in colors]

    # Plot the continuous legend
    fig, ax = figure_ax(height_scale=0.2, width_scale=0.75)
    fig.set_dpi(300)
    alpha = 0.8
    x1 = np.arange(0, max_dist_from_sen, 1)
    y = np.full(len(x1), 1)
    colors_list_alpha = [matplotlib.colors.to_rgba(c, alpha) for c in colors_list]
    for i, color in enumerate(colors_list_alpha):
        x = np.arange(i-0.5, i+1.5, 1)
        y = np.full(len(x), 0.5)
        ax.fill_between(x, y, color=color, edgecolor=None)

    ax.set_ylim(0, 1)
    ax.set_xlim( max_dist_from_sen + 0.5, -0.5)
    xticks = list(np.arange(0, max_dist_from_sen + 0.5, major_tick_space))
    xticks.reverse()
    ax.set_xticks(xticks)
    ax.spines['left'].set_visible(False)
    ax.yaxis.set_visible(False)
    format_ticks(ax, tickdirection='out', xminorspace=minor_tick_space)
    ax.set_xlabel("Buds before death")

    filetype = '.svg'
    savepath = os.path.join(os.getcwd(), f'plots\\gen_from_death_color_palette_max_dist_from_sen={max_dist_from_sen}{filetype}')
    fig.set_tight_layout(True)
    fig.savefig(savepath)
    print(f'Saved figure at\n{savepath}')