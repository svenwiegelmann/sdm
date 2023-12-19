"""
Plot Preferences Module for Matplotlib

Author: Sven Wiegelmann,
        Leibniz University Hannover,
        Institute of Electric Power Systems,
        Electric Energy Storage Systems Section,
        Appelstraße 9A,
        30167 Hannover,
        Germany
Version: 12.12.2023

Overview:
This module configures and customizes aesthetic aspects of matplotlib plots, 
offering functions for setting default preferences, adjusting plot elements, 
and enhancing visual presentation for publication. It includes advanced 
functionalities like custom tick locators and formatters.

Primary Functions:
- plot_pre: Configure default plot preferences, including color palette.
- hide_ticklabels: Remove tick labels from a specified axis.
- print_num_subplot: Annotate subplots with labels for clarity.
- pub_export: Save figures in PGF and PNG formats for publications.
- axis_add_custom_ticks: Add and label custom ticks on plot axes.

Secondary Functions:
- _int_sqrt: Calculate integer square root approximations.
- _set_lim: Extend plot limits symmetrically based on a factor.

Classes:
- AdditionalTickLocator: Adds additional ticks to an existing locator.
- AdditionalTickFormatter: Applies custom formatting to specific ticks.

Usage:
Designed for creating consistent, publication-quality plots. Users should 
invoke `plot_pre` at the start of scripts for uniform plot aesthetics. 
Familiarity with matplotlib is recommended for effective utilization.
"""


import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

import math
import numpy as np


#%%
def plot_pre():
    """
    Configure default aesthetic preferences for matplotlib plots.

    This function customizes various aesthetic aspects of matplotlib plots,
    including colors, line styles, and figure size. It includes configurations
    for LaTeX text rendering and serif fonts, but these are currently commented
    out to ensure compatibility across different environments. The primary focus
    is on setting general plot parameters and defining a custom color palette 
    to ensure a consistent and visually appealing style across all plots.

    The configuration includes:
    - Setting figure and axes background colors.
    - Updating plot parameters like font sizes and legend styles.
    - Defining a custom color palette 'farbe' for use in plots.
    - Specifying a range of line styles 'lstyle' for plot lines.
    - Setting a standard figure size 'fig_size'.
    - LaTeX settings for text rendering are present but commented out.

    Returns:
    tuple: Contains the following elements:
           - farbe (list): List of color codes for the custom palette.
           - z (int): An integer for layering of plot elements, initialized as -1.
           - g (dict): An empty dictionary for potential future use.
           - fig_size (np.array): Default figure size.
           - lstyle (list): List of line styles.

    Note:
    The function should be called at the start of scripts involving plotting
    to ensure uniformity in the visual style of the plots. Uncomment the LaTeX
    settings if LaTeX compatibility is confirmed in the environment.
    """
    
    # LaTeX renderings
    plt.rc('text', usetex=True)
    # plt.rc('font', family='serif')
    # plt.rc('pgf', texsystem='pdflatex')
    # plt.rc('pgf', rcfonts=False)
    # plt.rc('mathtext', fontset='cm')
    
    plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0), # RGB, Alpha
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.0),
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),
    
    # defaults
    'font.size': 11.0,              # 10.0
    'axes.titlesize':'large',
    'axes.labelsize': 'medium',
    'figure.titlesize': 'large',
    'xtick.labelsize':'medium',
    'ytick.labelsize':'medium',
    'legend.fontsize': 'medium',
    'legend.title_fontsize': None
    })
    
    # font_scalings = {
    # 'xx-small' : 0.579,
    # 'x-small'  : 0.694,
    # 'small'    : 0.833,
    # 'medium'   : 1.0,
    # 'large'    : 1.200,
    # 'x-large'  : 1.440,
    # 'xx-large' : 1.728,
    # 'larger'   : 1.2,
    # 'smaller'  : 0.833,
    # None       : 1.0}
    
    # reset rcParams
    # mpl.rcParams.update(mpl.rcParamsDefault)
        
    farbe = [[0, 0, 0],                      # 00 - black
              [85/255, 85/255, 85/255],      # 01 - dark grey <>
              [238/255, 119/255, 51/255],    # 02 - orange <>
              [0, 68/255, 136/255],          # 03 - blue <>
              [204/255, 51/255, 17/255],     # 04 - red <>
              [34/255, 136/255, 51/255],     # 05 - green <>
              [170/255, 68/255, 153/255],    # 06 - purple <>
              [221/255, 170/255, 51/255],    # 07 - yellow <>
              [238/255, 51/255, 119/255],    # 08 - magenta <>
              [102/255, 102/255, 51/255],    # 09 - dark yellow <>
              [51/255, 187/255, 238/255],    # 10 - cyan <>
              [0, 153/255, 136/255],         # 11 - teal <>
              [51/255, 34/255, 136/255],     # 12 - indigo <>
              [153/255, 153/255, 51/255],    # 13 - olive <>
              [136/255, 34/255, 85/255],     # 14 - wine <>
              [68/255, 187/255, 153/255]]*10 # 15 - mint <>
    
    lstyle =   [(0, (1, 1)),                # 0  densely dotted
                (0, (5, 1)),                # 1  densely dashed
                (0, (3, 1, 1, 1)),          # 2  densely dashdotted
                (0, (3, 1, 1, 1, 1, 1)),    # 3  densely dashdotdotted
    
                (0, (1, 5)),                # 4  dotted
                (0, (5, 5)),                # 5  dashed
                (0, (3, 5, 1, 5)),          # 6  dashdotted
                (0, (3, 10, 1, 10, 1, 10)), # 7  dashdotdotted
                
                (0, (1, 10)),               # 8  loosely dotted
                (0, (5, 10)),               # 9  loosely dashed
                (0, (3, 10, 1, 10)),        # 10 loosely dashdotted
                (0, (3, 5, 1, 5, 1, 5))]    # 11 loosely dashdotdotted
                
    z = -1
    g = {}
    
    fig_size = np.array([4.75]*2)
    
    return farbe, z, g, fig_size, lstyle


#%%
def _int_sqrt(num):
    """
    Calculate integer square root approximations for a given number.

    This function computes a pair of integers (x, y) such that the product x * y 
    is greater than or equal to the input number 'num'. It finds the integer 
    square root of 'num' and adjusts the values of x and y to meet the 
    aforementioned condition. This function is particularly useful for grid or 
    matrix size calculations where the area must be at least as large as a 
    certain value.

    Parameters:
    - num (int): The number for which to find integer square root approximations.

    Returns:
    tuple: A tuple (x, y) where x and y are integer approximations satisfying x * y >= num.

    Examples:
    >>> _int_sqrt(10)
    (4, 3)  # 4 * 3 = 12, which is greater than 10

    Notes:
    - This function is useful in scenarios where exact square root values are not
      essential, but an approximate product that exceeds a threshold is required.
    """
    y = math.isqrt(num)
    x = y
    if y == 0:
        return x,y
    elif math.sqrt(num) % y != 0:
        y += 1
        if y*x < num:
            x += 1
    return x,y


#%%
def _set_lim(xmin, xmax, dx):
    """
    Calculate extended plot limits based on a given range and extension factor.

    This internal function adjusts the minimum (xmin) and maximum (xmax) values of a 
    plot range, extending this range by a specified factor (dx). The extension is 
    symmetrically applied to both sides of the range, effectively 'padding' the plot 
    limits. This is particularly useful for ensuring plot elements near the edges 
    are clearly visible and not clipped by the plot boundaries.

    Parameters:
    - xmin (float): The minimum value of the original plot range.
    - xmax (float): The maximum value of the original plot range.
    - dx (float): The factor by which to extend the plot range, applied symmetrically.

    Returns:
    list: A list containing the new minimum and maximum values of the extended plot range.

    Example:
    >>> _set_lim(0, 10, 0.1)
    [-1.0, 11.0]  # Extends the range from 0-10 to -1 to 11

    Note:
    - Designed for internal module use, this function is beneficial for plot 
      adjustments where extra space around the primary data range enhances clarity.
    """
    return [xmin,xmax] + np.multiply(np.diff([xmin,xmax])*dx,[-1,1])


#%%
def hide_ticklabels(ax):
    """
    Hide all tick labels on a specified matplotlib axis.

    This function removes the visibility of both major tick marks and their 
    corresponding labels on the given axis. It is useful in creating plots where 
    the tick labels are not necessary for interpretation or where a minimalistic 
    visual style is desired.

    Parameters:
    - ax (matplotlib.axes.Axes): The Axes object on which tick labels will be hidden.

    Returns:
    None: Modifies the Axes object in place, hiding its tick labels.

    Example:
    >>> fig, ax = plt.subplots()
    >>> hide_ticklabels(ax)
    # The ax object now has hidden tick labels

    Note:
    - This function is particularly useful in multi-plot layouts where axis 
      labels might be redundant or in cases where the data presentation 
      benefits from a simplified axis appearance.
    """
    for tick in ax.get_major_ticks():
        tick.tick1line.set_visible(False)
        tick.tick2line.set_visible(False)
        tick.label1.set_visible(False)
        tick.label2.set_visible(False)


#%%
def print_num_subplot(ax, an_letter, 
                      bbox_args=dict(boxstyle="round", color='1', fc='1', ec='None', lw=0)):
    """
    Annotate a subplot with a label (e.g., 'a', 'b', 'c') in a matplotlib Axes object.

    This function adds a text annotation to a subplot, usually to label it as part of
    a larger figure composed of multiple subplots. The label is typically a letter 
    like 'a', 'b', 'c', etc. The annotation is placed at the top right corner of the 
    subplot.

    Parameters:
    - ax (matplotlib.axes.Axes): The Axes object to annotate.
    - an_letter (str): The letter or text used for the subplot label.
    - bbox_args (dict, optional): Dictionary of box style properties for the annotation. 
                                  Defaults to a rounded, white background box with no edge.

    Returns:
    None: The function adds an annotation to the Axes object but does not return any value.

    Example:
    >>> fig, ax = plt.subplots()
    >>> print_num_subplot(ax, 'a')
    # The ax object now has a label 'a' at the top right corner

    Note:
    - This function is particularly useful in figures with multiple subplots to 
      provide clear labels for each subplot, enhancing the readability of 
      complex figures.
    """
    ax.annotate(r'\textbf{{({})}}'.format(an_letter),
                (1,1),
                xycoords='axes fraction',
                xytext=(-8,-8),
                textcoords="offset points",
                size='x-large',
                c='k',
                ha='right',
                va='top',
                bbox=bbox_args,
                zorder=100)


#%%
def pub_export(fig, name, path='./img/', vers=''):
    """
    Save a matplotlib figure in both PGF and PNG formats for publication.

    This function exports a given matplotlib figure as a .pgf file and a .png file,
    suitable for use in publications. It allows for the inclusion of an optional 
    version string in the file name, enabling easy version control of figure files.
    The function automatically sets the backend to 'pgf' if it is not already set.

    Parameters:
    - fig (matplotlib.figure.Figure): The figure to be saved.
    - name (str): Base name for the saved figure files.
    - path (str, optional): Path prefix to be put in front of the file name. Defaults to './img/'.
    - vers (str, optional): Version suffix to be appended to the file name. Defaults to ''.

    Returns:
    None: The function saves the figure to the specified files but does not return any value.

    Example:
    >>> fig, ax = plt.subplots()
    >>> pub_export(fig, 'example_plot', path='./img/', vers='v1')
    # Saves './img/example_plot_v1.pgf' and './img/example_plot_v1.png'

    Note:
    - Ensure the figure is fully prepared and finalized before calling this function.
    - This function is ideal for preparing figures for publication, offering both
      vector (PGF) and raster (PNG) formats.
    """
    if mpl.get_backend() != 'pgf':
        mpl.use('pgf')
    
    fig.savefig(path + name + (('_'+vers) if vers else '') + '.pgf')
    fig.savefig(path + name + (('_'+vers) if vers else '') + '.png')


#%% Add custom Ticks
# https://stackoverflow.com/questions/22245949/adding-a-custom-tick-and-label

class AdditionalTickLocator(mticker.Locator):
    """
    A custom locator for adding additional ticks to an existing locator in Matplotlib.

    This locator works by chaining an existing locator and then adding custom ticks
    to the results produced by that locator.

    Attributes:
    - chain (mticker.Locator): The original locator to chain.
    - additional_ticks (np.ndarray): Array of additional tick positions.

    Methods:
    - tick_values(vmin, vmax): Returns the tick values including additional ticks.
    - nonsingular(v0, v1): Delegates to the chained locator's nonsingular method.
    - set_params(**kwargs): Sets parameters for the locator.
    - view_limits(vmin, vmax): Delegates to the chained locator's view_limits method.
    """
    def __init__(self, chain: mticker.Locator, ticks) -> None:
        super().__init__()
        assert chain is not None
        self._chain = chain
        self._additional_ticks = np.asarray(list(ticks))

    def _add_locs(self, locs):
        locs = np.unique(np.concatenate([
            np.asarray(locs),
            self._additional_ticks
        ]))
        return locs

    def tick_values(self, vmin, vmax):
        locs = self._chain.tick_values(vmin, vmax)
        return self._add_locs(locs)

    def __call__(self):
        # this will call into chain's own tick_values,
        # so we also add ours here
        locs = self._chain.__call__()
        return self._add_locs(locs)

    def nonsingular(self, v0, v1):
        return self._chain.nonsingular(v0, v1)
    def set_params(self, **kwargs):
        return self._chain.set_params(**kwargs)
    def view_limits(self, vmin, vmax):
        return self._chain.view_limits(vmin, vmax)


class AdditionalTickFormatter(mticker.Formatter):
    """
    A custom formatter for applying special formatting to custom ticks in Matplotlib.

    This formatter chains an existing formatter and applies custom formatting to specific
    tick values.

    Attributes:
    - chain (mticker.Formatter): The original formatter to chain.
    - additional_ticks (dict): Dictionary mapping tick positions to their custom labels.

    Methods:
    - __call__(x, pos=None): Returns the formatted label for tick value x.
    - format_data_short(value): Formats a value for short representation.
    - get_offset(): Returns the offset for the formatter.
    - set_locs(locs): Sets the locations for the formatter.
    """
    def __init__(self, chain: mticker.Formatter, ticks) -> None:
        super().__init__()
        assert chain is not None
        self._chain = chain
        self._additional_ticks = ticks

    def __call__(self, x, pos=None):
        if x in self._additional_ticks:
            return self._additional_ticks[x]
        res = self._chain.__call__(x, pos)
        return res

    def format_data_short(self, value):
        if value in self._additional_ticks:
            return self.__call__(value)
        return self._chain.format_data_short(value)

    def get_offset(self):
        return self._chain.get_offset()
    
    def _set_locator(self, locator):
        self._chain._set_locator(locator)

    def set_locs(self, locs):
        self._chain.set_locs(locs)


def axis_add_custom_ticks(axis, ticks):
    """
    Add custom ticks and labels to a specified axis in a Matplotlib plot.
    
    This function uses the AdditionalTickLocator and AdditionalTickFormatter classes
    to add custom tick marks and their labels to a given axis.
    
    Parameters:
    - axis (matplotlib.axis.Axis): The axis to which custom ticks will be added.
    - ticks (dict): A dictionary mapping custom tick positions (keys) to their labels (values).
    
    Example:
    >>> fig, ax = plt.subplots()
    >>> custom_ticks = {1.5: '1.5', 2.5: '2.5'}
    >>> axis_add_custom_ticks(ax.xaxis, custom_ticks)
    # Adds custom ticks at positions 1.5 and 2.5 with respective labels.
    
    Note:
    - This function is useful for plots where standard tick marks do not suffice or
      where specific data points need to be highlighted with custom annotations.
    """
    locator = axis.get_major_locator()
    formatter = axis.get_major_formatter()
    axis.set_major_locator(AdditionalTickLocator(locator, ticks.keys()))
    axis.set_major_formatter(AdditionalTickFormatter(formatter, ticks))