"""
System Design Method Functions (SDM_fcn.py)

Author:     Sven Wiegelmann,
            Leibniz University Hannover,
            Institute of Electric Power Systems,
            Electric Energy Storage Systems Section,
            Appelstraße 9A,
            30167 Hannover,
            Germany
Version:    12.12.2023

Overview:
This module offers a comprehensive suite of functions for conceptual design and 
analysis of cell-based energy systems. The System Design Method (SDM) is built upon 
the combined Extended Ragone Plot (ERP) and Constraint Satisfaction Problem (CSP) 
frameworks. It facilitates efficient system design and optimization, focusing on 
energy storage and management. The framework is implemented using a matplotlib-based 
slider widget, allowing users to interactively manipulate inputs and visualize the 
impact on the system design.

Primary Functions:
- calc_P_UI
- calc_design
- calc_intersection_curve
- compile_EP
- findIntersection
- get_paths
- hyperbel
- interp_dp
- inv_linear
- linear (overloaded function)
- plot_EP_limit
- plot_constraint_intersection_trajectory
- plot_dss
- plot_iso_intersection_curves
- plot_limit_loci
- recalc_design
- recalc_limits
- replot_EP_limits
- update_eta

Secondary Functions:
- _get_fcn_angle_loglog
- _get_largest_smaller
- _get_smallest_larger
- _plot_iso_CPapp
- _plot_isochrones
- _print_design
- _print_limit_values
- _set_intersection_limits
- _text_EP_units

Required Packages and Modules:
- matplotlib (for plotting and slider widget implementation)
- numpy (for numerical calculations)
- [Other necessary packages/modules]

Usage:
Designed for energy engineering researchers and practitioners, this module aids in 
the design and evaluation of cell-based energy systems, offering interactive and 
advanced analysis tools.

Note:
Users should possess a basic understanding of energy system concepts and Python 
programming for effective utilization of this tool. Familiarity with matplotlib 
and its widget functionalities is advantageous. Refer to individual function 
docstrings for detailed information on their usage and parameters.
"""


import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
import matplotlib.collections as clt
from matplotlib.patches import Rectangle, Polygon
from tol_colors import tol_cmap, tol_cset

from scipy import interpolate, optimize
import numpy as np
import copy

import plot_preferences


#%%
class UpdatablePatchCollection(clt.PatchCollection):
    """
    A subclass of matplotlib.collections.PatchCollection that allows for dynamic 
    updating of patch paths in a collection.

    This class is designed to facilitate the modification of patches in a plot, 
    enabling the paths of these patches to be updated as needed. This feature 
    is particularly useful in scenarios where the graphical elements need to be 
    dynamically changed or updated based on new data or user interactions.

    Methods:
    - __init__: Initialize the UpdatablePatchCollection with a set of patches.
    - get_paths: Updates and returns the current paths of the patches.

    Attributes:
    - patches (list): A list of patch objects to be included in the collection.

    Example Usage:
    ```python
    patches = [Rectangle(...), Circle(...), ...]
    updatable_collection = UpdatablePatchCollection(patches, ...)
    ax.add_collection(updatable_collection)
    # Later updates
    patches[0].set_width(new_width)
    updatable_collection.get_paths()  # Update the collection with new paths
    ```

    Inherits From:
    - matplotlib.collections.PatchCollection
    """
    def __init__(self, patches, *args, **kwargs):
        """
        Initialize the UpdatablePatchCollection with a set of patches.

        Parameters:
        - patches (list): A list of matplotlib.patches objects.
        - *args: Variable length argument list.
        - **kwargs: Arbitrary keyword arguments.
        """
        self.patches = patches
        clt.PatchCollection.__init__(self, patches, *args, **kwargs)

    def get_paths(self):
        """
        Update the paths of the patches in the collection to reflect any changes.

        Returns:
        - list: The updated list of paths derived from the current patches.
        """
        self.set_paths(self.patches)
        return self._paths


#%%
def compile_EP(tmp_res, calc_EPsys=False, prefix_str='cell_model[1,1].'):
    """
    Compiles continuous arrays of data, including energy and power, from 
    simulation results into a structured format for Ragone plot analysis. 
    This function processes results to categorize and aggregate data based on 
    conditions like voltage, Crate, and temperature limits. It supports both 
    system-level (E_sys) and cell-level (E_cell) analysis, compiling data for 
    different states of energy storage, crucial for forming a Ragone plot.

    Parameters:
    - tmp_res (dict): Dictionary with simulation results. Each key corresponds 
      to a unique result set with detailed sub-keys.
    - calc_EPsys (bool, optional): If True, calculates system-level data 
      (E_sys). Defaults to False for cell-level data (E_cell).
    - prefix_str (str, optional): String identifier for the specific sub-model
      in result keys. Defaults to 'cell_model[1,1].'.

    Returns:
    - tuple: Containing three elements:
        1. Dictionary with keys representing different states of energy storage, 
           each containing continuous arrays of data like energy and power 
           under various conditions.
        2. List of tuples with indices and corresponding states of energy from 
           input results.
        3. Minimum result number (rn_min) where 'time' is a key in results.

    Note:
    - Assumes specific key patterns and structures in `tmp_res`.
    - Uses NumPy for numerical operations, handling NaN values for error 
      conditions.
    - Categorizes data based on specified conditions, essential for Ragone 
      plot analysis.
    """
    
    d = {}
    
    tmp_ind = [(rn,tmp_res[rn]['_nP'],tmp_res[rn]['_nSoC']) for rn in range(len(tmp_res))]
    rn_min = np.nanmin([rn if 'time' in tmp_res[rn].keys() else np.nan for rn in range(len(tmp_res))]).astype(int)
    
    for nSoC in range(max(tmp_ind)[2]+1):
        e0 = []
        e1 = []
        e2 = []
        e3 = []
        e5 = []
        p = []
        soc = []
        for nP in range(max(tmp_ind)[1]+1):
            # if nP == 20:
            #     continue
            # print(nP,nSoC)
            rn = [item[0] for item in tmp_ind if item[1]==nP and item[2]==nSoC][0]
            # print('=== {} ==='.format(rn))
            try:
                if calc_EPsys:
                    tmp_E = abs(tmp_res[rn]['E_sys'][-1])
                    tmp_P = abs(tmp_res[rn]['P_sys'][0])
                else:
                    tmp_E = abs(tmp_res[rn][prefix_str+'E_cell'][-1])
                    tmp_P = abs(tmp_res[rn][prefix_str+'P_cell'][0])
                
                # print(e1, not e1)
                # print(res[rn][prefix_str+'U_cell'][-1])
                cond1 = tmp_res[rn][prefix_str+'U_cell']-tmp_res[rn]['U_max'][0]>0
                cond2 = tmp_res[rn][prefix_str+'U_cell']-tmp_res[rn]['U_min'][0]<0
                cond3 = tmp_res[rn][prefix_str+'Crate']-tmp_res[rn]['Crate_max'][0]>0
                # cond4 = tmp_res[rn][prefix_str+'Crate']-tmp_res[rn]['Crate_min'][0]>0
                cond5 = tmp_res[rn][prefix_str+'T_cell']-tmp_res[rn]['T_max'][0]>0
                # cond6 = tmp_res[rn][prefix_str+'T_cell']-tmp_res[rn]['T_min'][0]<0
                
                
                if not e0 or not any(ele == 0 for ele in e0):
                    e0.append(tmp_E)
                    if cond1[0]: #res[rn][prefix_str+'U_cell'][0]-res[rn]['U_max'][0]<=0:    # if U(t=0) < U_max
                        e1.append(tmp_E)
                        e2.append(np.nan)
                        e3.append(np.nan)
                        e5.append(np.nan)
                        
                    elif cond2.any() and cond3.any() and cond5.any():
                        e1.append(np.nan)
                        if np.where(cond2)[0].min() < np.where(cond3)[0].min() and np.where(cond2)[0].min() < np.where(cond5)[0].min():
                            e2.append(tmp_E)
                            e3.append(np.nan)
                            e5.append(np.nan)
                        elif np.where(cond3)[0].min() < np.where(cond2)[0].min() and np.where(cond3)[0].min() < np.where(cond5)[0].min():
                            e2.append(np.nan)
                            e3.append(tmp_E)
                            e5.append(np.nan)
                        else:
                            e2.append(np.nan)
                            e3.append(np.nan)
                            e5.append(tmp_E)
                    
                    elif cond2.any() and cond3.any():
                        e1.append(np.nan)
                        e5.append(np.nan)
                        if np.where(cond2)[0].min() < np.where(cond3)[0].min():
                            e2.append(tmp_E)
                            e3.append(np.nan)
                        else:
                            e2.append(np.nan)
                            e3.append(tmp_E)                        
                    elif cond2.any() and cond5.any():
                        e1.append(np.nan)
                        e3.append(np.nan)
                        if np.where(cond2)[0].min() < np.where(cond5)[0].min():
                            e2.append(tmp_E)
                            e5.append(np.nan)
                        else:
                            e2.append(np.nan)
                            e5.append(tmp_E)
                    elif cond3.any() and cond5.any():
                        e1.append(np.nan)
                        e2.append(np.nan)
                        if np.where(cond3)[0].min() < np.where(cond5)[0].min():
                            e3.append(tmp_E)
                            e5.append(np.nan)
                        else:
                            e3.append(np.nan)
                            e5.append(tmp_E)
            
                    elif cond2.any():
                        e1.append(np.nan)
                        e2.append(tmp_E)
                        e3.append(np.nan)
                        e5.append(np.nan)
                    elif cond3.any():
                        e1.append(np.nan)
                        e2.append(np.nan)
                        e3.append(tmp_E)
                        e5.append(np.nan)
                    elif cond5.any():
                        e1.append(np.nan)
                        e2.append(np.nan)
                        e3.append(np.nan)
                        e5.append(tmp_E)
                   
                    else:
                        e1.append(np.nan)
                        e2.append(np.nan)
                        e3.append(np.nan)
                        e5.append(np.nan)
                else:
                    e0.append(np.nan)
                    e1.append(np.nan)
                    e2.append(np.nan)
                    e3.append(np.nan)
                    e5.append(np.nan)
            except:
                e0.append(np.nan)
                e1.append(np.nan)
                e2.append(np.nan)
                e3.append(np.nan)
                e5.append(np.nan)
                
            if nSoC == 0:
                try:
                    p.append(tmp_P)
                except:
                    p.append(np.nan)
                try:
                    soc.append(tmp_res[rn][prefix_str+'SoC'][0])
                except:
                    soc.append(np.nan)
            d[nSoC] = [e0,e1,e2,e3,e5] # '.': U(t=0)<Umax, 'o': U(t=0)<Umax & U(t=end)>Umin, 'v': U(t=0)<Umax & Crate(t=end)<Cratemax [beendet durch Cratemax], 's': U(t=0)<Umax & T(t=end)<Tmax [beendet durch Tmax]
            if nSoC == 0:
                d['p'] = p
                d['soc'] = soc
    return d, tmp_ind, rn_min



#%%
def recalc_limits(res, Umax, Umin, Cratemax, Tmax, SoCmin=-0.1, calc_EPsys=False, 
                  prefix_str='cell_model[1,1].'):
    """
    Adjusts the limits of various parameters in simulation results and 
    recalculates data arrays based on these new limits. Handles both 
    system-level and cell-level data.

    Parameters:
    - res (list): List of dictionaries, each a set of simulation results.
    - Umax (float): Maximum voltage limit.
    - Umin (float): Minimum voltage limit.
    - Cratemax (float): Maximum Crate limit.
    - Tmax (float): Maximum temperature limit.
    - SoCmin (float, optional): Minimum state of charge limit. Defaults
      to -0.1.
    - calc_EPsys (bool, optional): If True, recalculates system-level 
      data. False for cell-level data.
    - prefix_str (str, optional): Substring for specific sub-model in 
      result keys. Defaults to 'cell_model[1,1].'.

    Returns:
    - list: Updated simulation results with new limits.

    Note:
    - Deep copies `res` to avoid mutating original data.
    - Uses NumPy for numerical operations and handling NaN values.
    - Involves conditional checks for different scenarios.
    """

    loc_res = copy.deepcopy(res)
    for rn,_ in enumerate(loc_res):
        # if rn == 8:
        #     0
        try:
            if calc_EPsys:
                tmp_E = abs(loc_res[rn]['E_sys'][-1])
                tmp_P = abs(loc_res[rn]['P_sys'][0])
            else:
                tmp_E = abs(loc_res[rn][prefix_str+'E_cell'][-1])
                tmp_P = abs(loc_res[rn][prefix_str+'P_cell'][0])
            
            loc_res[rn]['U_min'] = np.array([Umin]*2)
            loc_res[rn]['U_max'] = np.array([Umax]*2)
            loc_res[rn]['Crate_max'] = np.array([Cratemax]*2)
            loc_res[rn]['T_max'] = np.array([Tmax]*2)
            loc_res[rn]['SoC_min'] = np.array([SoCmin]*2)

            
            cond1 = loc_res[rn][prefix_str+'U_cell']-loc_res[rn]['U_max'][0]>0
            cond2 = loc_res[rn][prefix_str+'U_cell']-loc_res[rn]['U_min'][0]<0
            cond3 = loc_res[rn][prefix_str+'Crate']-loc_res[rn]['Crate_max'][0]>0    
            # cond4 = loc_res[rn][prefix_str+'Crate']-loc_res[rn]['Crate_min'][0]>0
            cond5 = loc_res[rn][prefix_str+'T_cell']-loc_res[rn]['T_max'][0]>0
            # cond6 = loc_res[rn][prefix_str+'T_cell']-loc_res[rn]['T_min'][0]<0
            # cond7 = loc_res[rn][prefix_str+'SoC']-loc_res[rn]['SoC_max'][0]<0
            cond8 = loc_res[rn][prefix_str+'SoC']-loc_res[rn]['SoC_min'][0]<0
            
            # ind_s ("-1" overlapping due to symbols)
            if cond1.any():
                ind_s = np.where(cond1)[0].max()-1 if np.where(cond1)[0].max() > 0 else 0 
            else:
                ind_s = 0
 
            # ind_e ("+1" overlapping due to symbols)
            if cond2.any() and cond3.any() and cond5.any() and cond8.any():
                ind_e = min(np.where(cond2)[0].min(),
                            np.where(cond3)[0].min(),
                            np.where(cond5)[0].min(),
                            np.where(cond8)[0].min())+1

            elif cond2.any() and cond3.any() and cond5.any():
                ind_e = min(np.where(cond2)[0].min(),
                            np.where(cond3)[0].min(),
                            np.where(cond5)[0].min())+1            
            elif cond2.any() and cond3.any() and cond8.any():
                ind_e = min(np.where(cond2)[0].min(),
                            np.where(cond3)[0].min(),
                            np.where(cond8)[0].min())+1
            elif cond2.any() and cond5.any() and cond8.any():
                ind_e = min(np.where(cond2)[0].min(),
                            np.where(cond5)[0].min(),
                            np.where(cond8)[0].min())+1
            elif cond3.any() and cond5.any() and cond8.any():
                ind_e = min(np.where(cond3)[0].min(),
                            np.where(cond5)[0].min(),
                            np.where(cond8)[0].min())+1

            elif cond2.any() and cond3.any():
                ind_e = min(np.where(cond2)[0].min(),
                            np.where(cond3)[0].min())+1
            elif cond2.any() and cond5.any():
                ind_e = min(np.where(cond2)[0].min(),
                            np.where(cond5)[0].min())+1
            elif cond2.any() and cond8.any():
                ind_e = min(np.where(cond2)[0].min(),
                            np.where(cond8)[0].min())+1
            elif cond3.any() and cond5.any():
                ind_e = min(np.where(cond3)[0].min(),
                            np.where(cond5)[0].min())+1
            elif cond3.any() and cond8.any():
                ind_e = min(np.where(cond3)[0].min(),
                            np.where(cond8)[0].min())+1
            elif cond5.any() and cond8.any():
                ind_e = min(np.where(cond5)[0].min(),
                            np.where(cond8)[0].min())+1
            
            elif cond2.any():
                ind_e = np.where(cond2)[0].min()+1
            elif cond3.any():
                ind_e = np.where(cond3)[0].min()+1
            elif cond5.any():
                ind_e = np.where(cond5)[0].min()+1
            elif cond8.any():
                ind_e = np.where(cond8)[0].min()+1            
            else:
                ind_e = len(loc_res[rn]['time'])
            
            for key in loc_res[rn]:
                if isinstance(loc_res[rn][key],np.ndarray) and key not in loc_res[rn]['_settings']['ParameterValues']:
                    if key == prefix_str+'E_cell' or key == 'E_cell[1,1]' or key == 'E_sys': # integral quantity
                        # print(key)
                        loc_res[rn][key] = loc_res[rn][key][ind_s:ind_e] - loc_res[rn][key][ind_s]
                    else:
                        loc_res[rn][key] = loc_res[rn][key][ind_s:ind_e]
        except:
            0
    return loc_res


#%%
def calc_design(x_data, y_data, m, b=0, show_plt=True, show_annotation=True, 
                swapxy=False):
    """
    Calculates and optionally plots design points for a battery system based 
    on given data and linear parameters.

    Interpolates EP curves of a battery system and calculates intersection 
    with a linear EP line of an application. Optionally plots these curves 
    and their intersection points.

    Parameters:
    - x_data (list/array): X-coordinates of data points.
    - y_data (list/array): Y-coordinates of data points.
    - m (float): Slope of the linear function for application's EP line.
    - b (float, optional): Y-intercept of linear function. Defaults to 0.
    - show_plt (bool, optional): If True, plots EP curves and lines. 
      Defaults to True.
    - show_annotation (bool, optional): If True and plotting is enabled, 
      shows annotations on the plot. Defaults to True.
    - swapxy (bool, optional): If True, swaps x and y axes in the plot. 
      Useful for alternate representations. Defaults to False.

    Returns:
    - tuple: Intersection x-coordinates (x_int), y-coordinates (y_int), 
      and the linear function (func1).

    Note:
    - Uses `scipy.interpolate.interp1d` for interpolation and 
      `scipy.optimize.fsolve` for finding intersection points.
    - Matplotlib used for plotting; behavior changes with `show_plt` 
      and `show_annotation` flags.
    """
    
    # Calculate EP-Curves of Battery System
    funcEP = interpolate.interp1d(x_data,y_data,fill_value='extrapolate')

    # Calculate EP-Line of Application
    def linear(m,b):
        return lambda x: m*x+b
    
    func1 = linear(m,b)
    x_ref = [0,max(x_data)]

    # Calculate Design Point(s)
    h = lambda x: func1(x) - funcEP(x)
    
    x_int = optimize.fsolve(h,np.average(np.array(x_data)[~np.isnan(y_data)]))
    y_int = funcEP(x_int)
    
    tmp_obj = []
    if show_plt:
        # tmp_obj.append(plt.plot(x_data,y_data,color='r',marker='.',markersize=4,zorder=2))
        tmp_obj.append(plt.plot(x_ref, [func1(i) for i in x_ref],color='r'))
        
        for n,_ in enumerate(x_int):
            # tmp_obj.append(plt.plot([0,x_int[n],x_int[n]],[y_int[n],y_int[n],0],color='r',lw=1,ls='dotted',zorder=2))
            # tmp_obj.append(plt.plot(x_int[n],y_int[n],marker='o',color='r',markersize=4.5,zorder=2))
            if show_annotation:
                tmp_x = y_int[n] if swapxy else x_int[n] 
                tmp_y = x_int[n] if swapxy else y_int[n]
                label = '({},{})'.format(round(tmp_x,1),round(tmp_y,1))
                bbox_args = dict(boxstyle="round",color='r',fc='1',lw=1)
                tmp_obj.append(plt.annotate(label, # this is the text
                              (tmp_x,tmp_y), # these are the coordinates to position the label
                              textcoords="offset points", # how to position the text
                              xytext=(5,5), # distance from text to points (x,y)
                              ha='left',
                              va='bottom',
                              bbox=bbox_args))
        
        if swapxy:
            for i in tmp_obj:
                tmp_x = i[0].get_xdata()
                tmp_y = i[0].get_ydata()
                i[0].set_xdata(tmp_y)
                i[0].set_ydata(tmp_x)
               
    return x_int, y_int, func1



#%%
def recalc_design(res, Umax, Umin, Cratemax, Tmax, m, calc_EPsys=False,
                  show_plt=True, show_annotation=True):
    """
    Recalculates design points for a battery system based on updated parameter 
    limits and a specified linear function.

    Uses `recalc_limits` to adjust parameter limits in simulation results, 
    `compile_EP` to compile energy data, and `calc_design` to calculate and 
    optionally plot design points.

    Parameters:
    - res (list): List of dictionaries, each a set of simulation results.
    - Umax (float): Maximum voltage limit.
    - Umin (float): Minimum voltage limit.
    - Cratemax (float): Maximum Crate limit.
    - Tmax (float): Maximum temperature limit.
    - m (float): Slope of linear function for application's EP line.
    - calc_EPsys (bool, optional): If True, uses system-level data in 
      calculations. False for cell-level data.
    - show_plt (bool, optional): If True, plots design points. Defaults 
      to True.
    - show_annotation (bool, optional): If True and plotting enabled, 
      shows annotations. Defaults to True.

    Returns:
    - None: Does not return a value but optionally plots design points.

    Note:
    - Integrates `recalc_limits`, `compile_EP`, and `calc_design` for a 
      streamlined process in recalculating and visualizing battery system 
      design points.
    - Uses Matplotlib for plotting if enabled.
    """
    res_recalc = recalc_limits(res, Umax, Umin, Cratemax, Tmax)
    d_recalc, tmp_ind_recalc, rn_min_recalc = compile_EP(res_recalc,calc_EPsys=calc_EPsys)
    
    calc_design(d_recalc['p'],d_recalc[0][0],m=m,show_plt=show_plt,show_annotation=show_annotation)
    


#%%
def _text_EP_units(t):
    """
    Converts time in minutes to a formatted string in hours, minutes, seconds, 
    or milliseconds.

    This function takes time in minutes and converts it into a more readable 
    format. It formats time into hours and minutes, minutes and seconds, or 
    seconds and milliseconds, depending on the input value.

    Parameters:
    - t (float): Time in minutes.

    Returns:
    - str: A string representing formatted time. Format varies based on 't':
        - If 't' >= 60, formatted as hours and minutes (e.g., '2 h', '2:30 h').
        - If 1 <= 't' < 60, as minutes and seconds (e.g., '30 min', '30:45 min').
        - If 't' < 1, as seconds and milliseconds (e.g., '45 s', '45:300 s').

    Note:
    - Rounds the smallest unit to the nearest whole number. E.g., if 't' is 
      1.5 minutes, it's formatted as '1:30 min'.
    - Uses raw string notation ('r' before the string) for proper formatting, 
      especially with backslashes or special characters.
    """
    hours, mins = divmod(t, 60)
    if hours > 0:
        text = (r'{} h'.format(int(hours)) if mins == 0 else
                r'{}:{:02d} h'.format(int(hours), int(round(mins, 0))))
    else:
        mins, secs = divmod(mins * 60, 60)
        if mins > 0:
            text = (r'{} min'.format(int(mins)) if secs == 0 else
                    r'{}:{:02d} min'.format(int(mins), int(round(secs, 0))))
        else:
            secs, msecs = divmod(secs * 1000, 1000)
            if secs > 0:
                text = (r'{} s'.format(int(secs)) if msecs == 0 else
                        r'{}:{:03d} s'.format(int(mins), int(round(msecs, 0))))
            else:
                text = (r'{} ms'.format(int(round(msecs, 0))))
    return text



#%%
def _plot_isochrones(ax, iso_mins=[10/60, 1, 2.5, 5, 7.5, 10, 15, 30, 60, 300], 
                   pos=0.95, num_subplot=False):
    """
    Plots isochrones (equal discharge duration lines) on a matplotlib axis.

    Draws dotted lines representing discharge duration of a storage system at
    various rates on the given axis. Also annotates these lines with their 
    corresponding time durations, formatted via `_text_EP_units`.

    Parameters:
    - ax (matplotlib.axes.Axes): Axis for plotting the isochrones.
    - iso_mins (list, optional): Times in minutes for isochrones. Defaults to 
      [10/60, 1, 2.5, 5, 7.5, 10, 15, 30, 60, 300].
    - pos (float, optional): Factor to determine annotation position on each 
      isochrone. Defaults to 0.95.
    - num_subplot (bool, optional): If True, adjusts annotation position for 
      better visibility in subplots. Defaults to False.

    Returns:
    - None: Modifies the provided matplotlib axis without returning a value.

    Note:
    - Calculates angle of isochrones based on axis limits for annotation 
      positioning.
    - Rotates annotations to align with isochrone angles.
    - Employs NumPy for numerical calculations and matplotlib for plotting.
    """
    # Function implementation remains unchanged

    xlim = np.array(ax.get_xlim())
    ylim = np.array(ax.get_ylim())
    for t in iso_mins:
        xvec = xlim
        yvec = xlim*(t/60)
        ax.plot(xvec,yvec,color=[0.5]*3,ls='dotted',lw=1,alpha=0.4,zorder=0)
        
        # angle
        dx = xvec[1] - xvec[0]
        dy = yvec[1] - yvec[0]
        angle = np.rad2deg(np.arctan2(dy, dx))
        
        lim_dx = xlim[1] - xlim[0]
        lim_dy = ylim[1] - ylim[0]
        lim_angle = np.rad2deg(np.arctan2(lim_dy, lim_dx))
        
        # position
        x = xvec[1]*pos if yvec[1] <= ylim[1] else np.interp(ylim[1], yvec, xvec)*pos
        y = yvec[1]*pos if yvec[1] <= ylim[1] else ylim[1]*pos
        
        if num_subplot:
            if angle <= lim_angle*1.1 and angle >= lim_angle*0.9:
                x = x*0.935
                y = y*0.935
                
        # annotation
        text = _text_EP_units(t)
            
        ax.annotate(text,
                    (x,y),
                    ha='right',
                    va='center',
                    color=[0.5]*3,
                    alpha=1/3,
                    size='smaller',
                    transform_rotates_text=True,
                    rotation=angle,
                    rotation_mode='anchor',
                    backgroundcolor=[1]*3,
                    zorder=0)
    
        
        
#%%
def _plot_iso_CPapp(ax, c, iso_Pcell=np.arange(0, 601, 50), del_nP=[], pos=0.975, 
                   num_subplot=False):
    """
    Plots iso-lines for power-based inequality constraint (CPapp) in Graphical 
    Synthesis on a matplotlib axis.

    Draws dotted lines for cell power variations (iso_Pcell) and annotates them,
    illustrating the Constraint Satisfaction Problem for CPapp using a hyperbolic
    function.

    Parameters:
    - ax (matplotlib.axes.Axes): Axis for CPapp iso-lines.
    - c (dict): Constants/parameters for hyperbolic function.
    - iso_Pcell (array): Cell power values for iso-lines. Defaults 0-600 (50).
    - del_nP (list, optional): Indices to exclude in power values. Default empty.
    - pos (float, optional): Position factor for annotations. Default 0.975.
    - num_subplot (bool, optional): Adjusts annotations for subplots. Default False.

    Returns:
    - None: Modifies axis without returning value.

    Note:
    - For Graphical Synthesis in Constraint Satisfaction.
    - Uses 'hyperbel' for iso-line calculations.
    - Employs '_get_fcn_angle_loglog' for angle calculations in annotations.
    - NumPy for calculations, matplotlib for plotting, 'plot_preferences' for
      annotation size.
    """
    xlim = np.array(ax.get_xlim())
    ylim = np.array(ax.get_ylim())
            
    for nP,valP in enumerate(iso_Pcell):
        if nP == 0:
            continue

        f_hyp = hyperbel((c['req']['Pac_req']/c['pec']['eta_op']/valP))
        ax.plot(xlim,f_hyp(xlim),c=[0.5]*3,ls='dotted',lw=0.75,alpha=0.4,zorder=1 if 500<valP<600 else 0)
        
        xvec = xlim
        yvec = f_hyp(xvec)
        
        # position            
        if yvec[1] >= np.exp(np.log(ylim[0])+np.diff(np.log(ylim))*(1-pos)):
            x = np.exp(np.log(xlim[0])+np.diff(np.log(xlim))*(pos))
            y = f_hyp(x)
        else:
            y = np.exp(np.log(ylim[0])+np.diff(np.log(ylim))*(1-pos))
            x = f_hyp(y)

        # angle
        angle = _get_fcn_angle_loglog(ax,f_hyp)
        
        # text
        text = r'{}{:.0f} W'.format('$P_\mathrm{{cell,max}}^\mathrm{{dis}}$ = ' if nP==1 else '',valP)
        
        if (valP<=400 or divmod(nP,2)[1]==0) and nP not in del_nP:
            ax.annotate(text,
                        (x,y),
                        textcoords="offset points",
                        xytext=(0,0),
                        ha='right',
                        va='center',
                        color=[0.5]*3,
                        alpha=1/2,
                        size=0.8*plot_preferences.plt.rcParams['font.size'],#'smaller',
                        rotation=angle,
                        rotation_mode='anchor',
                        backgroundcolor=[1]*3,
                        zorder=0)



#%%
def _get_smallest_larger(seq, value):
    """
    Identifies the smallest value in a sequence that is greater than or equal to a
    specified value, along with its index.

    This function searches through a list or dictionary to find the smallest value
    that is greater than or equal to the given value. It returns both the index (or
    key in the case of a dictionary) and the value itself.

    Parameters:
    - seq (list or dict): The sequence to be searched. It can be either a list or a
      dictionary. If it's a dictionary, the function operates on its values.
    - value (numeric): The value against which the elements in the sequence are
      compared.

    Returns:
    - tuple: A tuple containing the index (or key) and the value of the smallest
      entry greater than or equal to the specified value. Returns None if no such
      entry exists.

    Note:
    - If the sequence is a dictionary, the function returns the key associated with
      the found value, not the positional index.
    - For a list, the function returns the index of the found value.
    - In case no value in the sequence is greater than or equal to the specified
      value, the function returns None.
    """
    try:
        if type(seq) == dict:
            tmp_seq = list(seq.values()) 
            val_min = min(x for x in tmp_seq if x >= value)
            id_min = tmp_seq.index(val_min)
            id_min = list(seq.keys())[id_min]
        else:
            val_min = min(x for x in seq if x >= value)
            id_min = seq.index(val_min)
        return id_min, val_min
    except ValueError:
        return None
    
def _get_largest_smaller(seq, value):
    """
    Identifies the largest value in a sequence that is less than or equal to a
    specified value, along with its index.

    This function searches a list or dictionary for the largest value that is less
    than or equal to the given value. It returns both the index (or key, in the case
    of a dictionary) and the value itself.

    Parameters:
    - seq (list or dict): The sequence to search through. Can be either a list or a
      dictionary. If it's a dictionary, the function operates on its values.
    - value (numeric): The value to compare against the elements in the sequence.

    Returns:
    - tuple: A tuple containing the index (or key) and the value of the largest
      entry less than or equal to the specified value. Returns None if no such
      entry exists.

    Note:
    - For a dictionary, the function returns the key associated with the found
      value, not the positional index.
    - For a list, the function returns the index of the found value.
    - If no value in the sequence is less than or equal to the given value, the
      function returns None.
    """
    try:
        if type(seq) == dict:
            tmp_seq = list(seq.values()) 
            val_max = max(x for x in tmp_seq if x <= value)
            id_max = tmp_seq.index(val_max)
            id_max = list(seq.keys())[id_max]
        else:
            val_max = max(x for x in seq if x <= value)
            id_max = seq.index(val_max)
        return id_max, val_max
    except ValueError:
        return None
    


#%%
def interp_dp(tmp_res,value,
              calc_EPsys=False,prefix_str='cell_model[1,1].'):
    """
    Interpolates data points in a simulation result set based on a specified
    energy value.

    This function interpolates data points in the `tmp_res` dictionary to
    determine values corresponding to a given energy value. It handles different
    cases based on whether the value is below the minimum, above the maximum, or
    within the dataset's range. This function supports both system-level and
    cell-level data interpolation.

    Parameters:
    - tmp_res (dict): A dictionary containing simulation results. Each key
      corresponds to a different parameter, with values as data arrays.
    - value (float): The energy value for which data points are interpolated.
    - calc_EPsys (bool, optional): If True, interpolates system-level data.
      Defaults to False, which interpolates cell-level data.
    - prefix_str (str, optional): The string identifying the specific sub-model
      in the result keys. Defaults to 'cell_model[1,1]'.

    Returns:
    - tuple: A tuple containing two elements:
        1. A dictionary (`dp`) with interpolated values for each parameter in
           `tmp_res`.
        2. An integer (`case_dp`) indicating the interpolation case: 0 for a
           value less than the minimum, 1 for a value greater than the maximum,
           and 2 for a value within the range.

    Note:
    - The function utilizes `_get_largest_smaller` and `_get_smallest_larger` to
      find adjacent entries for accurate interpolation.
    - It handles different data types in `tmp_res`, interpolating only numerical
      arrays and excluding specific keys as specified in
      `tmp_res['_settings']['ParameterValues']`.
    - The function returns the first or last data point if the specified energy
      value is outside the dataset's range.
    """
    dp = {}
    try:
        if calc_EPsys:
            tmp_E = tmp_res['E_sys'].tolist()
        else:
            # tmp_E = tmp_res['E_cell[1,1]'].tolist()
            tmp_E = tmp_res[prefix_str+'E_cell'].tolist()
    
        if value <= tmp_E[0]: # case 0
            case_dp = 0
            for key in tmp_res:
                try:
                    dp[key] = tmp_res[key][0]
                except:
                    0
        elif value >= tmp_E[-1]: # case 1
            case_dp = 1
            for key in tmp_res:
                try:
                    dp[key] = tmp_res[key][-1]
                except:
                    0
            # dp = {}
        else: # case 2
            case_dp = 2
            id_min,valx_min = _get_largest_smaller(tmp_E,value)
            id_max,valx_max = _get_smallest_larger(tmp_E,value)
            for key in tmp_res:
                try:
                    if isinstance(tmp_res[key], np.ndarray) and key not in tmp_res['_settings']['ParameterValues']:
                        valy_min = tmp_res[key][id_min]
                        valy_max = tmp_res[key][id_max]
                        dp[key] = np.interp(value,[valx_min,valx_max],[valy_min,valy_max])
                    else:
                        dp[key] = tmp_res[key][-1]
                except:
                    0
    except:
        0
    return dp, case_dp



#%%
def calc_P_UI(tmp_P,tmp_E,Umin,Imax):
    """
    Calculates the specific switching point in power and energy on a Ragone plot
    based on voltage and current limits.

    This function determines the switching point (P_UI, E_UI) on a Ragone plot
    where the battery system transitions from being voltage-limited ('U' for
    minimum cell power) to current-limited ('I' for maximum current limit). It
    calculates the power at this switching point using the minimum voltage (Umin)
    and maximum current (Imax) limits, and interpolates the corresponding energy
    value.

    Parameters:
    - tmp_P (list or array): An array of power values, typically representing a
      power profile.
    - tmp_E (list or array): An array of corresponding energy values, typically
      representing an energy profile.
    - Umin (float): The minimum voltage limit.
    - Imax (float): The maximum current limit.

    Returns:
    - tuple: A tuple containing two elements:
        1. Power at the switching point (P_UI).
        2. Interpolated energy value at this power level (E_UI).

    Note:
    - The function employs '_get_largest_smaller' and '_get_smallest_larger' to
      accurately find power values close to P_UI for interpolation.
    - It then interpolates these values to find the corresponding energy value in
      'tmp_E'.
    - This calculation is particularly relevant in battery system analyses, where
      understanding the transition between voltage and current limits is crucial.
    """
    P_UI = Umin * Imax

    id_min,valx_min = _get_largest_smaller(tmp_P,P_UI)
    id_max,valx_max = _get_smallest_larger(tmp_P,P_UI)
    valy_min = tmp_E[id_min]
    valy_max = tmp_E[id_max]
    
    E_UI = np.interp(P_UI,[valx_min,valx_max],[valy_min,valy_max])
    
    return P_UI, E_UI



#%%
def calc_intersection_curve(g, res, set_EP, zfig=2, nSoC=0, print_error=True, 
                            prefix_str='cell_model[1,1].'):
    """
    Calculate the intersection curve for a specified EP (energy-per-power) set
    point on a Ragone plot.

    This function computes the intersection points between a linear EP line and 
    battery system performance curves at a given state of charge (SoC). It 
    interpolates data to construct an intersection curve and find boundary points.

    Parameters:
    - g (dict): Dictionary with matplotlib graph objects.
    - res (list of dicts): List containing dictionaries of simulation results.
    - set_EP (float): EP set point for calculating the intersection curve.
    - zfig (int, optional): The figure index for plotting. Defaults to 2.
    - nSoC (int, optional): State of charge for calculation. Defaults to 0.
    - print_error (bool, optional): Flag to print error messages. Defaults to True.
    - prefix_str (str, optional): Substring identifying specific sub-model in
      result keys. Defaults to 'cell_model[1,1]'.

    Returns:
    - dict: Dictionary with intersection curve data. Keys correspond to parameters,
      values are data arrays.

    Notes:
    - Utilizes helper functions for boundary finding and interpolation.
    - Intersection points calculation based on Stack Overflow methodology.
    - Assumes specific graph objects in 'g' for plot limits and scales.
    - Error handling includes optional printing of messages when a design point
      can't be interpolated.
    """
    
    # Calculate EP-Line of Application
    def linear(m,b):
        return lambda x: m*x+b
    def inv_linear(m,b):
        return lambda y: (y-b)/m
    
    # Find Line-Line Intersection (https://stackoverflow.com/questions/20677795/how-do-i-compute-the-intersection-point-of-two-lines)
    def findIntersection(x1,y1,x2,y2,x3,y3,x4,y4):
        px= ( (x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4) ) / ( (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4) ) 
        py= ( (x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4) ) / ( (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4) )
        return [px, py]
    
    
    f_EP = linear(set_EP,0)
    inv_f_EP = inv_linear(set_EP,0)
    
    loc_res = {}
    loc_Plist = {}

    # find all entries in res at nSoC
    for rn,_ in enumerate(res):
        if res[rn]['_nSoC'] == nSoC:
            loc_res[rn] = res[rn]
            loc_Plist[rn] = res[rn][prefix_str+'P_cell'][0]
  
    # interpolation between P_set
    pp = g['ax_{}_{}'.format(zfig,0)].get_xlim()
    ee = g['ax_{}_{}'.format(zfig,0)].get_ylim()
    pp = np.linspace(pp[0],pp[-1] if f_EP(pp[-1])<=ee[-1] else inv_f_EP(ee[-1]),int(5e2))
    
    # evaluate P_set
    # pp= np.array([loc_Plist[n] for n,i in enumerate(loc_Plist)])
    
    dp = {}
    dp_intersection = {}
    for n,i in enumerate(pp):      
        set_Psys = i
        set_Esys = f_EP(set_Psys)
        dp[n] = {}
        try:
            rn_l,val_l = _get_largest_smaller(loc_Plist, set_Psys) # --> left boundary 
            res_l = loc_res[rn_l]
            dp_l,case_dp_l = interp_dp(res_l,set_Esys,prefix_str=prefix_str)
            # print('[Info] lower - P = {:.3f}, case_dp = {}'.format(set_Psys,case_dp_l))
            
            if case_dp_l == 2:
                if set_Psys in loc_Plist.values():
                    dp[n] = dp_l
                else:
                    rn_u,val_u = _get_smallest_larger(loc_Plist, set_Psys) # --> right boundary
                    res_u = loc_res[rn_u]
                    dp_u,case_dp_u = interp_dp(res_u,set_Esys,prefix_str=prefix_str) #if set_Esys < tmp_EE[-1] else {}
                    # print('[Info] upper - P = {:.3f}, case_dp = {}'.format(set_Psys,case_dp_u))
                    
                    EP_intersect = findIntersection(val_l,res_l[prefix_str+'E_cell'][-1],
                                                    val_u,res_u[prefix_str+'E_cell'][-1],
                                                    val_l,f_EP(val_l),
                                                    val_u,f_EP(val_u))
                    if set_Psys <= EP_intersect[0]:
                        # interpolation between dp_l & dp_u
                        for key in dp_u:
                            if key in dp_l:
                                dp[n][key] = np.interp(set_Psys,[dp_l[prefix_str+'P_cell'],dp_u[prefix_str+'P_cell']],[dp_l[key],dp_u[key]])
                            else:
                                continue
                    
            # Error --> print
            if not any(dp[n]):
                raise ValueError
                        
        except:
            if print_error:
                print('[Info] Design Point n = {} not available: P/E = {:.3f}/{:.3f}'.format(n,set_Psys,set_Esys))
            
        # dp_intersectionectory
        if any(dp[n]) and not any(dp_intersection): # start
            for key in dp[n]:
                dp_intersection[key] = np.append([np.nan]*n,dp[n][key])
        else:
            for key in dp_intersection:
                dp_intersection[key] = np.append(dp_intersection[key],dp[n][key] if any(dp[n]) else np.nan)

    dp_intersection['valx_P'] = pp
    dp_intersection['valy_E'] = f_EP(dp_intersection['valx_P'])
    
    return dp_intersection


#%%
def _set_intersection_limits(g, EPline_x, EPline_y, zfig, zax, replot=False, **kwargs):
    """
    Set intersection limits for plotting based on specified EP line coordinates.

    This function adjusts plot limits and optionally replots intersection markers 
    based on the maximum x-coordinate from a non-NaN value in the EP line. It 
    interprets the intersection of the EP line with a boundary line, setting plot 
    data accordingly.

    Parameters:
    - g (dict): Dictionary with matplotlib graph objects.
    - EPline_x (list): x-coordinates of the EP line.
    - EPline_y (list): y-coordinates of the EP line.
    - zfig (int): Figure index used for plotting.
    - zax (int): Axis index used for plotting.
    - replot (bool, optional): Flag to control replotting. Defaults to False.
    - **kwargs: Additional keyword arguments. Can include 'x' to specify a custom
      x-coordinate for intersection.

    Returns:
    - tuple: A tuple (x, yint) representing the x-coordinate for intersection and 
      interpolated y-coordinate.

    Notes:
    - Utilizes `np.interp` for interpolation.
    - Handles replotting and plot adjustments based on 'replot' flag.
    - In case of exceptions, the function returns zero without error.
    """
    id_max = max([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
    x = kwargs['x'] if 'x' in kwargs else EPline_x[id_max]
    yint = np.interp(x,EPline_x,EPline_y)
    try:
        farbe = plot_preferences.plot_pre()[0]
        if replot:
            g['pl_{}_{}_{}'.format(zfig,zax,110)][0].set_data([x]*2,[-1e20,1e20])
            g['pl_{}_{}_{}'.format(zfig,zax,111)][0].set_data(x,yint)
        else:
            g['pl_{}_{}_{}'.format(zfig,zax,110)] = g['ax_{}_{}'.format(zfig,zax)].plot([x]*2,[-1e20,1e20],color=farbe[2],ls='dotted',lw=1,zorder=100)
            g['pl_{}_{}_{}'.format(zfig,zax,111)] = g['ax_{}_{}'.format(zfig,zax)].plot(x,yint,color=farbe[2],ls='None',marker='o',markersize=4.5,zorder=100)
        
        if zax == 1:
            g['an_{}_{}_{}'.format(zfig,zax,110)].xy = (x,g['ax_{}_{}'.format(zfig,zax)].get_ylim()[1])
    except:
        0
    return x, yint



#%%
def hyperbel(m):
    """
    Creates a hyperbolic function based on a specified multiplier.
    
    This function generates a hyperbolic function of the form y = m/x, where 'm' is the provided multiplier and 'x' is the variable. The generated function can then be used to calculate y for any given x.
    
    Parameters:
    - m (float): The multiplier used in the hyperbolic function.
    
    Returns:
    - function: A lambda function representing the hyperbolic function y = m/x.
    
    Note:
    - The returned function is a lambda function, which can be used directly to compute y for any given x.
    - This function is often used in scientific and engineering calculations where a hyperbolic relationship between two quantities is needed.
    """
    return lambda x: m*x**(-1)



#%%
def plot_dss(g, c, zfig, zax, print_opt=False, **kwargs):
    """
    Plot the Design Solution Space (DSS) with constraints for an energy storage
    system.

    This function visualizes the DSS for an energy storage system by plotting 
    constraint satisfaction problem solutions graphically. It computes and plots 
    constraint lines (or boundary lines), including hyperbolic and linear 
    constraints, to demarcate feasible design regions. The function also marks 
    optimal design points and, optionally, prints detailed design information. 
    Future iterations of the function will integrate a real optimization problem.

    Parameters:
    - g (dict): Dictionary with matplotlib graph objects.
    - c (dict): Dictionary containing cell and system design parameters.
    - zfig (int): Figure index for plotting.
    - zax (int): Axis index for plotting.
    - print_opt (bool, optional): If True, prints optimal design details. 
      Defaults to False.
    - **kwargs: Additional keyword arguments for plot customization.

    Returns:
    - dict: Dictionary containing location data of plotted design points.

    Notes:
    - Uses numpy for array operations and matplotlib for plotting.
    - Defines and utilizes hyperbolic functions to plot design constraints.
    - Handles both initial plotting and updating existing plots.
    - Dynamically adjusts annotations and color bars based on design data.
    - Currently solves the optimization problem graphically, with plans to
      integrate a computational optimization approach in the future.
    """
    xx = np.linspace(1e-4,1e4,int(1e4))
    yy = xx
    
    xmax,_,f_xmax = calc_design([xx.min(),xx.max()],[c['pec']['Udc_max']]*2,m=c['cell']['Udc_max'],show_plt=False,show_annotation=False)
    xmin,_,f_xmin = calc_design([xx.min(),xx.max()],[c['pec']['Udc_min']]*2,m=c['cell']['Udc_min'],show_plt=False,show_annotation=False)
    ymax,_,f_ymax = calc_design([yy.min(),yy.max()],[c['pec']['n']*c['pec']['Idc_max']]*2,m=c['cell']['Idc_DIS_max'],show_plt=False,show_annotation=False,swapxy=True)
     
    xmax = xmax[0]
    xmin = xmin[0]
    ymax = ymax[0]
     
    f_hyp1 = hyperbel((c['req']['Pac_req']/c['pec']['eta_op'])/c['cell']['Pdc_DIS_max'])
    f_hyp2 = hyperbel(c['pec']['n']*c['pec']['Pdc_max']/c['cell']['Pdc_DIS_max'])

    xxmax = np.floor(xmax)
    xxmin = np.ceil(xmin)
    yymax = np.floor(ymax)

    xmax_P = f_hyp2(ymax)
    xxmax_P = np.floor(xmax_P)
    ymax_P = f_hyp1(xmax_P)
    yymax_P = np.ceil(ymax_P)
    
    loc = {}
    loc['x'] = np.linspace(max(xmin,f_hyp1(ymax)),xmax,int(1e4)) if max(xmin,f_hyp1(ymax)) <= xmax else np.array([])
    loc['y'] = f_hyp1(loc['x'])
    loc['y_P'] = f_hyp2(loc['x'])
    loc['xx'] = [*range(int(np.ceil(loc['x'][0])),int(np.floor(loc['x'][-1])+1))] if loc['x'].any() else [] 
    loc['yy'] = [*range(int(np.ceil(loc['y'][-1])),int(np.floor(loc['y'][0])+1))] if loc['y'].any() else []
    loc['yy_P'] = [*range(int(np.ceil(loc['y_P'][-1])),int(np.floor(loc['y_P'][0])+1))] if loc['y_P'].any() else []

    points = np.array([loc['x'],loc['y']]).T
    points_P = np.array([loc['x'],loc['y_P']]).T
    if points.any(): # = max(xmin,f_hyp1(ymax)) <= xmax
        if xmax_P < f_hyp1(ymax): #case1 = errorPmax
            points2 = np.array([[0,0],[0,0]])
            points = np.array([[0,0],[0,0]])
        elif f_hyp1(ymax) <= xmin and xmax_P <= xmin: #case2
            points2 = np.vstack([points,points_P[::-1]])
            points = np.vstack([points,[points[-1,0],0],[points[0,0],0],[points[0,0],points[-1,1]],[0,points[-1,1]],[0,points_P[0,1]],[points[0,0],points_P[0,1]]])
        elif f_hyp1(ymax) <= xmin and xmin <= xmax_P <= xmax: #case3
            points_P = points_P[points_P[:,-1]<ymax] if points_P.any() else points_P
            points2 = np.vstack([points,points_P[::-1],[f_hyp2(ymax),ymax],[xmin,ymax]])
            points = np.vstack([points,[points[-1,0],0],[points[0,0],0],[points[0,0],points[-1,1]],[0,points[-1,1]],[0,ymax],[points[0,0],ymax]])
        elif f_hyp1(ymax) <= xmin and xmax <= xmax_P: #case4
            points2 = np.vstack([points,[points[-1,0],ymax],[points[0,0],ymax]])
            points = np.vstack([points,[points[-1,0],0],[points[0,0],0],[points[0,0],points[-1,1]],[0,points[-1,1]],[0,ymax],[points[0,0],ymax]])
        elif xmin <= f_hyp1(ymax) <= xmax and xmax <= xmax_P: #case5
            points2 = np.vstack([points,[points[-1,0],points[0,1]]])
            points = np.vstack([points,[points[-1,0],0],[points[0,0],0],[points[0,0],points[-1,1]],[0,points[-1,1]],[0,points[0,1]]])
        else: #case6
            points_P = points_P[points_P[:,-1]<points[0,-1]] if points_P.any() else points_P
            points2 = np.vstack([points,points_P[::-1],[f_hyp2(points[0,-1]),points[0,-1]]])
            points = np.vstack([points,[points[-1,0],0],[points[0,0],0],[points[0,0],points[-1,1]],[0,points[-1,1]],[0,points[0,1]]])
        loc['yy'] = [np.ceil(points2[:,1].min()).astype(int),np.floor(points2[:,1].max()).astype(int)]
    else: #case0
        points2 = np.array([[0,0],[0,0]])
        points = np.array([[0,0],[0,0]])


    # Create Plots
    farbe = plot_preferences.plot_pre()[0]
    cmap = tol_cmap('rainbow_PuRd').reversed()
    cmap.set_bad('None')
    
    if 'pl_{}_{}_{}'.format(zfig,zax,0) not in g:
        g['pl_{}_{}_{}'.format(zfig,zax,0)] = g['ax_{}_{}'.format(zfig,zax)].plot(xx,f_hyp1(xx),color=farbe[0],zorder=3)
        g['pl_{}_{}_{}'.format(zfig,zax,1)] = g['ax_{}_{}'.format(zfig,zax)].plot(xx,f_hyp2(xx),color=farbe[7],zorder=3)
        
        # try:
        #     g['pl_{}_{}_{}'.format(zfig,zax,2)] = g['ax_{}_{}'.format(zfig,zax)].plot([0,xmin],[f_hyp1(xmin)]*2,color='k',ls='dotted',marker='x')
        #     g['pl_{}_{}_{}'.format(zfig,zax,3)] = g['ax_{}_{}'.format(zfig,zax)].plot([0,xmin],[f_hyp2(xmin)]*2,color='k',ls='dotted',marker='x')
        #     g['pl_{}_{}_{}'.format(zfig,zax,4)] = g['ax_{}_{}'.format(zfig,zax)].plot([0,xmax],[f_hyp1(xmax)]*2,color='k',ls='dotted',marker='x')
        #     g['pl_{}_{}_{}'.format(zfig,zax,5)] = g['ax_{}_{}'.format(zfig,zax)].plot([0,xmax],[f_hyp2(xmax)]*2,color='k',ls='dotted',marker='x')
        #     g['pl_{}_{}_{}'.format(zfig,zax,6)] = g['ax_{}_{}'.format(zfig,zax)].plot([f_hyp1(ymax)]*2,[0,ymax],color='k',ls='dotted',marker='x')
        #     g['pl_{}_{}_{}'.format(zfig,zax,7)] = g['ax_{}_{}'.format(zfig,zax)].plot([f_hyp2(ymax)]*2,[0,ymax],color='k',ls='dotted',marker='x')
        # except:
        #     0
        
        g['pl_{}_{}_{}'.format(zfig,zax,10)] = g['ax_{}_{}'.format(zfig,zax)].plot([xmax,xmax],[0, 1e20],color=farbe[4],marker='^',zorder=3) 
        g['pl_{}_{}_{}'.format(zfig,zax,11)] = g['ax_{}_{}'.format(zfig,zax)].plot([xmin,xmin],[0, 1e20],color=farbe[5],marker='^',zorder=3) 
        g['pl_{}_{}_{}'.format(zfig,zax,12)] = g['ax_{}_{}'.format(zfig,zax)].plot([0, 1e20],[ymax,ymax],color=farbe[3],marker='v',zorder=3)      
        
        g['pl_{}_{}_{}'.format(zfig,zax,100)] = [Polygon(points, closed=True),Polygon(points2,closed=True)]
        polycollection = UpdatablePatchCollection(g['pl_{}_{}_{}'.format(zfig,zax,100)],color=[farbe[7],farbe[1]],zorder=0,alpha=[0,1/6],edgecolor=None)
        g['pl_{}_{}_{}'.format(zfig,zax,101)] = g['ax_{}_{}'.format(zfig,zax)].add_collection(polycollection)  
            
    else:
        g['pl_{}_{}_{}'.format(zfig,zax,0)][0].set_data(xx,f_hyp1(xx))
        g['pl_{}_{}_{}'.format(zfig,zax,1)][0].set_data(xx,f_hyp2(xx))
        
        try:
            g['pl_{}_{}_{}'.format(zfig,zax,2)][0].set_data([0,xmin],[f_hyp1(xmin)]*2)
            g['pl_{}_{}_{}'.format(zfig,zax,3)][0].set_data([0,xmin],[f_hyp2(xmin)]*2)
            g['pl_{}_{}_{}'.format(zfig,zax,4)][0].set_data([0,xmax],[f_hyp1(xmax)]*2)
            g['pl_{}_{}_{}'.format(zfig,zax,5)][0].set_data([0,xmax],[f_hyp2(xmax)]*2)
            g['pl_{}_{}_{}'.format(zfig,zax,6)][0].set_data([f_hyp1(ymax)]*2,[0,ymax])
            g['pl_{}_{}_{}'.format(zfig,zax,7)][0].set_data([f_hyp2(ymax)]*2,[0,ymax])
        except:
            0
   
        g['pl_{}_{}_{}'.format(zfig,zax,10)][0].set_data([xmax,xmax],[0,1e8])
        g['pl_{}_{}_{}'.format(zfig,zax,11)][0].set_data([xmin,xmin],[0,1e8])
        g['pl_{}_{}_{}'.format(zfig,zax,12)][0].set_data([0,1e8],[ymax,ymax])
        
        g['pl_{}_{}_{}'.format(zfig,zax,100)][0].xy = points
        g['pl_{}_{}_{}'.format(zfig,zax,100)][1].xy = points2


    minmax_range = lambda i: range(min(i),max(i)+1)
    # loc['yy'] = [np.ceil(points2[:,1].min()).astype(int),np.floor(points2[:,1].max()).astype(int)] if points.any() else []
    loc['xy_list'] = [[x,y,x*y if g['pl_{}_{}_{}'.format(zfig,zax,100)][1].contains_point((x,y)) else np.nan] for x in loc['xx'] for y in (minmax_range(loc['yy']+loc['yy_P']) if f_hyp1(ymax) <= xmax_P else [])]
    loc['xy_opt'] = min(loc['xy_list'], key = lambda t:t[2] if np.mod(t[0],c['mod']['xS'])==0 and np.mod(t[1],c['mod']['yP'])==0 and ~np.isnan(t[2]) else np.inf) if any(loc['xy_list']) else (np.nan,)*3
    loc['xy_opt'] = loc['xy_opt'] if ~np.isnan(loc['xy_opt'][2]) else (np.nan,)*3   # ~np.isnan(t[2]) funktioniert nicht in lambda t!


    # scatter
    sclist = np.array(loc['xy_list'])
    if sclist.size:
        # whole-number / modularity
        scx_mod = (np.mod(sclist[:,0],c['mod']['xS'])==0)
        scy_mod = (np.mod(sclist[:,1],c['mod']['yP'])==0)
        scxy_mod = scx_mod*scy_mod
        
        scx,scy,scz = sclist[:,0],sclist[:,1],sclist[:,2]
        
        if 'pl_{}_{}_{}'.format(zfig,zax,20) not in g:
            g['pl_{}_{}_{}'.format(zfig,zax,20)] = g['ax_{}_{}'.format(zfig,zax)].scatter(scx[~np.isnan(scz)*~scxy_mod],scy[~np.isnan(scz)*~scxy_mod],c='0',s=1,marker='.',zorder=10)
            g['pl_{}_{}_{}'.format(zfig,zax,21)] = g['ax_{}_{}'.format(zfig,zax)].scatter(scx[scxy_mod],scy[scxy_mod],c=scz[scxy_mod],cmap=cmap,s=75,marker='.',zorder=10) #cm.get_cmap('plasma_r')
        else:
            g['pl_{}_{}_{}'.format(zfig,zax,20)].set_offsets(np.vstack([scx[~np.isnan(scz)*~scxy_mod],scy[~np.isnan(scz)*~scxy_mod]]).T)
            g['pl_{}_{}_{}'.format(zfig,zax,21)].set_offsets(np.vstack([scx[scxy_mod],scy[scxy_mod]]).T)
            g['pl_{}_{}_{}'.format(zfig,zax,21)].set_array(scz[scxy_mod])
        
        if all(~np.isnan(scz)==False):
            g['pl_{}_{}_{}'.format(zfig,zax,21)].set_clim(0,1)
        else:
            g['pl_{}_{}_{}'.format(zfig,zax,21)].set_clim(plot_preferences._set_lim(min(scz[~np.isnan(scz)]),max(scz[~np.isnan(scz)]),0.005))
        
        loc['xy_list'] = sclist[~np.isnan(sclist[:,2])].astype(int).tolist()
    
    else:
        scxy_mod = np.array([False])
        scx,scy,scz = [],[],[]
        
        if 'pl_{}_{}_{}'.format(zfig,zax,20) not in g:
            g['pl_{}_{}_{}'.format(zfig,zax,20)] = g['ax_{}_{}'.format(zfig,zax)].scatter([],[],c='0',s=1,marker='.',zorder=10)
            g['pl_{}_{}_{}'.format(zfig,zax,21)] = g['ax_{}_{}'.format(zfig,zax)].scatter([],[],c=[],cmap=cmap,s=75,marker='.',zorder=10)
        else:
            g['pl_{}_{}_{}'.format(zfig,zax,20)].set_offsets(np.vstack([scx,scy]).T)
            g['pl_{}_{}_{}'.format(zfig,zax,21)].set_offsets(np.vstack([scx,scy]).T)

        
    
    if 'cbar_{}_{}'.format(zfig,zax) not in g:
        # scz_norm = [] if not sclist.size else (scz - np.min(scz))/np.ptp(scz)    
        g['cbar_{}_{}'.format(zfig,zax)] = g['fig_{}'.format(zfig)].colorbar(g['pl_{}_{}_{}'.format(zfig,zax,21)],ax=g['ax_{}_{}'.format(zfig,zax)],location='top')
        # g['cbar_{}_{}'.format(zfig,zax)].ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
        g['cbar_{}_{}'.format(zfig,zax)].ax.xaxis.set_major_locator(MaxNLocator(nbins=10,integer=True)) 
        g['cbar_{}_{}'.format(zfig,zax)].set_label(r'Number of cells, $n_{x_\mathrm{S}y_\mathrm{P}}$')
        
        g['cbar_{}_{}'.format(zfig,zax+1)] = g['cbar_{}_{}'.format(zfig,zax)].ax.secondary_xaxis('bottom')
        g['cbar_{}_{}'.format(zfig,zax+1)].zorder-=1
        

    if ~np.isnan(loc['xy_opt'][2]):
        clim = g['pl_{}_{}_{}'.format(zfig,zax,21)].get_clim()
        opt_share = (loc['xy_opt'][2]-clim[0])/(clim[1]-clim[0])
        if 'pl_{}_{}_{}'.format(zfig,zax,22) not in g:
            g['pl_{}_{}_{}'.format(zfig,zax,22)] = g['ax_{}_{}'.format(zfig,zax)].plot(loc['xy_opt'][0],loc['xy_opt'][1],c='k',marker='o',ms=7,mfc=cmap(opt_share),mew=2,zorder=10) # ms=6.5
        else:
            g['pl_{}_{}_{}'.format(zfig,zax,22)][0].set_data(loc['xy_opt'][0],loc['xy_opt'][1])
            g['pl_{}_{}_{}'.format(zfig,zax,22)][0].set_mfc(cmap(opt_share))

        if 'cpl_{}_{}_{}'.format(zfig,zax,0) not in g:
            g['cpl_{}_{}_{}'.format(zfig,zax,0)] = g['cbar_{}_{}'.format(zfig,zax)].ax.axvline(x=loc['xy_opt'][2], c='k')
        else:
            g['cpl_{}_{}_{}'.format(zfig,zax,0)].set_xdata([loc['xy_opt'][2]]*2)  
  
        g['cbar_{}_{}'.format(zfig,zax)].ax.set_visible(True)
        # g['cbar_{}_{}'.format(z,2)].ax.set_visible(True)
        g['cbar_{}_{}'.format(zfig,zax+1)].set_ticks([loc['xy_opt'][2]])
        g['cbar_{}_{}'.format(zfig,zax+1)].set_xticklabels([r'$n_{{x_{{\mathrm{{S}}}}y_{{\mathrm{{P}}}},\mathrm{{opt}}}} = {}$'.format(loc['xy_opt'][2])])

    else:
        if 'pl_{}_{}_{}'.format(zfig,zax,22) not in g:
            # clim = g['pl_{}_{}_{}'.format(zfig,zax,21)].get_clim()
            # opt_share = (loc['xy_opt'][2]-clim[0])/(clim[1]-clim[0])
            # g['pl_{}_{}_{}'.format(zfig,zax,22)] = g['ax_{}_{}'.format(zfig,zax)].plot(loc['xy_opt'][0],loc['xy_opt'][1],c='k',marker='o',ms=6.5,mfc='None',mew=2,zorder=10)
            g['pl_{}_{}_{}'.format(zfig,zax,22)] = g['ax_{}_{}'.format(zfig,zax)].plot([],[],c='k',marker='o',ms=6.5,mfc='None',mew=2,zorder=10)
        else:
            g['pl_{}_{}_{}'.format(zfig,zax,22)][0].set_data([],[])

        if 'cpl_{}_{}_{}'.format(zfig,zax,0) not in g:
            g['cpl_{}_{}_{}'.format(zfig,zax,0)]= g['cbar_{}_{}'.format(zfig,zax)].ax.axvline(x=np.nan, c='k')
        else:
            g['cpl_{}_{}_{}'.format(zfig,zax,0)].set_xdata([np.nan]*2)
            
        g['cbar_{}_{}'.format(zfig,zax)].ax.set_visible(False)
        # g['cbar_{}_{}'.format(z,2)].ax.set_visible(False)
        g['cbar_{}_{}'.format(zfig,zax+1)].set_ticks([np.nan])
        g['cbar_{}_{}'.format(zfig,zax+1)].set_xticklabels([''])
    
    
    # Limit
    xlist = [xmin, xmax, f_hyp1(ymax), xmax_P]        
    xlim = kwargs['xlim'] if ('xlim' in kwargs and kwargs['xlim']) else plot_preferences._set_lim(min(xlist),max(xlist),0.2)
    g['ax_{}_{}'.format(zfig,zax)].set_xlim(xlim)
    g['ax_{}_{}'.format(zfig,zax)].xaxis.set_major_locator(MaxNLocator(integer=True))
    
    ylist = [f_hyp1(xmax), f_hyp1(xmin), ymax, f_hyp1(xmax_P)]
    ylim = kwargs['ylim'] if ('ylim' in kwargs and kwargs['ylim']) else plot_preferences._set_lim(min(ylist),max(ylist),0.2) 
    g['ax_{}_{}'.format(zfig,zax)].set_ylim(ylim)
    g['ax_{}_{}'.format(zfig,zax)].yaxis.set_major_locator(MaxNLocator(integer=True))

    g['ax_{}_{}'.format(zfig,zax)].set_aspect(1./g['ax_{}_{}'.format(zfig,zax)].get_data_ratio())
        
        
    # Annotation
    bbox_args = dict(boxstyle="round",color='None',fc='1',lw=0)
    if 'an_{}_{}_{}'.format(zfig,zax,0) not in g:
        g['an_{}_{}_{}'.format(zfig,zax,0)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'$\textbf{\textsf{C}}_{P_\mathrm{req}}$',
                                                      (f_hyp1(ylim[0]+np.diff(ylim)*0.96),ylim[0]+np.diff(ylim)*0.96),
                                                      textcoords="offset points",
                                                      xytext=(0,0),
                                                      c=farbe[0],
                                                      fontsize=11,
                                                      # size='smaller',
                                                      ha='center',
                                                      va='center',
                                                      bbox=bbox_args,
                                                      zorder=100)
        g['an_{}_{}_{}'.format(zfig,zax,1)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'$\textbf{\textsf{C}}_{P_\mathrm{max}}$',
                                                      (f_hyp2(ylim[0]+np.diff(ylim)*0.96),ylim[0]+np.diff(ylim)*0.96),
                                                      textcoords="offset points",
                                                      xytext=(0,0),
                                                      c=farbe[7],
                                                      fontsize=11,
                                                      # size='smaller',
                                                      ha='center',
                                                      va='center',
                                                      bbox=bbox_args,
                                                      zorder=100)
        g['an_{}_{}_{}'.format(zfig,zax,10)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'$\textbf{\textsf{C}}_{U_\mathrm{max}}$',
                                                      (xmax,ylim[0]+np.diff(ylim)*0.04),
                                                      textcoords="offset points",
                                                      xytext=(0,0),
                                                      c=farbe[4],
                                                      fontsize=11,
                                                      # size='smaller',
                                                      ha='center',
                                                      va='center',
                                                      bbox=bbox_args,
                                                      zorder=100)
        g['an_{}_{}_{}'.format(zfig,zax,11)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'$\textbf{\textsf{C}}_{U_\mathrm{min}}$',
                                                      (xmin,ylim[0]+np.diff(ylim)*0.04),
                                                      textcoords="offset points",
                                                      xytext=(0,0),
                                                      c=farbe[5],
                                                      fontsize=11,
                                                      # size='smaller',
                                                      ha='center',
                                                      va='center',
                                                      bbox=bbox_args,
                                                      zorder=100)
        g['an_{}_{}_{}'.format(zfig,zax,12)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'$\textbf{\textsf{C}}_{I_\mathrm{max}}$',
                                                      (xlim[0]+np.diff(xlim)*0.08,ymax),
                                                      textcoords="offset points",
                                                      xytext=(0,0),
                                                      c=farbe[3],
                                                      fontsize=11,
                                                      # size='smaller',
                                                      ha='center',
                                                      va='center',
                                                      bbox=bbox_args,
                                                      zorder=100)
    else:    
        g['an_{}_{}_{}'.format(zfig,zax,0)].xy = (f_hyp1(ylim[0]+np.diff(ylim)*0.96),ylim[0]+np.diff(ylim)*0.96)
        g['an_{}_{}_{}'.format(zfig,zax,1)].xy = (f_hyp2(ylim[0]+np.diff(ylim)*0.96),ylim[0]+np.diff(ylim)*0.96)
        g['an_{}_{}_{}'.format(zfig,zax,10)].xy = (xmax,ylim[0]+np.diff(ylim)*0.04)
        g['an_{}_{}_{}'.format(zfig,zax,11)].xy = (xmin,ylim[0]+np.diff(ylim)*0.04)
        g['an_{}_{}_{}'.format(zfig,zax,12)].xy = (xlim[0]+np.diff(xlim)*0.08,ymax)

    # # Hide Ticklabels
    # plot_preferences.hide_ticklabels(g['ax_{}_{}'.format(zfig,zax)].xaxis)
    # plot_preferences.hide_ticklabels(g['ax_{}_{}'.format(zfig,zax)].yaxis)
    # plot_preferences.hide_ticklabels(g['cbar_{}_{}'.format(zfig,zax)].ax.yaxis)

    g['fig_{}'.format(zfig)].canvas.draw()
    
    if print_opt:
        _print_design(c,loc)
    return loc


#%%
def _print_design(c,loc):
    """
    Print design parameters and optimal design points from given data.
    
    This function outputs key design parameters and a list of design points to the
    console. It formats and displays minimum and maximum cell voltage, maximum 
    discharge current, maximum discharge power, and coordinates of the design 
    points. The optimal design point is highlighted at the end.
    
    Parameters:
    - c (dict): Dictionary containing cell design parameters.
    - loc (dict): Dictionary containing location data of design points.
    
    Returns:
    None - This function prints output to the console.
    
    Notes:
    - 'c' should contain keys 'Udc_min', 'Udc_max', 'Idc_DIS_max', 'Pdc_DIS_max'.
    - 'loc' should contain keys 'xy_list' and 'xy_opt', which hold design point
      coordinates and optimal design point data, respectively.
    """
    print('\n---')
    print('Udc_min = {:.3f} V'.format(c['cell']['Udc_min']))
    print('Udc_max = {:.3f} V'.format(c['cell']['Udc_max']))
    print('Idc_DIS_max = {:.3f} A'.format(c['cell']['Idc_DIS_max']))
    print('Pdc_DIS_max = {:.3f} W'.format(c['cell']['Pdc_DIS_max']))
    # print('xS = {}'.format(loc['xx']))
    # print('yP = {}'.format(loc['yy']))
    [print('xSyP = ({}, {})'.format(i[0],i[1])) for i in loc['xy_list']]
    print('---')
    print('Optimal Design Point:')
    print('xSyP_opt = {}'.format(tuple(loc['xy_opt'][:-1])))
    print('n_opt = {}'.format(loc['xy_opt'][-1]))
    
    
        
#%%
def _print_limit_values(g, c, zfig, zax, print_pec=False, xytext=(-95,0),
                        color='k', ha='left', va='top', 
                        bbox_args=dict(boxstyle="round", color='1', fc='1', 
                                       ec='None', lw=0)):
    """
    Annotate the plot with cell and power converter limit values.

    This function adds annotations to a plot, detailing the limit values of cells
    and power converters in an energy storage system. It lists maximum power, 
    voltage, and current limits for both cells and power converters. The function
    allows customization of annotation placement, color, and box styling.

    Parameters:
    - g (dict): Dictionary with matplotlib graph objects.
    - c (dict): Dictionary containing cell and power converter limit values.
    - zfig (int): Figure index for annotation.
    - zax (int): Axis index for annotation.
    - printinv (bool, optional): Flag to include power converter limits. 
      Defaults to False.
    - xytext (tuple, optional): Offset (in points) from the default annotation
      position. Defaults to (-95,0).
    - color (str, optional): Color of the annotation text. Defaults to 'k' (black).
    - ha (str, optional): Horizontal alignment of the annotation. Defaults to 'left'.
    - va (str, optional): Vertical alignment of the annotation. Defaults to 'top'.
    - bbox_args (dict, optional): Box style arguments for the annotation. Defaults
      to a rounded, white, and borderless box.

    Returns:
    None - This function adds annotations to an existing plot.

    Notes:
    - The annotations are formatted with LaTeX for readability.
    - Placement and styling of annotations are customizable for flexibility.
    """
    
    an_text = r'\textbf{Cell Limits}'
    an_text += '\n\n'
    an_text += r'$P_\mathrm{{cell,max}}^\mathrm{{dis}}$ = {:.2f} W'.format(c['cell']['Pdc_DIS_max'])
    an_text += '\n'
    an_text += r'$U_\mathrm{{cell,min}}$ = {:.2f} V'.format(c['cell']['Udc_min'])
    an_text += '\n'
    an_text += r'$U_\mathrm{{cell,max}}$ = {:.2f} V'.format(c['cell']['Udc_max'])
    an_text += '\n'
    an_text += r'$I_\mathrm{{cell,max}}^\mathrm{{dis}}$ = {:.2f} A'.format(c['cell']['Idc_DIS_max'])
    
    if print_pec:
        an_text += '\n\n\n'
        an_text += r'\textbf{Converter Limits}'
        an_text += '\n\n'
        an_text += r'$P_\mathrm{{pec,max}}^\mathrm{{DC}}$ = {:.0f} W'.format(c['pec']['Pdc_max'])
        an_text += '\n'
        an_text += r'$U_\mathrm{{pec,min}}^\mathrm{{DC}}$ = {:.0f} V'.format(c['pec']['Udc_min'])
        an_text += '\n'
        an_text += r'$U_\mathrm{{pec,max}}^\mathrm{{DC}}$ = {:.0f} V'.format(c['pec']['Udc_max'])
        an_text += '\n'
        an_text += r'$I_\mathrm{{pec,max}}^\mathrm{{DC}}$ = {:.0f} A'.format(c['pec']['Idc_max'])
    
    g['ax_{}_{}'.format(zfig,zax)].annotate(an_text,
                                        (1,0.875),
                                        xycoords='axes fraction',
                                        xytext=xytext,
                                        textcoords='offset points',
                                        size='smaller',
                                        c=color,
                                        ha=ha,
                                        va=va,
                                        bbox=bbox_args,
                                        zorder=100)
    
    
#%%
def plot_EP_limit(g, zfig, zax, res, n, nSoC=0, prefix_str='cell_model[1,1].', 
                  calc_EPsys=False, print_ann=False, plot_limit=[]):
    """
    Plot extended Ragone plot (ERP) curves, representing operational limits of 
    energy storage systems.

    This function generates curves on a Ragone plot (energy-per-power), extended 
    to include various operational limits such as minimum state of charge (SoC),
    voltage limits, current limits (Crate), and temperature limits. It visualizes
    how these parameters impact the energy storage system's performance, offering
    a detailed view of feasible operational regions under different constraints.

    Parameters:
    - g (dict): Dictionary with matplotlib graph objects.
    - zfig (int): Figure index for plotting.
    - zax (int): Axis index for plotting.
    - res (list): List of dictionaries containing simulation results.
    - n (int): Index specifying which result set from 'res' to use.
    - nSoC (int, optional): State of charge index for plotting. Defaults to 0.
    - prefix_str (str, optional): Sub-model string for result keys. Defaults to 
      'cell_model[1,1]'.
    - calc_EPsys (bool, optional): Flag to calculate system EP. Defaults to False.
    - print_ann (bool, optional): Flag to print annotations. Defaults to False.
    - plot_limit (list, optional): List of limits to plot. Defaults to empty list.

    Returns:
    None - Modifies the provided graph object with new plots.

    Notes:
    - The function enhances the traditional Ragone plot by including cell-based
      quantities, creating an ERP.
    - Dynamically recalculates and plots EP limits based on cell parameters.
    - Annotations are provided for clarity, highlighting specific limit values.
    """
    
    # Settings
    farbe = plot_preferences.plot_pre()[0]
    bbox_args = dict(boxstyle="round",color='None',fc='1',lw=0)

    n0 = 0

    # Limits
    limit = {'SoCmin':      np.arange(0,1,0.1).tolist(),
             'Umin':        np.arange(res[n]['U_min'][0]+0.3,res[n]['U_max'][0]-0.5+0.1,0.05).tolist(),
             'Umax':        np.arange(res[n]['U_min'][0],res[n]['U_max'][0]+0.1,0.05).T.tolist(),
             'Cratemax':    np.flip(np.arange(0.5,res[n]['Crate_max'][0]+0.5,0.5)).tolist(),
             'Tmax':        np.arange(25,res[n]['T_max'][0],0.5).tolist(),
             }
    for key in limit:
        if key not in plot_limit:
            limit[key] = []

    # SoCmin
    nSoCmin = 0
    for nSoCmin,u in enumerate(limit['SoCmin']):
        res_recalc = recalc_limits(res, prefix_str=prefix_str, Umax=res[n]['U_max'][0], Umin=res[n]['U_min'][0], Cratemax=res[n]['Crate_max'][0], Tmax=res[n]['T_max'][0], SoCmin=u)
        d_recalc, tmp_ind_recalc, rn_min_recalc = compile_EP(res_recalc, prefix_str=prefix_str, calc_EPsys=calc_EPsys)
  
        g['pl_{}_{}_{}_{}'.format(zfig,zax,4,n0+nSoCmin)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][0],color=farbe[8],marker='.',markersize=4,lw=1)
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,5,n0+nSoCmin)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][1],color=farbe[8],marker='v',markersize=4,ls='None')
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,6,n0+nSoCmin)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][2],color=farbe[8],marker='o',markersize=4,ls='None')
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,7,n0+nSoCmin)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][3],color=farbe[8],marker='d',markersize=4,ls='None')
      
        # Annotation
        if print_ann:
            n = 3+3*nSoCmin
            g['ax_{}_{}'.format(zfig,zax)].annotate(r'{:.0f}\,\%'.format(round(u*100,1)),
                                                    (d_recalc['p'][n],d_recalc[nSoC][0][n]),
                                                    textcoords="offset points",
                                                    xytext=(0,0) if nSoCmin>0 else (10.5,8.5),
                                                    c=farbe[8],
                                                    size='smaller',
                                                    ha='center',
                                                    va='center',
                                                    bbox=bbox_args,
                                                    zorder=100)
    n0 += nSoCmin
  
    # Umin
    nUmin = 0
    for nUmin,u in enumerate(limit['Umin']):
        res_recalc = recalc_limits(res, prefix_str=prefix_str, Umax=res[n]['U_max'][0], Umin=u, Cratemax=res[n]['Crate_max'][0], Tmax=res[n]['T_max'][0])
        d_recalc, tmp_ind_recalc, rn_min_recalc = compile_EP(res_recalc,prefix_str=prefix_str, calc_EPsys=calc_EPsys)
  
        g['pl_{}_{}_{}_{}'.format(zfig,zax,4,n0+nUmin)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][0],color=farbe[5],marker='.',markersize=4,lw=1)
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,5,n0+nUmin)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][1],color=farbe[5],marker='v',markersize=4,ls='None')
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,6,n0+nUmin)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][2],color=farbe[5],marker='o',markersize=4,ls='None')
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,7,n0+nUmin)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][3],color=farbe[5],marker='d',markersize=4,ls='None')
             
        # Annotation
        if print_ann:
            n = (3 if nUmin == 13 else
                  9 if nUmin == 12 else # 2.4 V
                  15 if nUmin == 11 else
                  20 if nUmin == 10 else # 2.3 V
                  22 if nUmin == 9 else
                  23 if nUmin == 8 else # 2.2 V
                  21 if nUmin == 7 else
                  18 if nUmin == 6 else # 2.1 V
                  17 if nUmin == 5 else
                  18 if nUmin == 4 else # 2.0 V
                  20 if nUmin == 3 else
                  23 if nUmin == 2 else # 1.9 V
                  27 if nUmin == 1 else
                  30) # 1.8 V
            g['ax_{}_{}'.format(zfig,zax)].annotate(r'{}\,V'.format(round(u,2)) if nUmin>2 else r'${}\,\dots\,{}\,V$'.format(round(u,2),round(res[n]['U_min'][0],2)) if nUmin==2 else '',
                                                    (d_recalc['p'][n],d_recalc[nSoC][0][n]),
                                                    textcoords="offset points",
                                                    xytext = ((5,5.5) if nUmin==13 else
                                                              (3,4.5) if nUmin==12 else
                                                              (1,1) if nUmin==11 else
                                                              # (0,-2) if nUmin==5 else
                                                              # (0,-3) if nUmin==4 else
                                                              (1,2) if nUmin==3 else
                                                              (27.5,14) if nUmin==2 else
                                                              (0,0)),
                                                    c=farbe[5],
                                                    size='smaller',
                                                    ha='center',
                                                    va='center',
                                                    bbox=bbox_args,
                                                    zorder=100)
    n0 += nUmin
 
    # Umax
    nUmax = 0
    for nUmax,u in enumerate(limit['Umax']):
        res_recalc = recalc_limits(res, prefix_str=prefix_str, Umax=u, Umin=res[n]['U_min'][0], Cratemax=res[n]['Crate_max'][0], Tmax=res[n]['T_max'][0])
        d_recalc, tmp_ind_recalc, rn_min_recalc = compile_EP(res_recalc, prefix_str=prefix_str, calc_EPsys=calc_EPsys)
  
        g['pl_{}_{}_{}_{}'.format(zfig,zax,4,n0+nUmax)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][0],color=farbe[10],marker='.',markersize=4,lw=1)
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,5,n0+nUmax)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][1],color=farbe[10],marker='v',markersize=4,ls='None')
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,6,n0+nUmax)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][2],color=farbe[10],marker='o',markersize=4,ls='None')
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,7,n0+nUmax)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][3],color=farbe[10],marker='d',markersize=4,ls='None')
    n0 += nUmax
  
    # Cratemax
    nCratemax = 0
    for nCratemax,u in enumerate(limit['Cratemax']):
        res_recalc = recalc_limits(res, prefix_str=prefix_str, Umax=res[n]['U_max'][0], Umin=res[n]['U_min'][0], Cratemax=u, Tmax=res[n]['T_max'][0])
        d_recalc, tmp_ind_recalc, rn_min_recalc = compile_EP(res_recalc, prefix_str=prefix_str, calc_EPsys=calc_EPsys)
  
        g['pl_{}_{}_{}_{}'.format(zfig,zax,4,n0+nCratemax+1)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][0],color=farbe[3],marker='.',markersize=4,alpha=1,lw=1)
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,5,n0+nCratemax+1)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][1],color=farbe[3],marker='v',markersize=4,ls='None')
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,6,n0+nCratemax+1)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][2],color=farbe[3],marker='o',markersize=4,ls='None')
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,7,n0+nCratemax+1)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][3],color=farbe[3],marker='d',markersize=4,ls='None')
  
        # Annotation
        if print_ann:
            f_int=interpolate.interp1d(d_recalc[nSoC][0],d_recalc['p'])
            g['ax_{}_{}'.format(zfig,zax)].annotate(r'{}\,C'.format(round(u,1)) if nCratemax>0 else '',
                                                    (f_int(3+6*nCratemax),3+6*nCratemax),
                                                    textcoords="offset points",
                                                    xytext=((6,0) if nCratemax==11 else #(0,0) without unit
                                                            (0,0)),
                                                    c=farbe[3],
                                                    size='smaller',
                                                    ha='center',
                                                    va='center',
                                                    bbox=bbox_args,
                                                    zorder=100)
    n0 += nCratemax
 
    # Tmax
    nTmax = 0
    for nTmax,u in enumerate(limit['Tmax']):
        if nTmax > 12:
            continue
        res_recalc = recalc_limits(res, prefix_str=prefix_str, Umax=res[n]['U_max'][0], Umin=res[n]['U_min'][0], Cratemax=res[n]['Crate_max'][0], Tmax=u)
        d_recalc, tmp_ind_recalc, rn_min_recalc = compile_EP(res_recalc, prefix_str=prefix_str, calc_EPsys=calc_EPsys)
  
        g['pl_{}_{}_{}_{}'.format(zfig,zax,4,n0+nTmax+1)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][0],color=farbe[2],marker='.',markersize=4,alpha=1,lw=1)
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,5,n0+nTmax+1)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][1],color=farbe[2],marker='v',markersize=4,ls='None')
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,6,n0+nTmax+1)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][2],color=farbe[2],marker='o',markersize=4,ls='None')
        # g['pl_{}_{}_{}_{}'.format(zfig,zax,7,n0+nTmax+1)] = g['ax_{}_{}'.format(zfig,zax)].plot(d_recalc['p'],d_recalc[nSoC][3],color=farbe[2],marker='d',markersize=4,ls='None')   
  
        # Annotation
        if print_ann:
            n = (1 if nTmax == 0 else # 25 °C
                  1 if nTmax == 1 else 
                  8 if nTmax == 2 else # 26 °C
                  15 if nTmax == 3 else 
                  19 if nTmax == 4 else # 27 °C
                  21 if nTmax == 5 else 
                  24 if nTmax == 6 else # 28 °C
                  26 if nTmax == 7 else
                  27 if nTmax == 8 else # 29 °C
                  28 if nTmax == 9 else
                  29 if nTmax == 10 else # 30 °C
                  29 if nTmax == 11 else 
                  31)
            g['ax_{}_{}'.format(zfig,zax)].annotate(r'{}\,$^\circ$C'.format(round(u,1)) if nTmax<12 else r'${}\,\dots\,{}\,^\circ$C'.format(round(u,2),round(res[n]['T_max'][0],2)) if nTmax>0 else '',
                                                    (d_recalc['p'][n],d_recalc[nSoC][0][n]),
                                                    textcoords="offset points",
                                                    xytext=((12.5,14) if nTmax==12 else
                                                            (-1,-4.75) if nTmax==11 else
                                                            (-3,-3) if nTmax==10 else
                                                            (-2,-4.5) if nTmax==9 else
                                                            (-2,-0.5) if nTmax==8 else
                                                            (-3,3.5) if nTmax==7 else
                                                            (-1,-2.5) if nTmax==6 else
                                                            (-1,0) if nTmax==4 else
                                                            (-1,2.5) if nTmax==2 else
                                                            (0,0)),
                                                    c=farbe[2],
                                                    size=0.8*plot_preferences.plt.rcParams['font.size'],#'smaller',
                                                    ha='center',
                                                    va='center',
                                                    bbox=bbox_args,
                                                    zorder=100)
    n0 += nTmax

    # Envelope
    d, tmp_ind, rn_min = compile_EP(res, prefix_str=prefix_str, calc_EPsys=calc_EPsys)
    
    g['pl_{}_{}_{}_{}'.format(zfig,zax,0,n0+nSoC+3)] = g['ax_{}_{}'.format(zfig,zax)].plot(d['p'],d[nSoC][0],color=farbe[nSoC],marker='.',markersize=6,alpha=1,lw=1.5,zorder=100)#,label='{}'.format(round(SoC_set[nSoC],1)),zorder=100)
    # g['pl_{}_{}_{}_{}'.format(zfig,zax,1,n0+nSoC+2)] = g['ax_{}_{}'.format(zfig,zax)].plot(d['p'],d[nSoC][1],color=farbe[nSoC],marker='v',markersize=6,ls='None',alpha=alph)
    # g['pl_{}_{}_{}_{}'.format(zfig,zax,2,n0+nSoC+2)] = g['ax_{}_{}'.format(zfig,zax)].plot(d['p'],d[nSoC][2],color=farbe[nSoC],marker='o',markersize=6,ls='None',alpha=alph)
    # g['pl_{}_{}_{}_{}'.format(zfig,zax,3,n0+nSoC+2)] = g['ax_{}_{}'.format(zfig,zax)].plot(d['p'],d[nSoC][3],color=farbe[nSoC],marker='d',markersize=6,ls='None',alpha=alph)
    
    


#%%
def _get_fcn_angle_loglog(ax, f):
    """
    Calculates the angle of a function on a logarithmic plot.

    Given a matplotlib axis with logarithmic scale and a function 'f', this function
    computes the angle of 'f' relative to the axis, which is useful for annotation 
    positioning.

    Parameters:
    - ax (matplotlib.axes.Axes): The matplotlib axis with a logarithmic scale.
    - f (function): The function for which the angle is to be calculated.

    Returns:
    - float: The angle of the function 'f' in degrees.

    Note:
    - Intended for use with logarithmic scales where both axes are in log scale.
    - Utilizes NumPy for logarithmic and trigonometric calculations.
    """
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    dx = np.log(xlim[1]) - np.log(xlim[0])
    dy = np.log(f(xlim[1])) - np.log(f(xlim[0]))
    
    return np.rad2deg(np.arctan2(dy/np.diff(np.log(ylim)), dx/np.diff(np.log(xlim))))[0]


    
#%%
def plot_limit_loci(g, zfig, zax, c, res, prefix_str='cell_model[1,1].',
                    dp_alpha=1, plot_new_dp=False, scale_log=True,
                    printapp=False, zoom_text='DSS', **kwargs):
    """
    Graphically synthesize constraint intersections of a combined CSP and ERP 
    for the design of energy storage systems.

    This function visualizes the intersection points of constraints within the 
    context of a Constraint Satisfaction Problem (CSP) combined with an extended 
    Ragone plot (ERP). It traces the Design Solution Space (DSS) by plotting 
    trajectories of constraint intersection points based on cell-based limit 
    quantities. The function allows dynamic adjustment of plot elements, enabling
    a comprehensive analysis of how different parameters influence the system's 
    operational feasibility and performance.

    Parameters:
    - g (dict): Dictionary with matplotlib graph objects.
    - zfig (int): Figure index for plotting.
    - zax (int): Axis index for plotting.
    - c (dict): Dictionary containing system design parameters.
    - res (list): List of dictionaries containing simulation results.
    - prefix_str (str, optional): Sub-model string for result keys. Defaults to 
      'cell_model[1,1]'.
    - dp_alpha (float, optional): Alpha value for data points. Defaults to 1.
    - plot_new_dp (bool, optional): Flag to plot new data points. Defaults to False.
    - scale_log (bool, optional): Flag to use logarithmic scaling. Defaults to True.
    - printapp (bool, optional): Flag to print application limits. Defaults to False.
    - zoom_text (str, optional): Text for zoom annotation. Defaults to 'DSS'.
    - **kwargs: Additional keyword arguments.

    Returns:
    function: Hyperbolic function representing the EP limit curve.

    Notes:
    - Essential in the graphical synthesis approach, facilitating the visualization
      of the DSS within the combined framework of CSP and ERP.
    - Adapts the underlying plot to reflect changes in ERP-based limit quantities
      and their impact on the system's DSS.
    """
    
    # Settings
    farbe = plot_preferences.plot_pre()[0]
    lstyle  = plot_preferences.plot_pre()[-1]
    bbox_args = dict(boxstyle="round",color='None',fc='1',lw=0)
    
    
    # Operation
    set_EP = c['req']['E/P_req']
    set_P = c['cell']['Pdc_DIS_max']
    
    dp_intersection = calc_intersection_curve(g,res,set_EP,prefix_str=prefix_str,print_error=False)
    EPline_x = dp_intersection['valx_P']
    EPline_y = dp_intersection[prefix_str+'P_cell']
    id_min = min([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
    id_max = max([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
    sPmaxset = EPline_x[id_min] if set_P < EPline_x[id_min] else (EPline_x[id_max] if set_P > EPline_x[id_max] else set_P)
   
    
    # Evolution of Limits Intersection Points (Loci)
    xmin = c['pec']['Udc_min']/dp_intersection[prefix_str+'U_cell']
    xmax = c['pec']['Udc_max']/dp_intersection['U_max']
    ymin = c['pec']['n']*c['pec']['Idc_max']/dp_intersection[prefix_str+'I_cell']
    f_hypP = hyperbel((c['req']['Pac_req']/c['pec']['eta_op'])/dp_intersection[prefix_str+'P_cell'])
    
    if 'pl_{}_{}_{}'.format(zfig,zax,0) not in g:
        g['pl_{}_{}_{}'.format(zfig,zax,0)] = g['ax_{}_{}'.format(zfig,zax)].plot(xmin,f_hypP(xmin),color=farbe[5],ls=lstyle[1],zorder=100)
        g['pl_{}_{}_{}'.format(zfig,zax,1)] = g['ax_{}_{}'.format(zfig,zax)].plot(f_hypP(ymin),ymin,color=farbe[3],ls=lstyle[1],zorder=100)
        g['pl_{}_{}_{}'.format(zfig,zax,2)] = g['ax_{}_{}'.format(zfig,zax)].plot(xmax,f_hypP(xmax),color=farbe[4],ls=lstyle[1],zorder=100)
    else:
        g['pl_{}_{}_{}'.format(zfig,zax,0)][0].set_data(xmin,f_hypP(xmin))
        g['pl_{}_{}_{}'.format(zfig,zax,1)][0].set_data(f_hypP(ymin),ymin)
        g['pl_{}_{}_{}'.format(zfig,zax,2)][0].set_data(xmax,f_hypP(xmax))

       
    # Evolution of Limits Intersection Points (Loci)
    _,sUminset = _set_intersection_limits(g,EPline_x,dp_intersection[prefix_str+'U_cell'],zfig,zax,x=sPmaxset,replot=True)
    _,sImaxset = _set_intersection_limits(g,EPline_x,dp_intersection[prefix_str+'I_cell'],zfig,zax,x=sPmaxset,replot=True)
    _,sUmaxset = _set_intersection_limits(g,EPline_x,dp_intersection['U_max'],zfig,zax,x=sPmaxset,replot=True)
   
    xxmin = c['pec']['Udc_min']/sUminset
    xxmax = c['pec']['Udc_max']/sUmaxset
    yymin = c['pec']['n']*c['pec']['Idc_max']/sImaxset
    ff_hypP = hyperbel((c['req']['Pac_req']/c['pec']['eta_op'])/sPmaxset)
    
        
    # Limit
    xlim = kwargs['xlim'] if ('xlim' in kwargs and kwargs['xlim']) else g['ax_{}_{}'.format(zfig,zax)].get_xlim()
    g['ax_{}_{}'.format(zfig,zax)].set_xlim(xlim)

    ylim = kwargs['ylim'] if ('ylim' in kwargs and kwargs['ylim']) else g['ax_{}_{}'.format(zfig,zax)].get_ylim()
    g['ax_{}_{}'.format(zfig,zax)].set_ylim(ylim)
    
    if scale_log:
        g['ax_{}_{}'.format(zfig,zax)].set_xscale('log')
        g['ax_{}_{}'.format(zfig,zax)].set_yscale('log')
    
    
    xvec = np.arange(xlim[0] if xlim[0]>0 else 1.0, xlim[1], 1.0)
    if 'pl_{}_{}_{}'.format(zfig,zax,120) not in g or plot_new_dp:
        g['pl_{}_{}_{}'.format(zfig,zax,120)] = g['ax_{}_{}'.format(zfig,zax)].plot(xvec,ff_hypP(xvec),color=farbe[0],alpha=dp_alpha,lw=1,zorder=100)           
        g['pl_{}_{}_{}'.format(zfig,zax,121)] = g['ax_{}_{}'.format(zfig,zax)].plot(xxmin,ff_hypP(xxmin),marker='>',ms=6.5,color=farbe[5],mfc='1',lw=1.25,alpha=dp_alpha,zorder=100)
        g['pl_{}_{}_{}'.format(zfig,zax,122)] = g['ax_{}_{}'.format(zfig,zax)].plot(ff_hypP(yymin),yymin,marker='>',ms=6.5,color=farbe[3],mfc='1',lw=1.25,alpha=dp_alpha,zorder=100)
        g['pl_{}_{}_{}'.format(zfig,zax,123)] = g['ax_{}_{}'.format(zfig,zax)].plot(xxmax,ff_hypP(xxmax),marker='<',ms=6.5,color=farbe[4],mfc='1',lw=1.25,alpha=dp_alpha,zorder=100)
        
        # Zoom Box
        xlim_zoom = g['ax_{}_{}'.format(3,0)].get_xlim()
        ylim_zoom = g['ax_{}_{}'.format(3,0)].get_ylim()
        g['pl_{}_{}_{}'.format(zfig,zax,124)] = g['ax_{}_{}'.format(zfig,zax)].add_patch(Rectangle((xlim_zoom[0],ylim_zoom[0]),xlim_zoom[-1]-xlim_zoom[0],ylim_zoom[-1]-ylim_zoom[0],ec='k',fc='None',ls='dashed',lw=1,alpha=dp_alpha,zorder=100))
        g['an_{}_{}_{}'.format(zfig,zax,124)] = g['ax_{}_{}'.format(zfig,zax)].annotate(zoom_text,
                                                                            (xlim_zoom[0],ylim_zoom[0]),
                                                                            textcoords="offset points",
                                                                            xytext=(0,-5),
                                                                            size='smaller',
                                                                            ha='center',
                                                                            va='top',
                                                                            alpha=dp_alpha,
                                                                            bbox=bbox_args,
                                                                            zorder=100)
        
    else:
        g['pl_{}_{}_{}'.format(zfig,zax,120)][0].set_data(xvec,ff_hypP(xvec))
        g['pl_{}_{}_{}'.format(zfig,zax,121)][0].set_data(xxmin,ff_hypP(xxmin))
        g['pl_{}_{}_{}'.format(zfig,zax,122)][0].set_data(ff_hypP(yymin),yymin)
        g['pl_{}_{}_{}'.format(zfig,zax,123)][0].set_data(xxmax,ff_hypP(xxmax))
        
        # Zoom Box
        xlim_zoom = g['ax_{}_{}'.format(3,0)].get_xlim()
        ylim_zoom = g['ax_{}_{}'.format(3,0)].get_ylim()
        g['pl_{}_{}_{}'.format(zfig,zax,124)].xy = (xlim_zoom[0],ylim_zoom[0])
        g['pl_{}_{}_{}'.format(zfig,zax,124)].set_width(xlim_zoom[-1]-xlim_zoom[0])
        g['pl_{}_{}_{}'.format(zfig,zax,124)].set_height(ylim_zoom[-1]-ylim_zoom[0])
        g['an_{}_{}_{}'.format(zfig,zax,124)].set_text('Zoom '+zoom_text)
        g['an_{}_{}_{}'.format(zfig,zax,124)].xy = (xlim_zoom[0],ylim_zoom[0])
    
    # Annotations   
    bbox_args = dict(boxstyle="round",color='None',fc='1',lw=0)
    if 'an_{}_{}_{}'.format(zfig,zax,10) not in g or plot_new_dp:
        
        angle = _get_fcn_angle_loglog(g['ax_{}_{}'.format(zfig,zax)],ff_hypP)
        ycoord = np.exp(np.log(ylim[0])+np.diff(np.log(ylim))*0.96)
        f_int = interpolate.interp1d(f_hypP(xmax),xmax) # interpolation with 'np.interp' caused problems!
        g['an_{}_{}_{}'.format(zfig,zax,10)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'$\mathcal{I}(\textbf{\textsf{C}}_{U_\mathrm{max}}\cap\textbf{\textsf{C}}_{P_\mathrm{req}}$)',
                                                      (f_int(ycoord),ycoord),
                                                      textcoords="offset points",
                                                      xytext=(0,0),
                                                      c=farbe[4],
                                                      # fontsize=11,
                                                      size='medium',
                                                      ha='center',
                                                      va='center',
                                                      # transform_rotates_text=True,
                                                      # rotation=angle,
                                                      # rotation_mode='anchor',
                                                      bbox=bbox_args,
                                                      zorder=100)
        
        ycoord = np.exp(np.log(ylim[0])+np.diff(np.log(ylim))*0.96)
        f_int = interpolate.interp1d(f_hypP(xmin),xmin)
        g['an_{}_{}_{}'.format(zfig,zax,11)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'$\mathcal{I}(\textbf{\textsf{C}}_{U_\mathrm{min}}\cap\textbf{\textsf{C}}_{P_\mathrm{req}}$)',
                                                      (f_int(ycoord),ycoord),
                                                      textcoords="offset points",
                                                      xytext=(0,0),
                                                      c=farbe[5],
                                                      # fontsize=11,
                                                      size='medium',
                                                      ha='center',
                                                      va='center',
                                                      # transform_rotates_text=True,
                                                      # rotation=angle,
                                                      # rotation_mode='anchor',
                                                      bbox=bbox_args,
                                                      zorder=100)
        
        ycoord = np.exp(np.log(ylim[0])+np.diff(np.log(ylim))*0.96)
        f_int = interpolate.interp1d(ymin,f_hypP(ymin))
        g['an_{}_{}_{}'.format(zfig,zax,12)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'$\mathcal{I}(\textbf{\textsf{C}}_{I_\mathrm{max}}\cap\textbf{\textsf{C}}_{P_\mathrm{req}}$)',
                                                      (f_int(ycoord),ycoord),
                                                      textcoords="offset points",
                                                      xytext=(0,0),
                                                      c=farbe[3],
                                                      # fontsize=11,
                                                      size='medium',
                                                      ha='center',
                                                      va='center',
                                                      # transform_rotates_text=True,
                                                      # rotation=angle,
                                                      # rotation_mode='anchor',
                                                      bbox=bbox_args,
                                                      zorder=100)
        
        
        xcoord = np.exp(np.log(xlim[0])+np.diff(np.log(xlim))*0.03)
        g['an_{}_{}_{}'.format(zfig,zax,13)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'$\textbf{\textsf{C}}_{P_\mathrm{req}}$',
                                                (xcoord,ff_hypP(xcoord)),
                                                textcoords="offset points",
                                                xytext=(0,0),
                                                c=farbe[0],
                                                # fontsize=11,
                                                size='medium',
                                                ha='left',
                                                va='center',
                                                transform_rotates_text=True,
                                                rotation=angle,
                                                rotation_mode='anchor',
                                                bbox=bbox_args,
                                                zorder=100)
        
        
        xcoord = np.exp(np.log(xlim[0])+np.diff(np.log(xlim))*0.975)
        g['an_{}_{}_{}'.format(zfig,zax,14)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'$P^\mathrm{{dis}}_\mathrm{{cell,max}}$ = {:.2f} W'.format(c['cell']['Pdc_DIS_max']),
                                                    (xcoord,ff_hypP(xcoord)),
                                                    textcoords="offset points",
                                                    xytext=(0,8),
                                                    c=farbe[2],
                                                    size='smaller',
                                                    ha='right',
                                                    va='center',
                                                    transform_rotates_text=True,
                                                    rotation=angle,
                                                    rotation_mode='anchor',
                                                    bbox=bbox_args,
                                                    zorder=100)
        
    else:
        ycoord = np.exp(np.log(ylim[0])+np.diff(np.log(ylim))*0.96)
        f_int = interpolate.interp1d(f_hypP(xmax),xmax)
        g['an_{}_{}_{}'.format(zfig,zax,10)].xy = (f_int(ycoord),ycoord)
        
        ycoord = np.exp(np.log(ylim[0])+np.diff(np.log(ylim))*0.96)
        f_int = interpolate.interp1d(f_hypP(xmin),xmin)
        g['an_{}_{}_{}'.format(zfig,zax,11)].xy = (f_int(ycoord),ycoord)
        
        ycoord = np.exp(np.log(ylim[0])+np.diff(np.log(ylim))*0.96)
        f_int = interpolate.interp1d(ymin,f_hypP(ymin))
        g['an_{}_{}_{}'.format(zfig,zax,12)].xy = (f_int(ycoord),ycoord)
        
        xcoord = np.exp(np.log(xlim[0])+np.diff(np.log(xlim))*0.03)
        g['an_{}_{}_{}'.format(zfig,zax,13)].xy = (xcoord,ff_hypP(xcoord))
        
        xcoord = np.exp(np.log(xlim[0])+np.diff(np.log(xlim))*0.97)
        g['an_{}_{}_{}'.format(zfig,zax,14)].set_text(r'$P^\mathrm{{dis}}_\mathrm{{cell,max}}$ = {:.2f} W'.format(c['cell']['Pdc_DIS_max']))
        g['an_{}_{}_{}'.format(zfig,zax,14)].xy = (xcoord,ff_hypP(xcoord))

        
    if printapp:
        an_text = r'\textbf{Application Limits}'
        an_text += '\n\n'
        an_text += r'$(E/P)_\mathrm{{req}}$ = '+_text_EP_units(set_EP*60)
        if 'an_{}_{}_{}'.format(zfig,zax,20) not in g:
            g['an_{}_{}_{}'.format(zfig,zax,20)] = g['ax_{}_{}'.format(zfig,zax)].annotate(an_text,
                                                (1,0.875),
                                                xycoords='axes fraction',
                                                xytext=(-97.5,0),
                                                textcoords='offset points',
                                                size='smaller',
                                                c=farbe[0],
                                                ha='left',
                                                va='top',
                                                bbox=bbox_args,
                                                zorder=100)
        else:
            g['an_{}_{}_{}'.format(zfig,zax,20)].set_text(an_text)
    else:
        if 'an_{}_{}_{}'.format(zfig,zax,20) in g:
            g['an_{}_{}_{}'.format(zfig,zax,20)].remove()
            g.pop('an_{}_{}_{}'.format(zfig,zax,20))
    
    g['fig_{}'.format(zfig)].canvas.draw()
    
    return ff_hypP



#%%
def plot_iso_intersection_curves(g, zfig, zax, c, res, 
                                 prefix_str='cell_model[1,1].',
                                 var='general', 
                                 iso_mins=[10/60, 1, 2.5, 5, 7.5, 10, 15, 30, 60, 300],
                                 print_ann=True):
    """
    Calculate and plot iso-curves based on data from an extended Ragone plot (ERP)
    for various cell-based quantities at specific E/P ratios.

    This function generates iso-curves that illustrate intersections of different 
    cell-based quantities at predetermined energy-per-power (E/P) ratios. These 
    curves are derived from data obtained from the ERP and provide insights into 
    the operational behavior and limits of an energy storage system under various 
    E/P scenarios. The function supports analysis of variables such as minimum 
    voltage, maximum current, and general constraints, offering a detailed 
    perspective on the system's performance.

    Parameters:
    - g (dict): Dictionary with matplotlib graph objects.
    - zfig (int): Figure index for plotting.
    - zax (int): Axis index for plotting.
    - c (dict): Dictionary containing system design parameters.
    - res (list): List of dictionaries containing simulation results.
    - prefix_str (str, optional): Sub-model string for result keys. Defaults to 
      'cell_model[1,1]'.
    - var (str, optional): Variable for iso-curve plotting. Options include 'Umin', 
      'Imax', or 'general'. Defaults to 'general'.
    - iso_mins (list, optional): List of E/P ratios for plotting iso-curves. 
      Defaults to [10/60, 1, 2.5, 5, 7.5, 10, 15, 30, 60, 300].
    - print_ann (bool, optional): Flag to print annotations. Defaults to True.

    Returns:
    None - Modifies the provided graph object with new plots.

    Notes:
    - Enables the analysis of cell-based operational parameters in relation to 
      various E/P ratios, elucidating the system's constraints and capabilities.
    - Annotations can be used to highlight significant points on the iso-curves 
      and enhance understanding of operational scenarios.
    """
    
    # Settings
    farbe = plot_preferences.plot_pre()[0]
    lstyle  = plot_preferences.plot_pre()[-1]
    bbox_args = dict(boxstyle="round",color='None',fc='1',lw=0)
    alph = 0.33

    # Operation
    for nEP,set_EP in enumerate(iso_mins):
        dp_intersection = calc_intersection_curve(g,res,set_EP,prefix_str=prefix_str,print_error=False)
        EPline_x = dp_intersection['valx_P']
        EPline_y = dp_intersection[prefix_str+'P_cell']
        id_min = min([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
        id_max = max([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
        
        if var=='Umin':
            g['pl_{}_{}_{}'.format(zfig,zax,nEP+200)] = g['ax_{}_{}'.format(zfig,zax)].plot(dp_intersection[prefix_str+'P_cell'],dp_intersection[prefix_str+'U_cell'],color=farbe[5],alpha=alph,lw=0.75,zorder=0)
        
            if print_ann:
                g['iso_{}_{}_{}'.format(zfig,zax,nEP)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'{}'.format(_text_EP_units(set_EP*60)),
                                                    (dp_intersection[prefix_str+'P_cell'][id_max],dp_intersection[prefix_str+'U_cell'][id_max]),
                                                    textcoords="offset points",
                                                    xytext=(0,-5),
                                                    c=farbe[5],
                                                    size='smaller',
                                                    ha='center',
                                                    va='top',
                                                    alpha=alph,
                                                    bbox=bbox_args,
                                                    zorder=0)
            
        elif var=='Imax':
            g['pl_{}_{}_{}'.format(zfig,zax,nEP+200)] = g['ax_{}_{}'.format(zfig,zax)].plot(dp_intersection[prefix_str+'P_cell'],dp_intersection[prefix_str+'I_cell'],color=farbe[3],alpha=alph,lw=0.75,zorder=0)

            if print_ann:
                g['iso_{}_{}_{}'.format(zfig,zax,nEP)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'{}'.format(_text_EP_units(set_EP*60)),
                                                    (dp_intersection[prefix_str+'P_cell'][id_max],dp_intersection[prefix_str+'I_cell'][id_max]),
                                                    textcoords="offset points",
                                                    xytext=(0,5),
                                                    c=farbe[3],
                                                    size='smaller',
                                                    ha='center',
                                                    va='top',
                                                    alpha=alph,
                                                    bbox=bbox_args,
                                                    zorder=0)
            
        elif var=='general':
            # Evolution of Limits Intersection Points (Loci)
            xmin = c['pec']['Udc_min']/dp_intersection[prefix_str+'U_cell']
            xmax = c['pec']['Udc_max']/dp_intersection['U_max']
            ymin = c['pec']['n']*c['pec']['Idc_max']/dp_intersection[prefix_str+'I_cell']
            f_hypP = hyperbel((c['req']['Pac_req']/c['pec']['eta_op'])/dp_intersection[prefix_str+'P_cell'])
            
            if 'pl_{}_{}_{}'.format(zfig,zax,nEP+200) not in g:                
                g['pl_{}_{}_{}'.format(zfig,zax,nEP+200)] = g['ax_{}_{}'.format(zfig,zax)].plot(xmin,f_hypP(xmin),color=farbe[5],alpha=alph,lw=0.75,zorder=0)
                g['pl_{}_{}_{}'.format(zfig,zax,nEP+len(iso_mins)+200)] = g['ax_{}_{}'.format(zfig,zax)].plot(f_hypP(ymin),ymin,color=farbe[3],alpha=alph,lw=0.75,zorder=0)
                g['pl_{}_{}_{}'.format(zfig,zax,nEP+2*len(iso_mins)+200)] = g['ax_{}_{}'.format(zfig,zax)].plot(xmax,f_hypP(xmax),color=farbe[4],alpha=alph,lw=0.75,zorder=0)
                
                if print_ann:
                    g['iso_{}_{}_{}'.format(zfig,zax,nEP)] = g['ax_{}_{}'.format(zfig,zax)].annotate(r'{}'.format(_text_EP_units(set_EP*60)),
                                                        (f_hypP(ymin)[~np.isnan(ymin)][-1],ymin[~np.isnan(ymin)][-1]),
                                                        textcoords="offset points",
                                                        xytext=(10,-2.5),
                                                        c=farbe[3],
                                                        size='smaller',
                                                        ha='center',
                                                        va='top',
                                                        alpha=alph,
                                                        bbox=bbox_args,
                                                        zorder=0)
            else:
                g['pl_{}_{}_{}'.format(zfig,zax,nEP+200)][0].set_data(xmin,f_hypP(xmin))
                g['pl_{}_{}_{}'.format(zfig,zax,nEP+len(iso_mins)+200)][0].set_data(f_hypP(ymin),ymin)
                g['pl_{}_{}_{}'.format(zfig,zax,nEP+2*len(iso_mins)+200)][0].set_data(xmax,f_hypP(xmax))
                
                if print_ann:
                    g['iso_{}_{}_{}'.format(zfig,zax,nEP)].xy = (f_hypP(ymin)[~np.isnan(ymin)][-1],ymin[~np.isnan(ymin)][-1])


#%%
def update_eta(g, c, zfig=5):
    """
    Update the efficiency (eta) plot in the graph based on the operational 
    requirements and maximum capabilities of the power converter.

    This function recalculates the efficiency (eta) of the power converter based 
    on the ratio of the required power output to the maximum power output 
    multiplied by the number of power converters (n). It then updates the 
    efficiency plot in the specified graph object to reflect the new efficiency 
    value. The function is designed to handle cases where the required power 
    exceeds the maximum capabilities by setting the ratio to 1.

    Parameters:
    - g (dict): Dictionary with matplotlib graph objects.
    - c (dict): Dictionary containing power converter and operational requirement 
      parameters.
    - zfig (int, optional): Figure index for updating the plot. Defaults to 5.

    Returns:
    float: The recalculated efficiency of the power converter at the given 
    operational point.

    Notes:
    - The function assumes the presence of specific plot and annotation objects 
      in the graph dictionary for updating.
    - Efficiency is calculated as a percentage and is updated on the graph 
      along with its corresponding annotation.
    """

    ratio_Pac = c['req']['Pac_req']/(c['pec']['Pac_max']*c['pec']['n'])
    ratio_Pac = 1 if ratio_Pac >= 1 else ratio_Pac 
    g['pl_{}_{}_{}'.format(zfig,0,10)][0].set_data([ratio_Pac]*2,[0,100])
    g['pl_{}_{}_{}'.format(zfig,0,11)][0].set_data([ratio_Pac],[c['pec']['fcn_eta(Umax)'](ratio_Pac)*100])
    g['an_{}_{}_{}'.format(zfig,0,11)].xy = (ratio_Pac,c['pec']['fcn_eta(Umax)'](ratio_Pac)*100)
    
    print('ratio_Pac = {:.3f}'.format(ratio_Pac))
    return c['pec']['fcn_eta(Umax)'](ratio_Pac)


#%%
def plot_constraint_intersection_trajectory(g, res, c, zfig=4,
                                            iso_mins=np.array([10/60, 1, 2.5, 5, 7.5, 10, 15, 30, 60, 300])/60,
                                            prefix_str='cell_model[1,1].'):
    """
    Adjust elements in the graphical synthesis for an energy storage system, 
    focusing on constraint intersection trajectories calculated from specific 
    energy-per-power (E/P) ratios.

    This function updates the graphical synthesis within the Constraint 
    Satisfaction Problem (CSP) framework using data from the extended Ragone 
    plot (ERP). It modifies the representation of constraint intersection 
    trajectories based on various E/P ratios, provided in minutes. The function 
    adjusts plots for voltage and current constraints and integrates the 
    visualization of constraint loci and iso-intersection curves, depicting the 
    system's operational requirements and cell characteristics.

    Parameters:
    - g (dict): Dictionary with matplotlib graph objects.
    - res (list): List of dictionaries containing simulation results.
    - c (dict): Dictionary containing system design parameters.
    - zfig (int, optional): Figure index for plotting. Defaults to 4.
    - iso_mins (np.array, optional): Array of E/P ratios in minutes for plotting 
      iso-intersection curves. Defaults to specific values ranging from 10/60 to 
      300/60.
    - prefix_str (str, optional): Sub-model string for result keys. Defaults to 
      'cell_model[1,1]'.

    Returns:
    None - Modifies the graphical synthesis with updated plots.

    Notes:
    - Plays a crucial role in the graphical synthesis process within the CSP 
      analysis, highlighting the intersections of various operational constraints.
    - Functions like `plot_limit_loci` and `plot_iso_intersection_curves` are 
      used for detailed visualization of operational constraints and behaviors 
      within the graphical synthesis.
    """

    # Operation
    set_EP = c['req']['E/P_req']
    set_P = c['cell']['Pdc_DIS_max']
    dp_intersection = calc_intersection_curve(g,res,set_EP,print_error=False,prefix_str=prefix_str)
    EPline_x = dp_intersection['valx_P']
    EPline_y = dp_intersection[prefix_str+'P_cell']
    id_min = min([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
    id_max = max([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
    sPmaxset = EPline_x[id_min] if set_P < EPline_x[id_min] else (EPline_x[id_max] if set_P > EPline_x[id_max] else set_P)


    # ax 0
    g['pl_{}_{}_{}'.format(zfig,0,0)][0].set_data(dp_intersection[prefix_str+'P_cell'],dp_intersection[prefix_str+'U_cell'])
    _,sUminset = _set_intersection_limits(g,EPline_x,dp_intersection[prefix_str+'U_cell'],zfig,0,x=sPmaxset,replot=True)
    
    # ax 1
    g['pl_{}_{}_{}'.format(zfig,1,0)][0].set_data(dp_intersection[prefix_str+'P_cell'],dp_intersection[prefix_str+'I_cell'])
    _,sImaxset = _set_intersection_limits(g,EPline_x,dp_intersection[prefix_str+'I_cell'],zfig,1,x=sPmaxset,replot=True)
    
    
    # ax 2
    plot_limit_loci(g,zfig,2,c,res,prefix_str=prefix_str,printapp=True)
    
    plot_iso_intersection_curves(g,zfig,2,c,res,prefix_str=prefix_str,iso_mins=iso_mins,var='general',print_ann=False)
        
    g['fig_{}'.format(zfig)].canvas.draw()
    
    
    
#%%
def replot_EP_limits(g, res, tmp_ind, Umax, Umin, Cratemax, Tmax, zfig=1, 
                     nPref=15, calc_EPsys=False, print_ann=False, 
                     prefix_str='cell_model[1,1].'):
    """
    Update and replot the Extended Ragone Plot (ERP) and measurement data
    based on new operational limits for an energy storage system.

    This function modifies the ERP and associated plots, reflecting revised
    limits such as maximum/minimum voltage, current rate (Crate), and 
    temperature. It adjusts the ERP and plots related to discharged cell 
    energy, ensuring accuracy in system performance visualization.

    Parameters:
    - g (dict): Dictionary with matplotlib graph objects.
    - res (list): Simulation results.
    - tmp_ind (list): Temporary indices for data selection.
    - Umax (float): Maximum voltage limit.
    - Umin (float): Minimum voltage limit.
    - Cratemax (float): Maximum current rate limit.
    - Tmax (float): Maximum temperature limit.
    - zfig (int, optional): Figure index for updating. Defaults to 1.
    - nPref (int, optional): Index for specific data points. Defaults to 15.
    - calc_EPsys (bool, optional): Flag for system EP calculation. Defaults to False.
    - print_ann (bool, optional): Flag for printing annotations. Defaults to False.
    - prefix_str (str, optional): Sub-model string for results. Defaults to 'cell_model[1,1]'.

    Returns:
    None - Modifies plots with updated ERP and data.

    Notes:
    - Crucial for maintaining ERP accuracy under new operational limits.
    - Aids in analyzing impacts of various limits on system behavior.
    """
    print("[Info] Submit New EP Limits")
    print("Umax =", Umax)
    print("Umin =", Umin)
    print("Cratemax =", Cratemax)
    print("Tmax =", Tmax)
   
    rn_min = 0
    Imax = Cratemax * res[rn_min]['_settings']['ParameterValues']['Q_n'][0]
   
    res_recalc = recalc_limits(res, Umax=Umax, Umin=Umin, Cratemax=Cratemax, Tmax=Tmax, prefix_str=prefix_str)
    d_recalc, tmp_ind_recalc, rn_min_recalc = compile_EP(res_recalc,calc_EPsys=calc_EPsys,prefix_str=prefix_str)
 
    for nSoC in range(max(tmp_ind)[2]+1):
        try:
            g['pl_{}_{}_{}_{}'.format(zfig,0,10,nSoC)][0].set_data(d_recalc['p'],d_recalc[nSoC][0])
            # g['pl_{}_{}_{}_{}'.format(zfig,0,11,nSoC)][0].set_data(d_recalc['p'],d_recalc[nSoC][1])
            # g['pl_{}_{}_{}_{}'.format(zfig,0,12,nSoC)][0].set_data(d_recalc['p'],d_recalc[nSoC][2])
            # g['pl_{}_{}_{}_{}'.format(zfig,0,13,nSoC)][0].set_data(d_recalc['p'],d_recalc[nSoC][3])
            # g['pl_{}_{}_{}_{}'.format(zfig,0,14,nSoC)][0].set_data(d_recalc['p'],d_recalc[nSoC][4])
            
            # Point
            g['pl_{}_{}_{}'.format(zfig,0,10)][0].set_data(d_recalc['p'][nPref],d_recalc[nSoC][0][nPref])
            
            # Umin
            res_recalc = recalc_limits(res, Umax=res[rn_min]['U_max'][0], Umin=Umin, Cratemax=res[rn_min]['Crate_max'][0], Tmax=res[rn_min]['T_max'][0],prefix_str=prefix_str)
            d_recalc, tmp_ind_recalc, rn_min_recalc = compile_EP(res_recalc,calc_EPsys=calc_EPsys,prefix_str=prefix_str)
            g['pl_{}_{}_{}_{}'.format(zfig,0,20,nSoC)][0].set_data(d_recalc['p'],d_recalc[nSoC][0])
            
            # Imax
            res_recalc = recalc_limits(res, Umax=res[rn_min]['U_max'][0], Umin=res[rn_min]['U_min'][0], Cratemax=Cratemax, Tmax=res[rn_min]['T_max'][0],prefix_str=prefix_str)
            d_recalc, tmp_ind_recalc, rn_min_recalc = compile_EP(res_recalc,calc_EPsys=calc_EPsys,prefix_str=prefix_str)
            g['pl_{}_{}_{}_{}'.format(zfig,0,21,nSoC)][0].set_data(d_recalc['p'],d_recalc[nSoC][0])
            
            # Umin & Imax
            res_recalc = recalc_limits(res, Umax=res[rn_min]['U_max'][0], Umin=Umin, Cratemax=Cratemax, Tmax=res[rn_min]['T_max'][0],prefix_str=prefix_str)
            d_recalc, tmp_ind_recalc, rn_min_recalc = compile_EP(res_recalc,calc_EPsys=calc_EPsys,prefix_str=prefix_str)

            rel_offset = 0.0025
            if any([Umin>res[rn_min]['U_min'][0],Cratemax<res[rn_min]['Crate_max'][0]]):
                g['pl_{}_{}_{}_{}'.format(zfig,0,10,nSoC)][0].set_visible(0)
                g['pl_{}_{}_{}_{}'.format(zfig,0,22,nSoC)][0].set_data(np.subtract(d_recalc['p'],rel_offset*max(d_recalc['p'])),
                                                                    np.subtract(d_recalc[nSoC][0],rel_offset*max(d_recalc[nSoC][0])))
                g['pl_{}_{}_{}_{}'.format(zfig,0,22,nSoC)][0].set_c(plot_preferences.plot_pre()[0][7])
                g['pl_{}_{}_{}_{}'.format(zfig,0,22,nSoC)][0].set_lw(3.5)
                g['pl_{}_{}_{}_{}'.format(zfig,0,22,nSoC)][0].set_zorder(0)
            else:
                g['pl_{}_{}_{}_{}'.format(zfig,0,10,nSoC)][0].set_visible(1)
                g['pl_{}_{}_{}_{}'.format(zfig,0,22,nSoC)][0].set_data(d_recalc['p'],d_recalc[nSoC][0])
                g['pl_{}_{}_{}_{}'.format(zfig,0,22,nSoC)][0].set_c(plot_preferences.plot_pre()[0][0])
                g['pl_{}_{}_{}_{}'.format(zfig,0,22,nSoC)][0].set_lw(1.5)
                g['pl_{}_{}_{}_{}'.format(zfig,0,22,nSoC)][0].set_zorder(g['pl_{}_{}_{}_{}'.format(zfig,0,10,nSoC)][0].get_zorder()+1)
            
            if print_ann:
                # P_cell_dis
                g['an_{}_{}_{}'.format(zfig,1,100)].set_visible(1)
                g['an_{}_{}_{}'.format(zfig,2,100)].set_visible(1)
                # P_cell_dis_UI
                g['pl_{}_{}_{}'.format(zfig,0,100)][0].set_visible(1)
                g['pl_{}_{}_{}'.format(zfig,0,101)][0].set_visible(1)
                g['an_{}_{}_{}'.format(zfig,0,100)].set_visible(1)
            else:
                # P_cell_dis
                g['an_{}_{}_{}'.format(zfig,1,100)].set_visible(0)
                g['an_{}_{}_{}'.format(zfig,2,100)].set_visible(0)
                # P_cell_dis_UI
                g['pl_{}_{}_{}'.format(zfig,0,100)][0].set_visible(0)
                g['pl_{}_{}_{}'.format(zfig,0,101)][0].set_visible(0)
                g['an_{}_{}_{}'.format(zfig,0,100)].set_visible(0)
        except:
            0
       
        for nP in range(max(tmp_ind)[1]+1):
            rn = [item[0] for item in tmp_ind if item[1]==nP and item[2]==nSoC][0]
            try: # Umin
                res_recalc = recalc_limits(res, Umax=res[rn_min]['U_max'][0], Umin=Umin, Cratemax=res[rn_min]['Crate_max'][0], Tmax=res[rn_min]['T_max'][0],prefix_str=prefix_str)
                d_recalc, tmp_ind_recalc, rn_min_recalc = compile_EP(res_recalc,calc_EPsys=calc_EPsys,prefix_str=prefix_str)
                g['pl_{}_{}_{}_{}'.format(zfig,1,1,rn)][0].set_data(res_recalc[rn][prefix_str+'E_cell'],res_recalc[rn][prefix_str+'U_cell'])
                g['pl_{}_{}_{}'.format(zfig,1,0)][0].set_ydata([Umin]*2)
                g['an_{}_{}_{}'.format(zfig,1,0)].set_y(Umin if Umin>res[rn_min]['U_min'][0] else g['an_{}_{}_{}'.format(zfig,1,0)].xy[1])
                g['an_{}_{}_{}'.format(zfig,1,0)].set_c(plot_preferences.plot_pre()[0][5] if Umin>res[rn_min]['U_min'][0] else plot_preferences.plot_pre()[0][0])
            except:
                0
            try: # Cratemax
                res_recalc = recalc_limits(res, Umax=res[rn_min]['U_max'][0], Umin=res[rn_min]['U_min'][0], Cratemax=Cratemax, Tmax=res[rn_min]['T_max'][0],prefix_str=prefix_str)
                d_recalc, tmp_ind_recalc, rn_min_recalc = compile_EP(res_recalc,calc_EPsys=calc_EPsys,prefix_str=prefix_str)
                g['pl_{}_{}_{}_{}'.format(zfig,2,1,rn)][0].set_data(res_recalc[rn][prefix_str+'E_cell'],res_recalc[rn][prefix_str+'I_cell'])
                g['pl_{}_{}_{}'.format(zfig,2,0)][0].set_ydata([Imax]*2)
                g['an_{}_{}_{}'.format(zfig,2,0)].set_y(Imax if Cratemax<res[rn_min]['Crate_max'][0] else g['an_{}_{}_{}'.format(zfig,2,0)].xy[1])
                g['an_{}_{}_{}'.format(zfig,2,0)].set_c(plot_preferences.plot_pre()[0][3] if Cratemax<res[rn_min]['Crate_max'][0] else plot_preferences.plot_pre()[0][0])
            except:
                0
       
    # limits
    # g['pl_{}_{}_{}'.format(zfig,1,2)][0].set_data(xlim,[res_recalc[rn_min_recalc]['U_min'][0]]*2)
    # g['pl_{}_{}_{}'.format(zfig,1,3)][0].set_data(xlim,[res_recalc[rn_min_recalc]['U_max'][0]]*2)
    # g['pl_{}_{}_{}'.format(zfig,2,2)][0].set_data(xlim,[res_recalc[rn_min_recalc]['Crate_max'][0]]*2)
   
    plt.draw()    