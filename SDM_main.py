"""
System Design Method Interactive Tool (SDM_main.py)

Author:     Sven Wiegelmann,
            Leibniz University Hannover,
            Institute of Electric Power Systems,
            Electric Energy Storage Systems Section,
            Appelstraße 9A,
            30167 Hannover,
            Germany
Version:    12.12.2023

Overview:
SDM_main.py is the core interactive tool of the System Design Method (SDM) suite,
integrating functionalities from SDM_fcn and SDM_data. It employs matplotlib
widgets for interactive visualization and manipulation of energy storage system
parameters, offering a dynamic user experience in system design.

Features:
- Interactive Extended Ragone Plot (ERP) generation and analysis.
- Visualization of constraints and performance metrics in energy systems.
- Dynamic adjustment of system parameters such as power, energy, and efficiency.
- Integrated case studies for practical understanding of energy system designs.
- Graphical synthesis for decision-making in system design parameters.

Usage:
Designed for energy system researchers and practitioners, this tool requires
Python programming knowledge and familiarity with energy storage concepts.
Users can interactively modify and explore energy system designs.

Execution:
Run this script in a Python environment with required dependencies (matplotlib,
numpy, etc.). Ensure SDM_fcn.py and SDM_data.py are in the same directory.

Notes:
- Part of the SDM framework, intended for educational and research use in energy
  system design.
- Prior knowledge of the System Design Method principles is recommended for
  effective use.
- For detailed function information, refer to SDM_fcn.py and SDM_data.py docstrings.
"""

import matplotlib as mpl
mpl.use('Qt5Agg')
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider#, TextBox, Button, RadioButtons, RangeSlider

import numpy as np

import sys
sys.path.append('./fcn')
import os
cpath = os.path.dirname(os.path.abspath(__file__))

import time
from IPython import get_ipython

tic_all = time.time()
get_ipython().magic('reset -sf')
plt.close('all')

import plot_preferences
import SDM_fcn
import SDM_data


if __name__ == "__main__":
    # %%
    print('##################################################################')
    print('PREPROCESSING')
    """   
    This section initializes and configures the primary settings and parameters for 
    the System Design Method (SDM) analysis. It includes the following key steps:
    
    - Basic Settings: Establishes initial configurations such as data prefixes, 
      unit conversions, and simulation limits.
    - Data Loading: Utilizes SDM_data to load and prepare cell-based experimental 
      data for energy storage systems.
    - Parameter Evaluation: Leverages SDM_fcn functions to compile Extended Ragone 
      Plot (ERP) data and recalculates system limits based on specific parameters.
    - Component Settings: Defines settings for different case studies and components
      like application requirements, power electronics converters, and cell modules.
    
    This preprocessing is critical for setting up the system's foundation, ensuring 
    accurate and efficient analysis in subsequent stages of the SDM tool.
    """

    tic = time.time()

    # Basic Settings
    prefix_str = ''
    convert_unit_E = True
    sim_limits = False
    calc_EPsys = False
 
    num_SoC = 11
    SoC_set = np.linspace(1,0,num_SoC)
    

    # Input: Cell-based (Measurement) Data of Energy Storage System
    res = SDM_data.load_results(0,prefix_str=prefix_str)    
    
    # Evaluation of Parameters for ERP Generation 
    d, tmp_ind, rn_min = SDM_fcn.compile_EP(res,calc_EPsys=calc_EPsys,prefix_str=prefix_str)
    res_recalc = SDM_fcn.recalc_limits(res, Umax=res[rn_min]['U_max'][0], Umin=res[rn_min]['U_min'][0], Cratemax=res[rn_min]['Crate_max'][0], Tmax=res[rn_min]['T_max'][0],prefix_str=prefix_str)
    d_recalc, tmp_ind_recalc, rn_min_recalc = SDM_fcn.compile_EP(res_recalc,calc_EPsys=calc_EPsys,prefix_str=prefix_str)
    
    
    # Component settings based on General Design Topology
    def set_case(no):
        """
        Set system analysis parameters for different case studies in the System 
        Design Method (SDM). Adjusts settings based on the case study, including 
        application requirements, power electronics converter settings, and energy 
        storage cell configurations.
    
        Parameters:
        - no (int): Identifier for the specific case study setup.
    
        Functionality:
        - Initializes settings for various case studies, adapting system design 
          parameters.
        - Manages configurations for application requirements, power electronics 
          converters, and energy storage cells.
        - Dynamically updates design parameters to align with selected case scenarios.
    
        Note:
        The function is tailored for the SDM framework. Effectiveness depends on 
        the consistency and format of input data related to energy storage cells.
        """
    
        ## Initial Settings for Plot Generation
        if (no==-1):
            c['req'] = SDM_data.set_app(0)    
            c['pec'] = SDM_data.set_pec(0)
            c['cell'] = ({'Udc_max':        res[rn_min]['U_max'][0],
                          'Udc_min':        res[rn_min]['U_min'][0],
                          'Idc_DIS_max':    res[rn_min]['I_dis_max'][0],
                          'Pdc_DIS_max':    res[rn_min]['U_max'][0]*res[rn_min]['I_dis_max'][0]})
            c['mod'] =({'xS': 1,
                        'yP': 1})
        
        ## Case Studies
        if (no==0): # DSS Explanation
            c['req'] = SDM_data.set_app(0)
            # c['pec'] = set_pec(0)
            c['pec']['n'] = 2
            
            c['mod']['xS'] = 4
            c['mod']['yP'] = 2
            
            update_dp(c['req']['E/P_req'],270,500,1800,550000,100)
            
            # Hide Ticklabels
            plot_preferences.hide_ticklabels(g['ax_{}_{}'.format(3,0)].xaxis)
            plot_preferences.hide_ticklabels(g['ax_{}_{}'.format(3,0)].yaxis)
            
            plot_preferences.hide_ticklabels(g['cbar_{}_{}'.format(3,0)].ax.xaxis)
            g['cbar_{}_{}'.format(3,1)].set_xticklabels([r'$n_{{x_{{\mathrm{{S}}}}y_{{\mathrm{{P}}}},\mathrm{{opt}}}}$'])
        
        
        if (no==1): # DSS available
            c['req'] = SDM_data.set_app(1)
            c['pec'] = SDM_data._set_by_arg(SDM_data.set_pec,'name','SMA Solar SCS800')
            c['pec']['n'] = np.ceil(c['req']['Pac_req']/c['pec']['Pac_max'])
        
            c['pec']['eta_op'] = SDM_fcn.update_eta(g,c)
        
            c['mod']['xS'] = 2
            c['mod']['yP'] = 1
            
            update_dp(c['req']['E/P_req'],
                   c['pec']['Udc_min'],
                   c['pec']['Udc_max'],
                   c['pec']['Idc_max'],
                   c['pec']['Pdc_max'])
            
        if no in [2,3,4]: # no DSS available --> Limit Variation
            c['req'] = SDM_data.set_app(2)
            
            c['pec'] = SDM_data._set_by_arg(SDM_data.set_pec,'name','SMA Solar SCS800')
            c['pec']['n'] = 1 #np.ceil(c['req']['Pac_req']/c['pec']['Pac_max'])
        
            c['mod']['xS'] = 2
            c['mod']['yP'] = 2
            
            update_dp(c['req']['E/P_req'],
                   c['pec']['Udc_min'],
                   c['pec']['Udc_max'],
                   c['pec']['Idc_max'],
                   c['pec']['Pdc_max'],
                   0 if (no==2) else 138 if (no==3) else 133)

        if (no==5): # no DSS ever available
            c['req'] = SDM_data.set_app(3)
            
            c['pec'] = SDM_data._set_by_arg(SDM_data.set_pec,'name','SMA Solar SCS800')
            c['pec']['n'] = 1 #np.ceil(c['req']['Pac_req']/c['pec']['Pac_max'])
        
            c['mod']['xS'] = 1
            c['mod']['yP'] = 1
            
            update_dp(c['req']['E/P_req'],
                   c['pec']['Udc_min'],
                   c['pec']['Udc_max'],
                   c['pec']['Idc_max'],
                   c['pec']['Pdc_max'])
            
        if no in [6,7]: # Methodology Comparison
            c['req'] = SDM_data.set_app(3)
            
            c['pec'] = SDM_data._set_by_arg(SDM_data.set_pec,'name','SMA Solar SCS800')
            c['pec']['n'] = 1 #np.ceil(c['req']['Pac_req']/c['pec']['Pac_max'])
        
            c['mod']['xS'] = 2
            c['mod']['yP'] = 1
            
            Ecell_nom = 78.4
            Pcell_nom = 2.21 * 198
            
            c['req']['E/P_req'] = Ecell_nom/Pcell_nom
                        
            update_dp(c['req']['E/P_req'],
                    c['pec']['Udc_min'],
                    c['pec']['Udc_max'],
                    c['pec']['Idc_max'],
                    c['pec']['Pdc_max'])
    
    global c
    c = {}
    
    set_case(-1)
    

    #%%
    print('##################################################################')
    print('POSTPROCESSOR')
    """
    This section of the script is dedicated to visualizing and interactively exploring 
    the System Design Method (SDM) analysis results. Key activities include:
    
    - Generating Plots:
      - Charge/Discharge Curves: Displays voltage and current vs. energy.
      - Extended Ragone Plot (ERP): Shows intersections of energy and power parameters.
      - Design Solution Space (DSS): Visualizes the feasible design options within the 
        constraint satisfaction problem framework.
      - Graphical Synthesis: Determines the availability of a design solution space.
      - Efficiency Curve of Power Electric Converter (PEC): Illustrates efficiency 
        variations with power ratios.
    
    - Embedding Slider Widget:
      - An interactive slider widget is integrated for dynamic parameter adjustment.
      - It allows real-time exploration of the impact of energy-to-power ratio (E/P) 
        and maximum cell power on system design.
      - This feature enables users to engage in an intuitive and dynamic exploration 
        of the SDM results, enhancing the understanding of system performance.
    
    These components are instrumental in providing a comprehensive and interactive 
    understanding of the system's capabilities and potential design considerations.
    """

    # Plot Preferences
    farbe, z, g, fig_size, lstyle = plot_preferences.plot_pre()        

    #%%
    # Plot: Charge/Discharge Curves U(E), I(E), Derivation of the Extended Ragone Plot (ERP)
    # -------------------------------------------------------------------------  
    z = 1
        
    nPref = 15
    pub_cut_Umin = 1.95
    pub_cut_Cratemax = 4.25

    grid = plt.GridSpec(1,3, wspace=0.55, hspace=0.05)
    if 'fig_{}'.format(z) not in g:
        g['fig_{}'.format(z)] = plt.figure(z,figsize=[fig_size[-1]*2.95,fig_size[-1]*1.05])


    # ax 0
    if 'ax_{}_{}'.format(z,0) not in g:
        g['ax_{}_{}'.format(z,0)] = g['fig_{}'.format(z)].add_subplot(grid[:,2]) 
        # g['ax_{}_{}'.format(z,0)].grid()
    
    nSoC = 0
    nUmin = 0
    nUmax = 0
    nCratemax = 0
    nTmax = 0
    nSoCmin = 0
    n0 = 0
    alph = 1/3
   
    
    g['pl_{}_{}_{}_{}'.format(z,0,0,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d['p'],d[nSoC][0],color=farbe[nSoC],marker='.',markersize=4,label='{}'.format(round(SoC_set[nSoC],1)),lw=0.75,alpha=alph)
    # g['pl_{}_{}_{}_{}'.format(z,0,1,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d['p'],d[nSoC][1],color=farbe[nSoC],marker='v',markersize=4,ls='None',alpha=alph)
    # g['pl_{}_{}_{}_{}'.format(z,0,2,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d['p'],d[nSoC][2],color=farbe[nSoC],marker='o',markersize=4,ls='None',alpha=alph)
    # g['pl_{}_{}_{}_{}'.format(z,0,3,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d['p'],d[nSoC][3],color=farbe[nSoC],marker='d',markersize=4,ls='None',alpha=alph)
    # g['pl_{}_{}_{}_{}'.format(z,0,4,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d['p'],d[nSoC][4],color=farbe[nSoC],marker='s',markersize=4,ls='None',alpha=alph)
   
   
    g['pl_{}_{}_{}_{}'.format(z,0,10,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d_recalc['p'],d_recalc[nSoC][0],color=farbe[nSoC],marker='.',markersize=6,label='{}'.format(round(SoC_set[nSoC],1)),lw=1.5)
    # g['pl_{}_{}_{}_{}'.format(z,0,11,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d_recalc['p'],d_recalc[nSoC][1],color=farbe[nSoC],marker='v',markersize=6,ls='None')
    # g['pl_{}_{}_{}_{}'.format(z,0,12,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d_recalc['p'],d_recalc[nSoC][2],color=farbe[nSoC],marker='o',markersize=6,ls='None')
    # g['pl_{}_{}_{}_{}'.format(z,0,13,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d_recalc['p'],d_recalc[nSoC][3],color=farbe[nSoC],marker='d',markersize=6,ls='None')
    # g['pl_{}_{}_{}_{}'.format(z,0,14,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d_recalc['p'],d_recalc[nSoC][4],color=farbe[nSoC],marker='s',markersize=6,ls='None')
   
    try:
        g['pl_{}_{}_{}'.format(z,0,10)] = g['ax_{}_{}'.format(z,0)].plot(d_recalc['p'][nPref],d_recalc[nSoC][0][nPref],color=farbe[2],marker='o',markersize=6,lw=1.5,zorder=100)
    except:
        0
    
    
    # P_UI
    res[rn_min]['P_UI'],res[rn_min]['E_UI'] = SDM_fcn.calc_P_UI(d_recalc['p'],d[nSoC][0],res[rn_min]['U_min'][0],res[rn_min]['I_dis_max'][0])
    
    g['pl_{}_{}_{}'.format(z,0,100)] = g['ax_{}_{}'.format(z,0)].plot([res[rn_min]['P_UI']]*2,[0,1e3],color=farbe[12],lw=1,zorder=0)
    g['pl_{}_{}_{}'.format(z,0,101)] = g['ax_{}_{}'.format(z,0)].plot(res[rn_min]['P_UI'],res[rn_min]['E_UI'],color=farbe[12],marker='o',markersize=8,mfc='None',lw=1,zorder=100)
    
    bbox_args = dict(boxstyle="round",color='None',fc='1',lw=0)
    g['an_{}_{}_{}'.format(z,0,100)] = g['ax_{}_{}'.format(z,0)].annotate(r'$P^\mathrm{dis}_\mathrm{cell,UI}$',
                                      (res[rn_min]['P_UI'],0.95*res[rn_min]['lim_E_dis_max']),
                                      textcoords="offset points",
                                      xytext=(0,0),
                                      ha='center',
                                      va='center',
                                      c=farbe[12],
                                      bbox=bbox_args,
                                      zorder=100)
    
    
    # Umin
    res_recalc = SDM_fcn.recalc_limits(res, Umax=res[rn_min]['U_max'][0], Umin=pub_cut_Umin, Cratemax=res[rn_min]['Crate_max'][0], Tmax=res[rn_min]['T_max'][0],prefix_str=prefix_str)
    d_recalc, tmp_ind_recalc, rn_min_recalc = SDM_fcn.compile_EP(res_recalc,calc_EPsys=calc_EPsys,prefix_str=prefix_str)
    g['pl_{}_{}_{}_{}'.format(z,0,20,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d_recalc['p'],d_recalc[nSoC][0],color=farbe[5],marker='.',markersize=4,lw=1)
   
    # Imax
    res_recalc = SDM_fcn.recalc_limits(res, Umax=res[rn_min]['U_max'][0], Umin=res[rn_min]['U_min'][0], Cratemax=pub_cut_Cratemax, Tmax=res[rn_min]['T_max'][0],prefix_str=prefix_str)
    d_recalc, tmp_ind_recalc, rn_min_recalc = SDM_fcn.compile_EP(res_recalc,calc_EPsys=calc_EPsys,prefix_str=prefix_str)
    g['pl_{}_{}_{}_{}'.format(z,0,21,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d_recalc['p'],d_recalc[nSoC][0],color=farbe[3],marker='.',markersize=4,lw=1)
   
    # Umin & Imax
    res_recalc = SDM_fcn.recalc_limits(res, Umax=res[rn_min]['U_max'][0], Umin=pub_cut_Umin, Cratemax=pub_cut_Cratemax, Tmax=res[rn_min]['T_max'][0],prefix_str=prefix_str)
    d_recalc, tmp_ind_recalc, rn_min_recalc = SDM_fcn.compile_EP(res_recalc,calc_EPsys=calc_EPsys,prefix_str=prefix_str)
    g['pl_{}_{}_{}_{}'.format(z,0,22,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d_recalc['p'],d_recalc[nSoC][0],color=farbe[0],marker='.',markersize=6,lw=1.5)
   
   
    # Legende
    # if rn == lres[0]:
    g['leg_{}'.format(z)] = {}
    handles, labels = g['ax_{}_{}'.format(z,0)].get_legend_handles_labels()
    g['leg_{}'.format(z)].update(dict(zip(labels, handles)))
            
    # Labeling
    g['ax_{}_{}'.format(z,0)].set_xlabel('Cell power, $P^\mathrm{dis}_\mathrm{cell}$ in W')
    g['ax_{}_{}'.format(z,0)].set_xlim([0,70 if sim_limits else res[rn_min]['lim_P_dis_max']])
    g['ax_{}_{}'.format(z,0)].set_ylabel('Cell energy, $E^\mathrm{{dis}}_\mathrm{{cell}}$ in {}'.format('Wh' if convert_unit_E else 'Ws'))
    g['ax_{}_{}'.format(z,0)].set_ylim([0,10 if sim_limits else res[rn_min]['lim_E_dis_max']])
    g['ax_{}_{}'.format(z,0)].set_aspect(1./g['ax_{}_{}'.format(z,0)].get_data_ratio())
    SDM_fcn._plot_isochrones(g['ax_{}_{}'.format(z,0)],num_subplot=True)


    # ax 1
    if 'ax_{}_{}'.format(z,1) not in g:
        g['ax_{}_{}'.format(z,1)] = g['fig_{}'.format(z)].add_subplot(grid[:,0])
        # g['ax_{}_{}'.format(z,1)].grid()
    
    
    for nSoC in range(max(tmp_ind)[2]+1):
        if nSoC > 0:
            continue
        for nP in range(max(tmp_ind)[1]+1):
            rn = [item[0] for item in tmp_ind if item[1]==nP and item[2]==nSoC][0]
            try:
                if np.mod(nP,3):# "thinning"
                    continue
                alph = 1/3
                g['pl_{}_{}_{}_{}'.format(z,1,0,rn)] = g['ax_{}_{}'.format(z,1)].plot(res[rn][prefix_str+'E_cell'],res[rn][prefix_str+'U_cell'],color=farbe[nSoC],lw=0.75,alpha=alph)
                g['pl_{}_{}_{}_{}'.format(z,1,1,rn)] = g['ax_{}_{}'.format(z,1)].plot(res_recalc[rn][prefix_str+'E_cell'],res_recalc[rn][prefix_str+'U_cell'],color=farbe[nSoC],lw=1.5)
                if nP == nPref:
                    g['pl_{}_{}_{}_{}'.format(z,1,1,rn)][0].remove()
                    g['pl_{}_{}_{}_{}'.format(z,1,1,rn)] = g['ax_{}_{}'.format(z,1)].plot(res_recalc[rn][prefix_str+'E_cell'],res_recalc[rn][prefix_str+'U_cell'],color=farbe[2],lw=1.5,zorder=100)
            except:
                0
    
    xlim = [0,85]
    
    g['ax_{}_{}'.format(z,1)].plot(xlim,[res[rn_min]['U_min'][0]]*2,linestyle='dashed',color=[0]*3,lw=1,zorder=1)
    g['ax_{}_{}'.format(z,1)].plot(xlim,[res[rn_min]['U_max'][0]]*2,linestyle='dashed',color=[0]*3,lw=1,zorder=1)

    # Annotations
    bbox_args = dict(boxstyle="round",color='None',fc='1',lw=0)
    g['an_{}_{}_{}'.format(z,1,0)] = g['ax_{}_{}'.format(z,1)].annotate(r'$U_\mathrm{cell,min}$', # this is the text
                                        (0.1*xlim[-1],res[rn_min]['U_min'][0]), # these are the coordinates to position the label
                                        # textcoords="offset points", # how to position the text
                                        # xytext=(0,0), # distance from text to points (x,y)
                                        # va='bottom',
                                        ha='center',
                                        va='center',
                                        # fontsize='large',
                                        bbox=bbox_args,
                                        arrowprops = dict(arrowstyle="<|-",color=farbe[5]),
                                        zorder=102)
    # bbox_args = dict(boxstyle="round",color='None',fc='1',lw=0)
    g['an_{}_{}_{}'.format(z,1,1)] = g['ax_{}_{}'.format(z,1)].annotate(r'$U_\mathrm{cell,max}$',
                                        (0.1*xlim[-1],res[rn_min]['U_max'][0]),
                                        # textcoords="offset points",
                                        # xytext=(0,0),
                                        ha='center',
                                        va='center',
                                        # fontsize='large',
                                        bbox=bbox_args,
                                        arrowprops = dict(arrowstyle="<|-",color=farbe[5]),
                                        zorder=102)
    
    g['pl_{}_{}_{}'.format(z,1,0)] = g['ax_{}_{}'.format(z,1)].plot(xlim,[pub_cut_Umin]*2,color=farbe[5],lw=1,ls='dashed',zorder=101)

            
    # Labeling
    # g['ax_{}_{}'.format(z,i)].set_title('Simulation {}'.format(len(res)))
    # plot_preferences.hide_ticklabels(g['ax_{}_{}'.format(z,1)].xaxis)
    g['ax_{}_{}'.format(z,1)].set_xlabel('Cell energy, $E^\mathrm{{dis}}_\mathrm{{cell}}$ in {}'.format('Wh' if convert_unit_E else 'Ws'))
    g['ax_{}_{}'.format(z,1)].set_xlim(xlim)
    g['ax_{}_{}'.format(z,1)].set_ylabel('Cell voltage, $U_\mathrm{cell}$ in V')
    g['ax_{}_{}'.format(z,1)].set_aspect(1./g['ax_{}_{}'.format(z,1)].get_data_ratio())

    
    # ax 2
    if 'ax_{}_{}'.format(z,2) not in g:
        g['ax_{}_{}'.format(z,2)] = g['fig_{}'.format(z)].add_subplot(grid[:,1])
        # g['ax_{}_{}'.format(z,2)].grid()
    
    
    for nSoC in range(max(tmp_ind)[2]+1):
        if nSoC > 0:
            continue
        for nP in range(max(tmp_ind)[1]+1):
            rn = [item[0] for item in tmp_ind if item[1]==nP and item[2]==nSoC][0]
            try:
                if np.mod(nP,3):# "thinning"
                    continue
                alph = 1/3
                g['pl_{}_{}_{}_{}'.format(z,2,0,rn)] = g['ax_{}_{}'.format(z,2)].plot(res[rn][prefix_str+'E_cell'],res[rn][prefix_str+'I_cell'],color=farbe[nSoC],lw=0.75,alpha=0.33)
                g['pl_{}_{}_{}_{}'.format(z,2,1,rn)] = g['ax_{}_{}'.format(z,2)].plot(res_recalc[rn][prefix_str+'E_cell'],res_recalc[rn][prefix_str+'I_cell'],color=farbe[nSoC],lw=1.5)
    
                if nP == 0:
                    g['pl_{}_{}_{}_{}'.format(z,2,0,rn)][0].set_label('{}'.format(round(SoC_set[nSoC],1)))
                if nP == nPref:
                    g['pl_{}_{}_{}_{}'.format(z,2,1,rn)][0].remove()
                    g['pl_{}_{}_{}_{}'.format(z,2,1,rn)] = g['ax_{}_{}'.format(z,2)].plot(res_recalc[rn][prefix_str+'E_cell'],res_recalc[rn][prefix_str+'I_cell'],color=farbe[2],lw=1.5,zorder=100)                   
            except:
                0
    
    g['ax_{}_{}'.format(z,2)].plot(xlim,[0]*2,linestyle='dashed',color=[0]*3,lw=1,zorder=1)
    g['ax_{}_{}'.format(z,2)].plot(xlim,[res[rn_min]['I_dis_max'][0]]*2,linestyle='dashed',color=[0]*3,lw=1,zorder=1)
    
    # Annotations
    bbox_args = dict(boxstyle="round",color='None',fc='1',lw=0)
    g['an_{}_{}_{}'.format(z,2,0)] = g['ax_{}_{}'.format(z,2)].annotate(r'$I^\mathrm{dis}_\mathrm{cell,max}$',
                                      (0.1*xlim[-1],res[rn_min]['I_dis_max'][0]),
                                      # textcoords="offset points",
                                      # xytext=(2,0),
                                      ha='center',
                                      va='center',
                                      # fontsize='large',
                                      bbox=bbox_args,
                                      arrowprops = dict(arrowstyle="<|-",color=farbe[3]),
                                      zorder=102)
    # set zorder of 2Dline above annotation
    g['pl_{}_{}_{}_{}'.format(z,2,1,36)][0].set_zorder(g['an_{}_{}_{}'.format(z,2,0)].get_zorder()+1)

    
    # g['ax_{}_{}'.format(z,2)].plot(xlim,[res[rn_min]['Crate_max'][0]]*2,linestyle='dashed',color=[0]*3,lw=1.5,zorder=1)
    g['pl_{}_{}_{}'.format(z,2,0)] = g['ax_{}_{}'.format(z,2)].plot(xlim,[pub_cut_Cratemax*res[rn_min]['_settings']['ParameterValues']['Q_n'][0]]*2,color=farbe[3],ls='dashed',lw=1,zorder=101)

    # Labeling
    # g['ax_{}_{}'.format(z,i)].set_title('Simulation {}'.format(len(res)))
    g['ax_{}_{}'.format(z,2)].set_xlabel('Cell energy, $E^\mathrm{{dis}}_\mathrm{{cell}}$ in {}'.format('Wh' if convert_unit_E else 'Ws'))
    g['ax_{}_{}'.format(z,2)].set_xlim(xlim)
    g['ax_{}_{}'.format(z,2)].set_ylabel('Cell current, $I^\mathrm{dis}_\mathrm{cell}$ in A')
    g['ax_{}_{}'.format(z,2)].set_aspect(1./g['ax_{}_{}'.format(z,2)].get_data_ratio())
    
    
    # Annotations
    plot_preferences.print_num_subplot(g['ax_{}_{}'.format(z,0)],'c')
    plot_preferences.print_num_subplot(g['ax_{}_{}'.format(z,1)],'a')
    plot_preferences.print_num_subplot(g['ax_{}_{}'.format(z,2)],'b')

    g['an_{}_{}_{}'.format(z,1,100)] = g['ax_{}_{}'.format(z,1)].annotate(r'$P^\mathrm{dis}_\mathrm{cell} \uparrow$',
                                      (60,2.225),
                                      # textcoords="offset points",
                                      xytext=(30,1.925),
                                      ha='center',
                                      va='center',
                                      color=[0.5]*3,
                                      # fontsize='large',
                                      bbox=bbox_args,
                                      arrowprops = dict(arrowstyle="<|-",color=[0.5]*3,connectionstyle='angle3,angleA=5,angleB=75'),
                                      zorder=0)
    g['an_{}_{}_{}'.format(z,2,100)] = g['ax_{}_{}'.format(z,2)].annotate(r'$P^\mathrm{dis}_\mathrm{cell} \uparrow$',
                                      (80,80),
                                      # textcoords="offset points",
                                      xytext=(76.5,160),
                                      ha='center',
                                      va='center',
                                      color=[0.5]*3,
                                      # fontsize='large',
                                      bbox=bbox_args,
                                      arrowprops = dict(arrowstyle="<|-",color=[0.5]*3,connectionstyle='angle3,angleA=90,angleB=100'),
                                      zorder=0)

    grid.tight_layout(g['fig_{}'.format(z)])
    # g['fig_{}'.format(z)].canvas.draw()

    SDM_fcn.replot_EP_limits(g,res,tmp_ind,Umax=5.0,Umin=1.0,Cratemax=10,Tmax=400,zfig=1,nPref=nPref,calc_EPsys=calc_EPsys,print_ann=True,prefix_str=prefix_str)


    #%%
    # Plot: Extended Ragone Plot (ERP) - Variation of Limits
    # -------------------------------------------------------------------------    
    z = 2
    bbox_args = dict(boxstyle="round",color='1',fc='1',ec='None',lw=0)

    grid = plt.GridSpec(2,2, wspace=0.25, hspace=0.05)
    if 'fig_{}'.format(z) not in g:
            g['fig_{}'.format(z)] = plt.figure(z,figsize=[fig_size[-1]*2.05,fig_size[-1]])
       
    if 'ax_{}_{}'.format(z,0) not in g:
        g['ax_{}_{}'.format(z,0)] = g['fig_{}'.format(z)].add_subplot(grid[:,:-1]) 
        # g['ax_{}_{}'.format(z,0)].grid()
   
    SDM_fcn.plot_EP_limit(g,z,0, res, rn_min, prefix_str=prefix_str, plot_limit=['Umin','Cratemax'])
   
    # Labeling
    g['ax_{}_{}'.format(z,0)].set_xlabel('Cell power, $P^\mathrm{{dis}}_\mathrm{cell}$ in W')
    g['ax_{}_{}'.format(z,0)].set_xlim([0,70 if sim_limits else res[rn_min]['lim_P_dis_max']])
    g['ax_{}_{}'.format(z,0)].set_ylabel('Cell energy, $E^\mathrm{{dis}}_\mathrm{{cell}}$ in {}'.format('Wh' if convert_unit_E else 'Ws'))
    g['ax_{}_{}'.format(z,0)].set_ylim([0,10 if sim_limits else res[rn_min]['lim_E_dis_max']])
    g['ax_{}_{}'.format(z,0)].set_aspect(1./g['ax_{}_{}'.format(z,0)].get_data_ratio())
    SDM_fcn._plot_isochrones(g['ax_{}_{}'.format(z,0)],num_subplot=True)
   
    # SDM_fcn.recalc_design(res, Umax=3.5, Umin=2.0, Cratemax=5, Tmax=400, m=c['req']['E/P_req'], show_annotation=False)



    # Plot: Intersection Curves of Operational Quantities (by Limit Variation) 
    # -------------------------------------------------------------------------        
    dp_intersection = SDM_fcn.calc_intersection_curve(g,res,c['req']['E/P_req'],print_error=False,prefix_str=prefix_str)
    g['pl_{}_{}_{}'.format(z,0,100)] = g['ax_{}_{}'.format(z,0)].plot(dp_intersection['valx_P'],dp_intersection['valy_E'],color=farbe[2])
    SDM_fcn._set_intersection_limits(g,dp_intersection['valx_P'],dp_intersection['valy_E'],z,0)


    # ax1
    if 'ax_{}_{}'.format(z,1) not in g:
        g['ax_{}_{}'.format(z,1)] = g['fig_{}'.format(z)].add_subplot(grid[0,-1:]) 
        # g['ax_{}_{}'.format(z,1)].grid()
    
    try: 
        g['pl_{}_{}_{}'.format(z,1,100)] = g['ax_{}_{}'.format(z,1)].plot(dp_intersection['valx_P'],dp_intersection[prefix_str+'U_cell'],color=farbe[5])
        
        id_max = max([n for n,i in enumerate(dp_intersection[prefix_str+'U_cell']) if ~np.isnan(i)])
        yint = np.interp(dp_intersection['valx_P'][id_max],dp_intersection['valx_P'],dp_intersection[prefix_str+'U_cell'])
        g['pl_{}_{}_{}'.format(z,1,101)] = g['ax_{}_{}'.format(z,1)].plot(dp_intersection['valx_P'][id_max],dp_intersection[prefix_str+'U_cell'][id_max],color=farbe[0],marker='.',markersize=8)

        g['pl_{}_{}_{}'.format(z,1,110)] = g['ax_{}_{}'.format(z,1)].plot([dp_intersection['valx_P'][id_max]]*2,[-1e20,1e20],color=farbe[2],ls='dotted',lw=1)
        g['pl_{}_{}_{}'.format(z,1,111)] = g['ax_{}_{}'.format(z,1)].plot(dp_intersection['valx_P'][id_max],yint,color=farbe[2],ls='None',marker='o',markersize=4.5)     
    except:
        0
    
    plot_preferences.hide_ticklabels(g['ax_{}_{}'.format(z,1)].xaxis)
    g['ax_{}_{}'.format(z,1)].set_ylabel('Minimum cell voltage,\n$U_\mathrm{min}$ in V')
    ylim = plot_preferences._set_lim(res[rn_min]['U_min'][0],res[rn_min]['U_max'][0],0.1)
    g['ax_{}_{}'.format(z,1)].set_ylim(ylim)
    
    # Annotations
    bbox_args = dict(boxstyle="round",color='None',fc='1',lw=0)
    g['an_{}_{}_{}'.format(z,1,110)] = g['ax_{}_{}'.format(z,1)].annotate(r'$P^\mathrm{dis}_\mathrm{cell,max}$',
                                      (dp_intersection['valx_P'][id_max],ylim[1]),
                                      textcoords="offset points",
                                      xytext=(0,-8),
                                      ha='center',
                                      va='top',
                                      c=farbe[2],
                                      bbox=bbox_args,
                                      zorder=100)
    
    
    # ax1
    if 'ax_{}_{}'.format(z,2) not in g:
        g['ax_{}_{}'.format(z,2)] = g['fig_{}'.format(z)].add_subplot(grid[1,-1:],sharex=g['ax_{}_{}'.format(z,1)]) 
        # g['ax_{}_{}'.format(z,2)].grid()
    
    try:
        g['pl_{}_{}_{}'.format(z,2,100)] = g['ax_{}_{}'.format(z,2)].plot(dp_intersection['valx_P'],dp_intersection[prefix_str+'I_cell'],color=farbe[3])
        
        id_max = max([n for n,i in enumerate(dp_intersection[prefix_str+'I_cell']) if ~np.isnan(i)])
        yint = np.interp(dp_intersection['valx_P'][id_max],dp_intersection['valx_P'],dp_intersection[prefix_str+'I_cell'])
        g['pl_{}_{}_{}'.format(z,2,101)] = g['ax_{}_{}'.format(z,2)].plot(dp_intersection['valx_P'][id_max],dp_intersection[prefix_str+'I_cell'][id_max],color=farbe[0],marker='.',markersize=8)
        
        g['pl_{}_{}_{}'.format(z,2,110)] = g['ax_{}_{}'.format(z,2)].plot([dp_intersection['valx_P'][id_max]]*2,[-1e20,1e20],color=farbe[2],ls='dotted',lw=1)
        g['pl_{}_{}_{}'.format(z,2,111)] = g['ax_{}_{}'.format(z,2)].plot(dp_intersection['valx_P'][id_max],yint,color=farbe[2],ls='None',marker='o',markersize=4.5)        
    except:
        0
    
    g['ax_{}_{}'.format(z,2)].set_xlabel('Cell power, $P^\mathrm{dis}_\mathrm{cell}$ in W')
    g['ax_{}_{}'.format(z,2)].set_xlim([0,70 if sim_limits else res[rn_min]['lim_P_dis_max']])
    g['ax_{}_{}'.format(z,2)].set_ylabel('Maximum cell current,\n$I^\mathrm{dis}_\mathrm{cell,max}$ in A')
    ylim = plot_preferences._set_lim(0,30 if sim_limits else res[rn_min]['I_dis_max'][0],0.1)
    g['ax_{}_{}'.format(z,2)].set_ylim(ylim)
    grid.tight_layout(g['fig_{}'.format(z)])
    
    

    #%%
    # Plot: Design Solution Space (DSS) of the underlying Constraint Satisfaction Problem (CSP)
    # -------------------------------------------------------------------------
    z = 3
    grid = plt.GridSpec(1, 1, wspace=0.55, hspace=0.05)
    if 'fig_{}'.format(z) not in g:
        g['fig_{}'.format(z)] = plt.figure(z,figsize=[fig_size[-1],fig_size[-1]*1.2])
        
    if 'ax_{}_{}'.format(z,0) not in g:
        g['ax_{}_{}'.format(z,0)] = g['fig_{}'.format(z)].add_subplot(grid[0,-1]) 
        # g['ax_{}_{}'.format(z,1)].grid()

    SDM_fcn.plot_dss(g,c,zfig=z,zax=0)

    g['ax_{}_{}'.format(z,0)].set_xlabel('Number of cells in series, $x_\mathrm{S}$')
    g['ax_{}_{}'.format(z,0)].set_ylabel('Number of cells in parallel, $y_\mathrm{P}$')
    grid.tight_layout(g['fig_{}'.format(z)])



    #%%
    # Plot: Graphical Synthesis of the Results - Will a DSS be available?
    # -------------------------------------------------------------------------
    z = 4
    iso_mins = np.array([10/60,1,2.5,5,7.5,10,15,30,60,300])/60
    alph = 0.2
    
    grid = plt.GridSpec(2,2, wspace=0.25, hspace=0.05)
    if 'fig_{}'.format(z) not in g:
            g['fig_{}'.format(z)] = plt.figure(z,figsize=[fig_size[-1]*2.05,fig_size[-1]])
      
    # ax 0 
    if 'ax_{}_{}'.format(z,0) not in g:
        g['ax_{}_{}'.format(z,0)] = g['fig_{}'.format(z)].add_subplot(grid[0,0])
        # g['ax_{}_{}'.format(z,0)].grid()
   

    # Operation
    set_EP = c['req']['E/P_req']
    set_P = c['cell']['Pdc_DIS_max']
    dp_intersection = SDM_fcn.calc_intersection_curve(g,res,set_EP,print_error=False,prefix_str=prefix_str)
    EPline_x = dp_intersection['valx_P']
    EPline_y = dp_intersection[prefix_str+'P_cell']
    id_min = min([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
    id_max = max([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
    sPmaxset = EPline_x[id_min] if set_P < EPline_x[id_min] else (EPline_x[id_max] if set_P > EPline_x[id_max] else set_P)
    
    g['pl_{}_{}_{}'.format(z,0,0)] = g['ax_{}_{}'.format(z,0)].plot(dp_intersection[prefix_str+'P_cell'],dp_intersection[prefix_str+'U_cell'],color=farbe[5],zorder=100)
    calc_P,sUminset = SDM_fcn._set_intersection_limits(g,EPline_x,dp_intersection[prefix_str+'U_cell'],z,0,x=sPmaxset)

    SDM_fcn.plot_iso_intersection_curves(g,z,0,c,res,prefix_str=prefix_str,iso_mins=iso_mins,var='Umin',print_ann=False)
        
   
    # Labeling
    plot_preferences.hide_ticklabels(g['ax_{}_{}'.format(z,0)].xaxis)
    g['ax_{}_{}'.format(z,0)].set_xlim([0,70 if sim_limits else res[rn_min]['lim_P_dis_max']])
    g['ax_{}_{}'.format(z,0)].set_ylabel('Minimum cell voltage,\n$U_\mathrm{min}$ in V')
    ylim = plot_preferences._set_lim(res[rn_min]['U_min'][0],res[rn_min]['U_max'][0],0.1)
    g['ax_{}_{}'.format(z,0)].set_ylim(ylim)
    
    
    
    # ax 1 
    if 'ax_{}_{}'.format(z,1) not in g:
        g['ax_{}_{}'.format(z,1)] = g['fig_{}'.format(z)].add_subplot(grid[-1,0]) 
        # g['ax_{}_{}'.format(z,0)].grid()


    # Operation
    set_EP = c['req']['E/P_req']
    set_P = c['cell']['Pdc_DIS_max']
    dp_intersection = SDM_fcn.calc_intersection_curve(g,res,set_EP,print_error=False,prefix_str=prefix_str)
    EPline_x = dp_intersection['valx_P']
    EPline_y = dp_intersection[prefix_str+'P_cell']
    id_min = min([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
    id_max = max([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
    sPmaxset = EPline_x[id_min] if set_P < EPline_x[id_min] else (EPline_x[id_max] if set_P > EPline_x[id_max] else set_P)
    
    g['pl_{}_{}_{}'.format(z,1,0)] = g['ax_{}_{}'.format(z,1)].plot(dp_intersection[prefix_str+'P_cell'],dp_intersection[prefix_str+'I_cell'],color=farbe[3],zorder=100)
    calc_P,sImaxset = SDM_fcn._set_intersection_limits(g,EPline_x,dp_intersection[prefix_str+'I_cell'],z,1,x=sPmaxset)

    SDM_fcn.plot_iso_intersection_curves(g,z,1,c,res,prefix_str=prefix_str,iso_mins=iso_mins,var='Imax',print_ann=False)


    # Labeling
    g['ax_{}_{}'.format(z,1)].set_xlabel('Maximum cell power, $P^\mathrm{dis}_\mathrm{cell,max}$ in W')
    g['ax_{}_{}'.format(z,1)].set_xlim([0,70 if sim_limits else res[rn_min]['lim_P_dis_max']])
    g['ax_{}_{}'.format(z,1)].set_ylabel('Maximum cell current,\n$I^\mathrm{dis}_\mathrm{cell,max}$ in A')
    ylim = plot_preferences._set_lim(0,30 if sim_limits else res[rn_min]['I_dis_max'][0],0.1)
    g['ax_{}_{}'.format(z,1)].set_ylim(ylim)
    
    
    
    # ax 2
    if 'ax_{}_{}'.format(z,2) not in g:
        g['ax_{}_{}'.format(z,2)] = g['fig_{}'.format(z)].add_subplot(grid[:,1:])
        # g['ax_{}_{}'.format(z,2)].grid()
   
    xlim = [150,750]
    ylim = [5,1000]
    SDM_fcn.plot_limit_loci(g,z,2,c,res,xlim=xlim,ylim=ylim,prefix_str=prefix_str,scale_log=True,printapp=True)

    SDM_fcn.plot_iso_intersection_curves(g,z,2,c,res,prefix_str=prefix_str,iso_mins=iso_mins,var='general',print_ann=False)
    # SDM_fcn._plot_iso_CPapp(g['ax_{}_{}'.format(z,2)],c)

   
    # Labeling
    g['ax_{}_{}'.format(z,2)].set_xlabel('Number of cells in series, $x_\mathrm{S}$')
    g['ax_{}_{}'.format(z,2)].set_ylabel('Number of cells in parallel, $y_\mathrm{P}$')

    grid.tight_layout(g['fig_{}'.format(z)])


    #%%
    # Plot: Efficiency Curve of the Power Electric Converter (PEC)
    # -------------------------------------------------------------------------
    z = 5
    grid = plt.GridSpec(1, 1, wspace=0.55, hspace=0.05)
    if 'fig_{}'.format(z) not in g:
        scale_y = 1 #0.75
        g['fig_{}'.format(z)] = plt.figure(z,figsize=[fig_size[-1],fig_size[-1]*scale_y])
        
    if 'ax_{}_{}'.format(z,0) not in g:
        g['ax_{}_{}'.format(z,0)] = g['fig_{}'.format(z)].add_subplot(grid[:,:]) 
   
    for n,key in enumerate(['fcn_eta(Umin)','fcn_eta(Umax)']):
        g['pl_{}_{}_{}'.format(z,0,n)] = g['ax_{}_{}'.format(z,0)].plot(c['pec'][key].x,c['pec'][key].y*100,marker='None',ls='dashed' if n==0 else '-',color=farbe[1] if n==0 else farbe[0],label='$U^{{\mathrm{{DC}}}}_{{\mathrm{{pec,{}}}}}$'.format('min' if n==0 else 'max'))#,key.split('=')[-1][:-1]))
        if n==0:
            id_fUmin = c['pec'][key].y.argmax()
            xmax_Umin = c['pec'][key].x[id_fUmin]
            ymax_Umin = c['pec'][key].y[id_fUmin]
            g['ax_{}_{}'.format(z,0)].plot([xmax_Umin]*3,[-0.1,ymax_Umin*100,100.1],marker='.',color=farbe[1],ls='dotted',lw=1)
        else:
            id_fUmax = c['pec'][key].y.argmax()
            xmax_Umax = c['pec'][key].x[id_fUmax]
            ymax_Umax = c['pec'][key].y[id_fUmax]
            g['ax_{}_{}'.format(z,0)].plot([xmax_Umax]*3,[-0.1,ymax_Umax*100,100.1],marker='.',color=farbe[1],ls='dotted',lw=1)
    
    
    g['ax_{}_{}'.format(z,0)].annotate(r'$\eta_\mathrm{pec,max}$', xy=(xmax_Umin,ymax_Umin*100),textcoords="offset points",xytext=(0,5/scale_y),va='bottom',ha='center',bbox=bbox_args,color=farbe[1])
    
    # update eta_op
    ratio_Pac = c['req']['Pac_req']/(c['pec']['Pac_max']*c['pec']['n'])
    g['pl_{}_{}_{}'.format(z,0,10)] = g['ax_{}_{}'.format(z,0)].plot([ratio_Pac]*2,[-0.1,100.1],color=farbe[2],ls='dotted',lw=1)
    g['pl_{}_{}_{}'.format(z,0,11)] = g['ax_{}_{}'.format(z,0)].plot([ratio_Pac],[c['pec']['fcn_eta(Umax)'](ratio_Pac)*100],color=farbe[2],marker='o',markersize=4.5,ls='dotted',lw=1)
    g['an_{}_{}_{}'.format(z,0,11)] = g['ax_{}_{}'.format(z,0)].annotate(r'$\eta_\mathrm{pec,set}$', xy=(ratio_Pac,c['pec']['fcn_eta(Umax)'](ratio_Pac)*100),textcoords="offset points",xytext=(-5,-5/scale_y),va='top',ha='right',color=farbe[2])
    
    c['pec']['eta_op'] = SDM_fcn.update_eta(g,c,z)
   
   
    # Legende
    g['leg_{}'.format(z)] = {}
    handles, labels = g['ax_{}_{}'.format(z,0)].get_legend_handles_labels()
    g['leg_{}'.format(z)].update(dict(zip(labels, handles)))
            
    # Labeling
    if g['ax_{}_{}'.format(z,0)].get_legend() == None:
        leg = g['ax_{}_{}'.format(z,0)].legend(bbox_to_anchor=(1,0), loc='lower right', fancybox=True) #bbox_to_anchor=(0.025,0), loc='lower left'
        leg.get_title().set_multialignment('center')
   
    g['ax_{}_{}'.format(z,0)].set_xlabel(r'AC-side power ratio, $P^{\mathrm{AC}}_{\mathrm{pec,out}} / P^{\mathrm{AC}}_{\mathrm{pec,rated}}$ in -')
    g['ax_{}_{}'.format(z,0)].set_xlim([0,1])
    g['ax_{}_{}'.format(z,0)].set_ylabel(r'Efficiency, $\eta_{\mathrm{pec}}$ in \%')
    g['ax_{}_{}'.format(z,0)].set_ylim([90,100])
    
    bbox = g['ax_{}_{}'.format(z,0)].get_position()
    g['ax_{}_{}'.format(z,0)].set_position((bbox.x0, bbox.y0+0.05, bbox.width, bbox.height))



    #%%
    # Slider Widget
    # -------------------------------------------------------------------------
    z = 0
    if 'fig_{}'.format(z) not in g:
            g['fig_{}'.format(z)] = plt.figure(z,figsize=[fig_size[-1],fig_size[-1]/2])    
    
    # Slider Settings
    axcolor = 'lightgoldenrodyellow'
    axEPset = g['fig_{}'.format(z)].add_axes([0.15, 0.66, 0.7, 0.05], facecolor=axcolor, xscale='log')
    axEPset.add_artist(axEPset.xaxis)

    axPset = g['fig_{}'.format(z)].add_axes([0.15, 0.33, 0.7, 0.05], facecolor=axcolor)
    axPset.add_artist(axPset.xaxis)
        
    # axUmin = g['fig_{}'.format(z)].add_axes([0.15, 0.2, 0.7, 0.03], facecolor=axcolor)
    # axUmin.add_artist(axUmin.xaxis)
    # axUmax = g['fig_{}'.format(z)].add_axes([0.15, 0.3, 0.7, 0.03], facecolor=axcolor)
    # axUmax.add_artist(axUmax.xaxis)
    # axImax = g['fig_{}'.format(z)].add_axes([0.15, 0.4, 0.7, 0.03], facecolor=axcolor)
    # axImax.add_artist(axImax.xaxis)
    # axPmax = g['fig_{}'.format(z)].add_axes([0.15, 0.5, 0.7, 0.03], facecolor=axcolor)
    # axPmax.add_artist(axPmax.xaxis)    

    sEPset_0 = c['req']['E/P_req']
    sEPset = Slider(axEPset, r'$\mathrm{(E/P)_{req}}$', 0.001, 10, valinit=sEPset_0)#, orientation='vertical')
    sPset_0 = dp_intersection[prefix_str+'P_cell'][~np.isnan(dp_intersection[prefix_str+'P_cell'])][-1]
    sPset = Slider(axPset, r'$P^\mathrm{dis}_\mathrm{cell,max}$', g['ax_{}_{}'.format(2,0)].get_xlim()[0], g['ax_{}_{}'.format(2,0)].get_xlim()[1], valinit=sPset_0, color=farbe[2])#, orientation='vertical')



    # Update Functions for Slider Widget
    def replot_ERP(set_EP,set_P,zfig=2):
        """
        Dynamically updates and replots the Extended Ragone Plot (ERP) based on 
        user-adjusted slider values. This function is integral to the interactive 
        element of the System Design Method (SDM), allowing real-time visualization 
        of the effects of changing the energy-to-power ratio (E/P) and maximum cell 
        power on the ERP.
        
        Parameters:
        - set_EP (float): The energy-to-power ratio value adjusted by the user.
        - set_P (float): The maximum cell power value adjusted by the user.
        - zfig (int, optional): Figure index in the plotting environment. Defaults to 2.
        
        Functionality:
        - Recalculates intersection curves based on the new E/P ratio.
        - Updates the plots with new data points reflecting the adjusted parameters.
        - Resets the maximum cell power slider to align with the new E/P ratio.
        - Adjusts system design parameters based on the new settings and replots 
          relevant figures.
        
        Note:
        - This function is called upon slider value changes and is key for interactive
          analysis in the SDM framework.
        - Requires global variables like 'g', 'res', and 'sPset' to be predefined 
          and accessible.
        """

        tic = time.time()
        print("[Info] Calculation of E/P Intersection Curves")
        print("E/P_req = {:.5f} h".format(set_EP))
        
        dp_intersection = SDM_fcn.calc_intersection_curve(g,res,set_EP,print_error=False,prefix_str=prefix_str)
        EPline_x = dp_intersection['valx_P']
        EPline_y = dp_intersection[prefix_str+'P_cell']
        id_min = min([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
        id_max = max([n for n,i in enumerate(EPline_y) if ~np.isnan(i)])
        sPmaxset = EPline_x[id_min] if set_P < EPline_x[id_min] else (EPline_x[id_max] if set_P > EPline_x[id_max] else set_P)
        
        try:
            g['pl_{}_{}_{}'.format(zfig,0,100)][0].set_data(EPline_x,dp_intersection['valy_E'])
            _,sEmaxset = SDM_fcn._set_intersection_limits(g,EPline_x,dp_intersection['valy_E'],zfig,0,x=sPmaxset,replot=True)
    
            g['pl_{}_{}_{}'.format(zfig,1,100)][0].set_data(EPline_x,dp_intersection[prefix_str+'U_cell'])
            id_max = max([n for n,i in enumerate(dp_intersection[prefix_str+'U_cell']) if ~np.isnan(i)])
            g['pl_{}_{}_{}'.format(zfig,1,101)][0].set_data(EPline_x[id_max],dp_intersection[prefix_str+'U_cell'][id_max])
            _,sUminset = SDM_fcn._set_intersection_limits(g,EPline_x,dp_intersection[prefix_str+'U_cell'],zfig,1,x=sPmaxset,replot=True)
                            
            g['pl_{}_{}_{}'.format(zfig,2,100)][0].set_data(EPline_x,dp_intersection[prefix_str+'I_cell'])
            id_max = max([n for n,i in enumerate(dp_intersection[prefix_str+'I_cell']) if ~np.isnan(i)])
            g['pl_{}_{}_{}'.format(zfig,2,101)][0].set_data(EPline_x[id_max],dp_intersection[prefix_str+'I_cell'][id_max])
            _,sImaxset = SDM_fcn._set_intersection_limits(g,EPline_x,dp_intersection[prefix_str+'I_cell'],zfig,2,x=sPmaxset,replot=True)
         
            g['fig_{}'.format(zfig)].canvas.draw()
        except:
            0
        
        # reset sPset when changing sEPset
        sPset.eventson = False
        sPset.valmin = EPline_x[id_min]
        sPset.valmax = EPline_x[id_max]
        sPset.vline.set_data([sPset.valmax]*2,sPset.vline.get_data()[-1])
        # sPset.hline.set_data(sPset.hline.get_data()[0],[sPset.valmax]*2)
        sPset.set_val(sPmaxset)
        sPset.eventson = True
        
        
        try:
            c['cell']['Udc_min'] = sUminset
            c['cell']['Idc_DIS_max'] = sImaxset
            c['cell']['Pdc_DIS_max'] = sPmaxset
            update_dss(0)
        except:
            0
        
        print('[Info] tReplot = {:.3f}s'.format(time.time()-tic))


    def update_intersection(val):
        """
        Updates the energy-to-power ratio (E/P) and maximum cell power parameters based on 
        user input from the slider widget in the System Design Method (SDM) tool. It 
        triggers the replotting of the Extended Ragone Plot (ERP) to reflect these changes.
        
        Parameters:
        - val: The new value from the slider input.
        
        Functionality:
        - Adjusts system parameters for E/P ratio and maximum cell power in response to 
          slider changes.
        - Calls 'replot_ERP' to update the ERP visualization with the new parameters.
        
        Usage:
        - Integral to the interactive element of the SDM tool, allowing users to explore 
          the impact of parameter variations on system design and performance.
        
        Note:
        - This function is part of the interactive analysis features in the SDM tool.
        """
        c['req']['E/P_req'] = sEPset.val
        set_P = sPset.val
        replot_ERP(c['req']['E/P_req'],set_P)
         
    
    def update_dss(val):
        """
        Updates the Design Solution Space (DSS) visualization in response to parameter 
        changes in the System Design Method (SDM) analysis tool. It refreshes the DSS plot 
        and constraint intersection trajectory based on the current parameter values.
        
        Parameters:
        - val: The updated value triggering the function, typically from a UI element like 
          a slider or input field.
        
        Functionality:
        - Calls the 'plot_dss' function to redraw the DSS with updated parameters.
        - Invokes 'plot_constraint_intersection_trajectory' to update the visualization 
          of constraints based on the new parameter values.
        - Utilizes matplotlib's 'draw' function to render the updated plots on the canvas.
        
        Usage:
        - Essential for maintaining an interactive and responsive visualization environment 
          in the SDM tool, allowing users to see the immediate impact of parameter changes 
          on the DSS and constraints.
        
        Note:
        - This function plays a key role in the interactivity of the SDM tool, ensuring 
          dynamic visualization updates.
        """
        # c['pec']['Udc_min'] = sUinvmin.val
        # c['pec']['Udc_max'] = sUinvmax.val
        # c['pec']['Idc_max'] = sIinvmax.val
        # c['pec']['Pdc_max'] = sPinvmax.val
        
        SDM_fcn.plot_dss(g,c,zfig=3,zax=0,print_opt=True)
        SDM_fcn.plot_constraint_intersection_trajectory(g, res, c, prefix_str=prefix_str)

        plt.draw()
       
    
    def update_dp(req_EP,pec_Umin,pec_Umax,pec_Imax,pec_Pmax,cell_P=0):
        """
        Adjusts the values of slider widgets based on specified electrical parameters of 
        power electronic converters (PEC) and energy storage cells. This function aligns 
        slider positions with changes in PEC and cell parameters, enabling dynamic and 
        interactive analysis in the SDM tool.
        
        Parameters:
        - req_EP (float): The required energy-to-power ratio (E/P) setting for the system.
        - pec_Umin (float): The minimum voltage limit for the power electronic converter.
        - pec_Umax (float): The maximum voltage limit for the power electronic converter.
        - pec_Imax (float): The maximum current limit for the power electronic converter.
        - pec_Pmax (float): The maximum power limit for the power electronic converter.
        - cell_P (float, optional): The cell power value. If greater than 0, the slider 
          is set to this value; otherwise, the slider is reset. Defaults to 0.
        
        Functionality:
        - Sets the energy-to-power ratio slider to the 'req_EP' value.
        - Adjusts the cell power slider based on the 'cell_P' parameter.
        - (Optional) Adjusts additional sliders for PEC limits if uncommented.
        
        Usage:
        - Integral to the interactive aspect of the SDM tool, allowing for the adjustment 
          of system parameters through a graphical interface.
        
        Note:
        - This function is essential for ensuring that changes in PEC and cell parameters 
          are accurately reflected in the corresponding sliders, enhancing the interactive 
          user experience in system design exploration.
        """

        # Manipulate Sliders
        sEPset.set_val(req_EP)

        if cell_P>0:
            sPset.set_val(cell_P)
        else:
            sPset.reset()
            sPset.set_val(sPset.valmax)

        # sUinvmin.eventson = False
        # sUinvmin.set_val(pec_Umin)
        # sUinvmin.eventson = True

        # sUinvmax.eventson = False
        # sUinvmax.set_val(pec_Umax)
        # sUinvmax.eventson = True

        # sIinvmax.eventson = False
        # sIinvmax.set_val(pec_Imax)
        # sIinvmax.eventson = True
        
        # sPinvmax.set_val(pec_Pmax)


    # Update Slider Widget
    sEPset.on_changed(update_intersection)
    sPset.on_changed(update_intersection)
    
    # sUinvmin = Slider(axUmin, r'$U_\mathrm{inv,min}$', 0, 1500, valinit=c['pec']['Udc_min'])
    # sUinvmin_xlim = sUinvmin.ax.get_xlim()
    # sUinvmin_xy = sUinvmin.poly.get_xy()
    # sUinvmin_xy[np.where(sUinvmin_xy == sUinvmin_xlim[0])] = sUinvmin_xlim[-1]
    # sUinvmin.poly.set_xy(sUinvmin_xy)
    
    # sUinvmax = Slider(axUmax, r'$U_\mathrm{inv,max}$', 0, 1500, valinit=c['pec']['Udc_max'])
    # sIinvmax = Slider(axImax, r'$I_\mathrm{inv,max}$', 0, 2000, valinit=c['pec']['Idc_max'])
    # sPinvmax = Slider(axPmax, r'$P_\mathrm{inv,max}$', 0, 1500*2000, valinit=c['pec']['Pdc_max'])
    
    # sUinvmin.on_changed(update_dss)
    # sUinvmax.on_changed(update_dss)
    # sIinvmax.on_changed(update_dss)
    # sPinvmax.on_changed(update_dss)
    
    
    # Initial Update of Slider Positions
    sEPset.set_val(sEPset.val)


    #%%
    print('##################################################################')
    print('SYSTEM DESIGN METHOD: USER SETTINGS')

    set_case(1)