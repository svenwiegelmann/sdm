"""
System Design Method Interactive Tool (SDM_main.py)

Author:     Sven Wiegelmann,
            Leibniz University Hannover,
            Institute of Electric Power Systems,
            Electric Energy Storage Systems Section,
            Appelstraße 9A,
            30167 Hannover,
            Germany
Version:    18.12.2024

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
from matplotlib.widgets import Slider# TextBox, Button, RadioButtons, RangeSlider

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

import plot_preferences as plot_preferences
import SDM_fcn as SDM_fcn
import SDM_data as SDM_data


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
    global res
    load_res = SDM_data.load_results(0.0,True)
    
    if isinstance(load_res,list):
        load_res_keys = []
        load_res = SDM_data.adjust_res_keys(load_res,prefix_str=prefix_str)
        res = load_res
    else:
        load_res_keys = list(load_res.keys())
        for n,key in enumerate(load_res_keys):
            load_res[key] = SDM_data.adjust_res_keys(load_res[key],prefix_str=prefix_str)
            if n>0: # sort initial DIS processes of 'indU'-loops
                load_res[load_res_keys[n-1]][0] = load_res[key][0]
    
        key_res = 'indU=0'
        res = load_res[key_res]
    
    if isinstance(load_res,dict) and 'ocv' in load_res:
        load_ocv = load_res['ocv']
        load_res.pop('ocv',None)
        load_res_keys = list(load_res.keys())
    
    tmp_res = res
    
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
            c['pec'] = SDM_data.set_pec(0)
            c['pec']['n'] = 2
            c['pec']['Udc_min'] = 270
            c['pec']['Udc_max'] = 500
            c['pec']['Idc_max'] = 1675
            c['pec']['Pdc_max'] = 550000
            
            c['mod']['xS'] = 4
            c['mod']['yP'] = 2
            
            set_dp(set_EPreq=c['req']['E/P_req'],
                   set_Pmax=100)
            
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
            
            set_dp(set_EPreq=c['req']['E/P_req'])
            
        if no in [2,3,4]: # no DSS available --> Limit Variation
            c['req'] = SDM_data.set_app(2)
            
            c['pec'] = SDM_data._set_by_arg(SDM_data.set_pec,'name','SMA Solar SCS800')
            c['pec']['n'] = 1 #np.ceil(c['req']['Pac_req']/c['pec']['Pac_max'])
        
            c['mod']['xS'] = 2
            c['mod']['yP'] = 2
            
            set_dp(set_EPreq=c['req']['E/P_req'],
                   set_Pmax=0 if (no==2) else 138 if (no==3) else 133)

        if (no==5): # no DSS ever available
            c['req'] = SDM_data.set_app(3)
            
            c['pec'] = SDM_data._set_by_arg(SDM_data.set_pec,'name','SMA Solar SCS800')
            c['pec']['n'] = 1 #np.ceil(c['req']['Pac_req']/c['pec']['Pac_max'])
        
            c['mod']['xS'] = 1
            c['mod']['yP'] = 1
            
            set_dp(set_EPreq=c['req']['E/P_req'])
            
        if no in [6,7]: # Methodology Comparison
            c['req'] = SDM_data.set_app(3)
            
            c['pec'] = SDM_data._set_by_arg(SDM_data.set_pec,'name','SMA Solar SCS800')
            c['pec']['n'] = 1 #np.ceil(c['req']['Pac_req']/c['pec']['Pac_max'])
        
            c['mod']['xS'] = 2
            c['mod']['yP'] = 1
            
            Ecell_nom = 78.4
            Pcell_nom = 2.21 * 198
            
            c['req']['E/P_req'] = Ecell_nom/Pcell_nom
                        
            set_dp(set_EPreq=c['req']['E/P_req'])
    
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
    
    set_figsize = {'h1':[fig_size[-1]]*2,
                   'h2':[fig_size[-1]*2.05,fig_size[-1]],
                   'h3':[fig_size[-1]*2.95,fig_size[-1]*1.05]}

    alph = 1/3


    #%%
    # Plot: Charge/Discharge Curves U(E), I(E), Derivation of the Extended Ragone Plot (ERP)
    # -------------------------------------------------------------------------  
    nPref = 15

    pub_lim = {'Umax':      {'n':    1,
                             'val':  res[rn_min]['U_max'][0],
                             'type': 'upper',
                             'c':    farbe[10]},
               'Umin':      {'n':    2,
                             'val':  res[rn_min]['U_min'][0],
                             'type': 'lower',
                             'c':    farbe[5]},
               'Cratemax':  {'n':    3,
                             'val':  res[rn_min]['Crate_max'][0],
                             'type': 'upper',
                             'c':    farbe[3]},
               'Cratemin':  {'n':    4,
                             'val':  np.nan,
                             'type': 'lower',
                             'c':    'y'},
               'Tmax':      {'n':    5,
                             'val':  res[rn_min]['T_max'][0],
                             'type': 'upper',
                             'c':    farbe[2]},
               'Tmin':      {'n':    6,
                             'val':  np.nan,
                             'type': 'lower',
                             'c':    'y'},
               'SoCmax':    {'n':    7,
                             'val':  1.2,
                             'type': 'upper',
                             'c':    farbe[4]},
               'SoCmin':    {'n':    8,
                             'val':  -0.2,
                             'type': 'lower',
                             'c':    farbe[8]}}
    for key in pub_lim:
        pub_lim[key]['cut'] = pub_lim[key]['val']

    ## Calculate specific Operating Limits
    pub_lim['Umin']['cut'] = 1.95
    pub_lim['Cratemax']['cut'] = 4.25

    # # # pub_lim['Umax']['cut'] = 2.7
    # pub_lim['Umax'].update({'ocv': {'calc_SoC': True,
    #                                 'SoC':  SoC_chg,
    #                                 'pOCV': pOCV_chg}})


    z = 1
    grid = plt.GridSpec(1,3, wspace=0.55, hspace=0.05)
    if 'fig_{}'.format(z) not in g:
        g['fig_{}'.format(z)] = plt.figure(z,figsize=set_figsize['h3'])


    # ax 0
    if 'ax_{}_{}'.format(z,0) not in g:
        g['ax_{}_{}'.format(z,0)] = g['fig_{}'.format(z)].add_subplot(grid[:,2]) 
        # g['ax_{}_{}'.format(z,0)].grid()
    
    nSoC = 0
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
                                      # textcoords="offset points",
                                      # xytext=(0,0),
                                      ha='center',
                                      va='center',
                                      c=farbe[12],
                                      bbox=bbox_args,
                                      zorder=100)
    
    ### Indices: 20+ single limits; 40+ bidirectional limit interrelations; 99 overall limit interrelations    
    for key in pub_lim:
        if key in ['Cratemin','Tmin']:
            continue
        res_recalc = SDM_fcn.recalc_limits(res,
                                            Umax=pub_lim['Umax']['cut'] if key=='Umax' else pub_lim['Umax']['val'],
                                            Umin=pub_lim['Umin']['cut'] if key=='Umin' else pub_lim['Umin']['val'],
                                            Cratemax=pub_lim['Cratemax']['cut'] if key=='Cratemax' else pub_lim['Cratemax']['val'],
                                            Tmax=pub_lim['Tmax']['cut'] if key=='Tmax' else pub_lim['Tmax']['val'],
                                            SoCmax=pub_lim['SoCmax']['cut'] if key=='SoCmax' else pub_lim['SoCmax']['val'],
                                            SoCmin=pub_lim['SoCmin']['cut'] if key=='SoCmin' else pub_lim['SoCmin']['val'],
                                            prefix_str=prefix_str)
        d_recalc, tmp_ind_recalc, rn_min_recalc = SDM_fcn.compile_EP(res_recalc,calc_EPsys=calc_EPsys,prefix_str=prefix_str)
        g['pl_{}_{}_{}_{}'.format(z,0,20+pub_lim[key]['n'],nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d_recalc['p'],d_recalc[nSoC][0],color=pub_lim[key]['c'],marker='.',markersize=4,lw=1)
    

    # overall   
    res_recalc = SDM_fcn.recalc_limits(res,
                                        Umax=pub_lim['Umax']['cut'],
                                        Umin=pub_lim['Umin']['cut'],
                                        Cratemax=pub_lim['Cratemax']['cut'] ,
                                        Tmax=pub_lim['Tmax']['cut'],
                                        SoCmax=pub_lim['SoCmax']['cut'],
                                        SoCmin=pub_lim['SoCmin']['cut'],
                                        prefix_str=prefix_str)
    d_recalc, tmp_ind_recalc, rn_min_recalc = SDM_fcn.compile_EP(res_recalc,calc_EPsys=calc_EPsys,prefix_str=prefix_str)
    g['pl_{}_{}_{}_{}'.format(z,0,99,nSoC)] = g['ax_{}_{}'.format(z,0)].plot(d_recalc['p'],d_recalc[nSoC][0],color=farbe[0],marker='.',markersize=6,lw=1.5)


   
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
                if np.mod(nP,1):# "thinning"
                    continue
                alph = 1/3
                g['pl_{}_{}_{}_{}'.format(z,1,0,rn)] = g['ax_{}_{}'.format(z,1)].plot(res[rn][prefix_str+'E_cell'],res[rn][prefix_str+'U_cell'],color=farbe[nSoC],lw=0.75,alpha=alph)
                g['pl_{}_{}_{}_{}'.format(z,1,1,rn)] = g['ax_{}_{}'.format(z,1)].plot(res_recalc[rn][prefix_str+'E_cell'],res_recalc[rn][prefix_str+'U_cell'],color=farbe[nSoC],lw=1.5)
                if nP == nPref:
                    g['pl_{}_{}_{}_{}'.format(z,1,1,rn)][0].remove()
                    g['pl_{}_{}_{}_{}'.format(z,1,1,rn)] = g['ax_{}_{}'.format(z,1)].plot(res_recalc[rn][prefix_str+'E_cell'],res_recalc[rn][prefix_str+'U_cell'],color=farbe[2],lw=1.5,zorder=100)
            except:
                0
    
    xlim = [0,10 if sim_limits else res[rn_min]['lim_E_dis_max']]
    
    g['ax_{}_{}'.format(z,1)].plot(xlim,[pub_lim['Umin']['val']]*2,linestyle='dashed',color=[0]*3,lw=1,zorder=1)
    g['ax_{}_{}'.format(z,1)].plot(xlim,[pub_lim['Umax']['val']]*2,linestyle='dashed',color=[0]*3,lw=1,zorder=1)

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
                                        zorder=99)
    # bbox_args = dict(boxstyle="round",color='None',fc='1',lw=0)
    g['an_{}_{}_{}'.format(z,1,1)] = g['ax_{}_{}'.format(z,1)].annotate(r'$U_\mathrm{cell,max}$',
                                        (0.1*xlim[-1],res[rn_min]['U_max'][0]),
                                        # textcoords="offset points",
                                        # xytext=(0,0),
                                        ha='center',
                                        va='center',
                                        # fontsize='large',
                                        bbox=bbox_args,
                                        arrowprops = dict(arrowstyle="<|-",connectionstyle="angle,angleA=90,angleB=0,rad=0",color=farbe[4]),
                                        zorder=99)
    
    g['pl_{}_{}_{}'.format(z,1,0)] = g['ax_{}_{}'.format(z,1)].plot(xlim,[pub_lim['Umin']['val']]*2,color=pub_lim['Umin']['c'],lw=1,ls='dashed',zorder=98)
    g['pl_{}_{}_{}'.format(z,1,0)][0].set_visible(0)
    g['pl_{}_{}_{}'.format(z,1,1)] = g['ax_{}_{}'.format(z,1)].plot(xlim,[pub_lim['Umax']['val']]*2,color=pub_lim['Umax']['c'],lw=1,ls='dashed',zorder=98)
    g['pl_{}_{}_{}'.format(z,1,1)][0].set_visible(0)
            
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
                if np.mod(nP,1):# "thinning"
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
                                      zorder=99)
    # set zorder of 2Dline above annotation
    # g['pl_{}_{}_{}_{}'.format(z,2,1,36)][0].set_zorder(g['an_{}_{}_{}'.format(z,2,0)].get_zorder()+1)

    
    # g['ax_{}_{}'.format(z,2)].plot(xlim,[res[rn_min]['Crate_max'][0]]*2,linestyle='dashed',color=[0]*3,lw=1.5,zorder=1)
    g['pl_{}_{}_{}'.format(z,2,0)] = g['ax_{}_{}'.format(z,2)].plot(xlim,[pub_lim['Cratemax']['val']*res[rn_min]['_settings']['ParameterValues']['Q_n'][0]]*2,color=pub_lim['Cratemax']['c'],ls='dashed',lw=1,zorder=98)
    g['pl_{}_{}_{}'.format(z,2,0)][0].set_visible(0)

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
    
    xlim = g['ax_{}_{}'.format(z,1)].get_xlim()
    ylim = g['ax_{}_{}'.format(z,1)].get_ylim()
    g['an_{}_{}_{}'.format(z,1,100)] = g['ax_{}_{}'.format(z,1)].annotate(r'$P^\mathrm{dis}_\mathrm{cell} \uparrow$',
                                      (xlim[0]+0.71*(xlim[1]-xlim[0]),ylim[0]+0.52*(ylim[1]-ylim[0])),
                                      # textcoords="offset points",
                                      xytext=(xlim[0]+0.35*(xlim[1]-xlim[0]),ylim[0]+0.32*(ylim[1]-ylim[0])),
                                      ha='center',
                                      va='center',
                                      color=[0.5]*3,
                                      # fontsize='large',
                                      bbox=bbox_args,
                                      arrowprops = dict(arrowstyle="<|-",color=[0.5]*3,connectionstyle='angle3,angleA=5,angleB=75'),
                                      zorder=0)
    
    xlim = g['ax_{}_{}'.format(z,2)].get_xlim()
    ylim = g['ax_{}_{}'.format(z,2)].get_ylim()
    g['an_{}_{}_{}'.format(z,2,100)] = g['ax_{}_{}'.format(z,2)].annotate(r'$P^\mathrm{dis}_\mathrm{cell} \uparrow$',
                                      (xlim[0]+0.94*(xlim[1]-xlim[0]),ylim[0]+0.41*(ylim[1]-ylim[0])),
                                      # textcoords="offset points",
                                      xytext=(xlim[0]+0.9*(xlim[1]-xlim[0]),ylim[0]+0.78*(ylim[1]-ylim[0])),
                                      ha='center',
                                      va='center',
                                      color=[0.5]*3,
                                      # fontsize='large',
                                      bbox=bbox_args,
                                      arrowprops = dict(arrowstyle="<|-",color=[0.5]*3,connectionstyle='angle3,angleA=90,angleB=100'),
                                      zorder=0)

    grid.tight_layout(g['fig_{}'.format(z)])
    # g['fig_{}'.format(z)].canvas.draw()


    SDM_fcn.replot_EP_limits(g,res,tmp_ind,pub_lim=pub_lim,zfig=1,nPref=nPref,calc_EPsys=calc_EPsys,print_ann=False,prefix_str=prefix_str)



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
   
    SDM_fcn.plot_EP_limit(g, res, z, 0, rn_min, prefix_str=prefix_str, plot_limit=['Umin','Cratemax'])
   
    # Labeling
    g['ax_{}_{}'.format(z,0)].set_xlabel('Cell power, $P^\mathrm{{dis}}_\mathrm{cell}$ in W')
    g['ax_{}_{}'.format(z,0)].set_xlim([0,70 if sim_limits else res[rn_min]['lim_P_dis_max']])
    g['ax_{}_{}'.format(z,0)].set_ylabel('Cell energy, $E^\mathrm{{dis}}_\mathrm{{cell}}$ in {}'.format('Wh' if convert_unit_E else 'Ws'))
    g['ax_{}_{}'.format(z,0)].set_ylim([0,10 if sim_limits else res[rn_min]['lim_E_dis_max']])
    g['ax_{}_{}'.format(z,0)].set_aspect(1./g['ax_{}_{}'.format(z,0)].get_data_ratio())
    SDM_fcn._plot_isochrones(g['ax_{}_{}'.format(z,0)],num_subplot=True)
   


    # Plot: Intersection Curves of Operational Quantities (by Limit Variation) 
    # -------------------------------------------------------------------------        
    dp_intersection = SDM_fcn.calc_intersection_curve(g,res,c['req']['E/P_req'],print_error=False,prefix_str=prefix_str)
    g['pl_{}_{}_{}'.format(z,0,100)] = g['ax_{}_{}'.format(z,0)].plot(dp_intersection['valx_P'],dp_intersection['valy_E'],color=farbe[2],zorder=100)
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
    
    
    # ax2
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
    # Slider Widget Functions
    # -------------------------------------------------------------------------
    def update_plots(slider_name, val):
        """
        Handles updates to plots and system parameters in response to slider changes.
    
        Parameters:
        - slider_name (str): Name of the slider triggering the update.
        - val (float): The new value of the slider.
    
        Returns:
        - None:  Updates the global `state` dictionary, recalculates system parameters, 
          and refreshes plots across multiple figures.
    
        Note:
        - Updates the `state` dictionary to reflect new slider values.
        - Adjusts the energy-to-power ratio, maximum power, and cell-specific parameters.
        - Handles plot updates across multiple figures:
            - Fig. 1: ERP visualization (`U(E)`, `I(E)`, and ERP limits).
            - Fig. 2: ERP limit intersections (`Imax(Pmax)` and `Umin(Pmax)`).
            - Fig. 3: Constraint Satisfaction Problem (CSP) and Decision Support System (DSS) plots.
            - Fig. 4: Graphical synthesis including iso-intersection curves and trajectories.
        - Removes deprecated plot elements (e.g., ERP limits, intersection trajectories) 
        - Integrates feedback loops between slider values, system constraints, and plot elements.
          If not required, set bool_loc_ax_limits=True!
    
        Dependencies:
        - Relies on `plot_preferences` for plot management (e.g., removing outdated objects).
        - Uses `SDM_fcn` for recalculating and replotting ERP, CSP, and limit intersections.
        - Requires access to a global `state` dictionary for storing intermediate 
          values and slider configurations.
        """
        # Update slider value
        state["slider_values"][slider_name] = val
        
        # Update state['c']
        state['c']['req']['E/P_req'] = state['slider_values']['sEPset']
        state['c']['req']['Eac_req'] = state['c']['req']['E/P_req']*state['c']['req']['Pac_req']
        state['c']['cell']['Pdc_DIS_max'] = state['slider_values']['sPset']
        
        # === Fig. 2 | ERP, Limit Intersections: Imax(Pmax),Umin(Pmax) ===        
        SDM_fcn.replot_ERP(g,state['res'],state['c'],sPset,set_EP=state["slider_values"]["sEPset"],set_P=state["slider_values"]["sPset"],prefix_str=prefix_str)
        
        # === Fig. 3 | CSP/DSS ===
        SDM_fcn.plot_dss(g,state['c'],print_opt=True)
        
        # === Fig. 4 | Graphical Synthesis ===  
        SDM_fcn.plot_constraint_intersection_trajectory(g,state['res'],state['c'],prefix_str=prefix_str)
        
        # === Fig. 1 | U(E),I(E),ERP ===
        # Update with Feedback: c -> pub_lim (after SDM_fcn.replot_ERP!)
        pub_lim['Umin']['cut'] = c['cell']['Udc_min']
        pub_lim['Cratemax']['cut'] = c['cell']['Idc_DIS_max']/state['res'][0]['_settings']['ParameterValues']['Q_n'][0]
        
        # Coupling between different Limits forced by bool_loc_ax_limits=False
        SDM_fcn.replot_EP_limits(g,res,tmp_ind,pub_lim=pub_lim,zfig=1,nPref=nPref,calc_EPsys=calc_EPsys,print_ann=False,prefix_str=prefix_str,bool_loc_ax_limits=False)
    
    
    def set_dp(set_EPreq=0.2,set_Pmax=0):
        """
        Configures the actual design point by manipulating the energy-to-power
        ratio (EP) and maximum power (Pmax) sliders.
    
        Parameters:
        - set_EPreq (float, optional): Desired energy-to-power ratio for the design 
          point. Defaults to `0.2`.
        - set_Pmax (float, optional): Maximum power for the design point. If set to 
          `0`, the slider is reset to its maximum value. Defaults to `0`.
    
        Returns:
        - None: This function updates the slider values but does not return any value.
    
        Note:
        - Directly manipulates `sEPset` and `sPset` slider objects to set the design 
          point parameters.
        - Resets the `Pmax` slider if the input value is `0` and then sets it to its 
          maximum allowable value.
        """
        sEPset.set_val(set_EPreq)

        if set_Pmax>0:
            sPset.set_val(set_Pmax)
        else:
            sPset.reset()
            sPset.set_val(sPset.valmax)

    
    #%%
    # Plot: Slider Widget
    # -------------------------------------------------------------------------
    z = 0
    if 'fig_{}'.format(z) not in g:
            g['fig_{}'.format(z)] = plt.figure(z,figsize=[fig_size[-1],fig_size[-1]/2])    
    
    # Slider States
    state = {'res':  res,
              'c':    c,
              'slider_values': {# recent slider values
                                'sEPset':       c['req']['E/P_req'],
                                'sPset':        dp_intersection[prefix_str+'P_cell'][~np.isnan(dp_intersection[prefix_str+'P_cell'])][-1],
                                'sUcellmax':    c['cell']['Udc_max']}#pub_lim['Umax']['cut'],}
              }
    
    # Slider Settings
    axcolor = 'lightgoldenrodyellow'
    axEPset = g['fig_{}'.format(z)].add_axes([0.15, 0.66, 0.7, 0.05], facecolor=axcolor, xscale='log')
    axEPset.add_artist(axEPset.xaxis)

    axPset = g['fig_{}'.format(z)].add_axes([0.15, 0.33, 0.7, 0.05], facecolor=axcolor)
    axPset.add_artist(axPset.xaxis)
    
    sEPset = Slider(axEPset, r'$\mathrm{(E/P)_{req}}$', 0.001, 10, valinit=state['slider_values']['sEPset'])#, orientation='vertical')
    sPset = Slider(axPset, r'$P^\mathrm{dis}_\mathrm{cell,max}$', g['ax_{}_{}'.format(2,0)].get_xlim()[0], g['ax_{}_{}'.format(2,0)].get_xlim()[1], valinit=state['slider_values']['sPset'], color=farbe[2])#, orientation='vertical')
    
    # Update Slider Widget
    sEPset.on_changed(lambda val: update_plots("sEPset", val))
    sPset.on_changed(lambda val: update_plots("sPset", val))


    #%%
    print('##################################################################')
    print('SYSTEM DESIGN METHOD: USER SETTINGS')
    set_case(1)