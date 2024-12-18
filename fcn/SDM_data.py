"""
System Design Method Data Processing (SDM_data.py)

Author:     Sven Wiegelmann,
            Leibniz University Hannover,
            Institute of Electric Power Systems,
            Electric Energy Storage Systems Section,
            Appelstraße 9A,
            30167 Hannover,
            Germany
Version:    18.12.2024

Overview:
This module, part of the System Design Method (SDM) toolkit, offers essential 
functions for the loading, processing, and configuring of data in energy storage 
systems. It supports analysis and optimization within the SDM framework, 
featuring functionalities for loading experimental data, setting application and 
power electronics converter parameters, and auxiliary data manipulation and 
retrieval functions.

Primary Functions:
- load_results: Processes experimental data for energy storage systems.
- set_app: Configures application settings for energy scenarios.
- set_pec: Sets parameters for SMA Solar Power Electronics Converters.

Secondary Functions:
- _fcn_eta: Calculates efficiency curves for Power Converter models.
- _set_by_arg: Searches for configurations based on key-value pair matching.

Usage:
Designed for researchers and practitioners in energy engineering, this module 
facilitates integrating experimental data and predefined settings into the 
System Design Method. It is instrumental for the design and analysis of 
cell-based energy systems.

Note:
Users should ensure compatibility with data formats and structures. This module 
forms a part of a comprehensive toolkit for advanced energy system design and analysis.
"""

import pickle
import numpy as np
from scipy import interpolate

import xls2dict_digidata as xls2dict_digidata

#%%
def load_results(no,print_bat_list=False):
    """
    Loads and processes experimental data for specific energy storage system 
    experiments. This function retrieves data, derived parameters and plot
    limits for energy system analysis.
    
    Parameters:
    - no (int): Identifier for the specific dataset to be loaded.
    
    Returns:
    - list: A list of dictionaries with processed data from each experiment. 
      Includes experimental data, parameters, and plot limits for the 
      System Design Method.
    
    Note:
    - Function tailored to a specific data structure; changes may be needed for 
      different formats or sources.
    - Part of a toolkit for designing and analyzing energy storage systems.
    """
    
    # List of Battery Cells
    batteries = {0:    'LTO20130205-03',
                 1:    'LTO20130205-05',
                 2:    'MJ1-1119-01',
                 3:    'MJ1-1119-02',
                 4:    'SCiB-23Ah-01',
                 5:    'SCiB-23Ah-02',
                 6:    'HE2-1123-01',
                 7:    'HE2-1123-02',
                 8:    'M1B-1223-01',
                 9:    'M1B-1223-02',
                 10:   'SLPB8644143-45',
                 11:   'SLPB8644143-353',
                 12:    'HE2-1123-03',
                 13:    'HE2-1123-04',
                 # 14:   'LF50-0423-1',
                 # 15:   'LF50-0423-2',
                 }
    
    if print_bat_list:
        print('List of Battery Cells:')
        for key in batteries:
            print(key,'-',batteries[key])
    
    res = []
    no = float(no)
    
    ## dis [0 - 99]
    if (no==0.0):
        fn = 'Scienlab__LTO20130205__LTO20130205-03__PowerTemperatureTest_dis#3'
    if (no==1.0):
        fn = 'Scienlab__LTO20130205__LTO20130205-05__PowerTemperatureTest_dis#3'
    if (no==2.0):
        fn = 'Scienlab__18650 MJ1__MJ1-1119-01__PowerTemperatureTest_dis#2'
    if (no==3.0):
        fn = 'Scienlab__18650 MJ1__MJ1-1119-02__PowerTemperatureTest_dis#2'
    if (no==4.0):
        fn = 'Scienlab__SCiB 23Ah__SCiB-23Ah-01__PowerTemperatureTest_dis#1'
    if (no==5.0):
        fn = 'Scienlab__SCiB 23Ah__SCiB-23Ah-02__PowerTemperatureTest_dis#1'
    if (no==6.0):
        fn = 'Scienlab__ICR18650 HE2__HE2-1123-01__PowerTemperatureTest_dis#1'
    if (no==7.0):
        fn = 'Scienlab__ICR18650 HE2__HE2-1123-02__PowerTemperatureTest_dis#1'
    if (no==10.0):
        fn = 'Scienlab__SLPB8644143__SLPB8644143-45__PowerTemperatureTest_dis#1'
    if (no==10.1):
        fn = 'Scienlab__SLPB8644143__SLPB8644143-45__PowerTemperatureTest_dis#2'    
    if (no==11.0):
        fn = 'Scienlab__SLPB8644143__SLPB8644143-353__PowerTemperatureTest_dis#1'
    if (no==11.1):
        fn = 'Scienlab__SLPB8644143__SLPB8644143-353__PowerTemperatureTest_dis#2'
    # if (no==14.0):
    #     fn = 'Scienlab__LF50K-73103__LF50-0423-1__PowerTemperatureTest_dis#1'
    # if (no==15.0):
    #     fn = 'Scienlab__LF50K-73103__LF50-0423-2__PowerTemperatureTest_dis#1'

    ## chg [100 - 199]
    if (no==100.0):
        fn = 'Scienlab__LTO20130205__LTO20130205-03__PowerTemperatureTest_chg#1'
    if (no==101.0):
        fn = 'Scienlab__LTO20130205__LTO20130205-05__PowerTemperatureTest_chg#1'
    if (no==104.0):
        fn = 'Scienlab__SCiB 23Ah__SCiB-23Ah-01__PowerTemperatureTest_chg#1'
    if (no==105.0):
        fn = 'Scienlab__SCiB 23Ah__SCiB-23Ah-02__PowerTemperatureTest_chg#1'
        
    ## Variation of Umax [200 - 299]
    if (no==200.0): # 25°C
        fn = 'Scienlab__LTO20130205__LTO20130205-03__PowerTemperatureTest_dis-Umax-manual#1'
    if (no==200.1): # -10°C
        fn = 'Scienlab__LTO20130205__LTO20130205-03__PowerTemperatureTest_dis-Umax-manual#2__Tamb=-10°C'
    if (no==200.2): # 10°C
        fn = 'Scienlab__LTO20130205__LTO20130205-03__PowerTemperatureTest_dis-Umax-manual#3__Tamb=10°C'
    if (no==200.3): # 40°C
        fn = 'Scienlab__LTO20130205__LTO20130205-03__PowerTemperatureTest_dis-Umax-manual#3__Tamb=40°C'
    
    if (no==201.0): # 25°C
        fn = 'Scienlab__LTO20130205__LTO20130205-05__PowerTemperatureTest_dis-Umax-manual#1'
    if (no==201.1): # -10°C
        fn = 'Scienlab__LTO20130205__LTO20130205-05__PowerTemperatureTest_dis-Umax-manual#2__Tamb=-10°C'
    if (no==201.2): # 10°C
        fn = 'Scienlab__LTO20130205__LTO20130205-05__PowerTemperatureTest_dis-Umax-manual#3__Tamb=10°C'
    if (no==201.3): # 40°C
        fn = 'Scienlab__LTO20130205__LTO20130205-05__PowerTemperatureTest_dis-Umax-manual#3__Tamb=40°C'
    
    if (no==204.0):
        fn = 'Scienlab__SCiB 23Ah__SCiB-23Ah-01__PowerTemperatureTest_dis-Umax#1'
    if (no==205.0):
        fn = 'Scienlab__SCiB 23Ah__SCiB-23Ah-02__PowerTemperatureTest_dis-Umax#1'
    if (no==206.0):
        fn = 'Scienlab__ICR18650 HE2__HE2-1123-01__PowerTemperatureTest_dis-Umax#1'
    if (no==208.0):
        fn = 'Scienlab__ANR26650m1B__M1B-1223-01__PowerTemperatureTest_dis-Umax-inclOCV#1'
    if (no==210.0):
        fn = 'Scienlab__SLPB8644143__SLPB8644143-45__PowerTemperatureTest_dis-Umax-manual#1'
    if (no==211.0):
        fn = 'Scienlab__SLPB8644143__SLPB8644143-353__PowerTemperatureTest_dis-Umax-manual#1'
    if (no==212.0):
        fn = 'Scienlab__ICR18650 HE2__HE2-1123-03__PowerTemperatureTest_dis-Umax-manual#1'
    if (no==213.0):
        fn = 'Scienlab__ICR18650 HE2__HE2-1123-04__PowerTemperatureTest_dis-Umax-manual#1'
        
    res = pickle.load(open('{}/data/ess/{}.pickle'.format('.',fn),'rb'))
    
    
    #%% clean raw data
    ## dis
    if (fn=='Scienlab__LTO20130205__LTO20130205-03__PowerTemperatureTest_dis#3'):
        del_rn = {11: -3}
        
        for rn in del_rn:
            for key in res[rn]:
                if isinstance(res[rn][key],np.ndarray):
                    res[rn][key] = res[rn][key][:del_rn[rn]]
    elif (fn=='Scienlab__SLPB8644143__SLPB8644143-45__PowerTemperatureTest_dis#1'):
        del_rn = {1: -1,
                  4: -1,
                  12: -2}
        
        for rn in del_rn:
            for key in res[rn]:
                if isinstance(res[rn][key],np.ndarray):
                    res[rn][key] = res[rn][key][:del_rn[rn]]
    elif (fn=='Scienlab__SLPB8644143__SLPB8644143-45__PowerTemperatureTest_dis#2'):
        del_rn = {16: -2}
        
        for rn in del_rn:
            for key in res[rn]:
                if isinstance(res[rn][key],np.ndarray):
                    res[rn][key] = res[rn][key][:del_rn[rn]]
    elif (fn=='Scienlab__SLPB8644143__SLPB8644143-353__PowerTemperatureTest_dis#1'):
        del_rn = {8: -3,
                  20: -3,
                  24: -2}
        
        for rn in del_rn:
            for key in res[rn]:
                if isinstance(res[rn][key],np.ndarray):
                    res[rn][key] = res[rn][key][:del_rn[rn]]
    
    ## Variation of Umax
    if (fn=='Scienlab__SCiB 23Ah__SCiB-23Ah-01__PowerTemperatureTest_dis-Umax#1'):
        del_rn = {}
        del_rn['indU=0'] = {2: -1}
        
        for ref_key in del_rn:
            for rn in del_rn[ref_key]:
                for key in res[ref_key][rn]:
                    if isinstance(res[ref_key][rn][key],np.ndarray):
                        res[ref_key][rn][key] = res[ref_key][rn][key][:del_rn[ref_key][rn]]            
    if (fn=='Scienlab__ICR18650 HE2__HE2-1123-01__PowerTemperatureTest_dis-Umax#1'):
        del_rn = {}       
        # remove indices from power regimes
        for key in res:
            del_rn[key] = {0: 0}
        del_rn['indU=0'].update({12: -3}) 
        
        for ref_key in del_rn:
            for rn in del_rn[ref_key]:
                for key in res[ref_key][rn]:
                    if isinstance(res[ref_key][rn][key],np.ndarray):
                        res[ref_key][rn][key] = res[ref_key][rn][key][:del_rn[ref_key][rn]]
    
    ## Variation of Umax +OCV
    if (fn=='Scienlab__ANR26650m1B__M1B-1223-01__PowerTemperatureTest_dis-Umax-inclOCV#1'):
        del_rn = {}
        pop_rn = {}
        
        # remove indices from power regimes
        for key in res:
            if key != 'ocv':
                del_rn[key] = {0: 0} # each smallest

        del_rn['indU=2'].update({15: -4})
        
        # remove full power regimes
        for key in ['indU=2', 'indU=3']:
            pop_rn[key] = [k for k in range(len(res[key])) if k%2]

        pop_rn['indU=2'] += [43,44]
        pop_rn['indU=2'] = list(set(pop_rn['indU=2']))  # remove duplicates
                            
        for ref_key in del_rn:
            for rn in del_rn[ref_key]:
                for key in res[ref_key][rn]:
                    if isinstance(res[ref_key][rn][key],np.ndarray):
                        res[ref_key][rn][key] = res[ref_key][rn][key][:del_rn[ref_key][rn]]
        
        for ref_key in pop_rn:
            for rn in sorted(pop_rn[ref_key],reverse=True):
                res[ref_key].pop(rn)
            
            try: # adjust numbers after removing odd-numbered power regimes 
                for rn in range(len(res[ref_key])):
                    if ~(res[ref_key][rn]['ind_C_']%2):
                        res[ref_key][rn]['ind_C_'] = int(res[ref_key][rn]['ind_C_']/2)
            except:
                 0
    
    if (fn=='Scienlab__SLPB8644143__SLPB8644143-45__PowerTemperatureTest_dis-Umax-manual#1' or
        fn=='Scienlab__SLPB8644143__SLPB8644143-353__PowerTemperatureTest_dis-Umax-manual#1'):
        del_rn = {}
        pop_rn = {}
        
        # remove indices from power regimes
        if (fn=='Scienlab__SLPB8644143__SLPB8644143-45__PowerTemperatureTest_dis-Umax-manual#1'):
            del_rn['indU=0'] = {16: -2} 
        elif (fn=='Scienlab__SLPB8644143__SLPB8644143-353__PowerTemperatureTest_dis-Umax-manual#1'):
            for key in res:
                if key != 'ocv':
                    del_rn[key] = {0: 0} # each smallest
            del_rn['indU=0'].update({9: -3}) 
            
        # remove full power regimes
        for key in ['indU=0']:
            pop_rn[key] = [k for k in range(len(res[key])) if k%2]
        
        for ref_key in del_rn:
            for rn in del_rn[ref_key]:
                for key in res[ref_key][rn]:
                    if isinstance(res[ref_key][rn][key],np.ndarray):
                        res[ref_key][rn][key] = res[ref_key][rn][key][:del_rn[ref_key][rn]]
        
        for ref_key in pop_rn:
            for rn in sorted(pop_rn[ref_key],reverse=True):
                res[ref_key].pop(rn)
            
            try: # adjust numbers after removing odd-numbered power regimes 
                for rn in range(len(res[ref_key])):
                    if ~(res[ref_key][rn]['ind_C_']%2):
                        res[ref_key][rn]['ind_C_'] = int(res[ref_key][rn]['ind_C_']/2)
            except:
                 0
        
    return res


#%%
def adjust_res_keys(res,prefix_str='cell_model[1,1].',mode='dis'):
    """
    Adjusts keys and processes data within the results list from energy storage 
    system experiments. This includes renaming keys, calculating additional 
    parameters, and setting plot limits based on experimental settings.

    Parameters:
    - res (list): A list of dictionaries, each containing data from a single experiment.
    - prefix_str (str, optional): Prefix for keys in the result dictionaries to 
      indicate specific sub-model data. Defaults to 'cell_model[1,1].'.
    - mode (str, optional): Specifies the processing mode. Either 'dis' 
      (discharge) or 'chg' (charge). Defaults to 'dis'.

    Returns:
    - list: The same list of dictionaries with adjusted keys and additional 
      parameters for further analysis.

    Note:
    - Assumes a specific structure for the input data. Adaptations may be necessary 
      for different data formats or experimental setups.
    - Intended for use within a broader toolkit for energy storage system analysis.
    - Mode 'chg' or 'dis' determines the sign convention for currents, 
      power, and Crate. Raises a `ValueError` if `mode` is not 'dis' or 'chg'.
    """
    
    if mode not in ['dis','chg']:
        raise ValueError('[Error] Argument of mode not in ["dis","chg"].')
    
    for n,_ in enumerate(res):
        res[n]['_nP'] = res[n].pop('ind_C_')
        res[n]['_nSoC'] = 0
        res[n]['time'] = res[n]['Time_s']
        
        # Data Arrays
        res[n][prefix_str+'E_cell'] = res[n]['E_J']/3600
        # res[n][prefix_str+'P_cell'] = res[n].pop('P_W')
        # res[n][prefix_str+'P_cell'] = [np.average(res[n]['P_W'])]*len(res[n]['P_W'])

        ### CAUTION: Sign Error in Publication P1 due to I_cut > 0. Cell Limits
        ###          specified by Plot Annotations differ slightly. Graphical
        ###          Synthesis does not reveal a backward Intersection between
        ###          upper and lower Limits when decreasing Pcell. However,
        ###          apart from this, all Key Statements remain the same.
        ###          Plots can be reproduced using: 'tmp_Pmin = -tmp_Pmin'
        tmp_Pmin = abs(res[n]['_settings']['ParameterValues']['U_max'][0]*res[n]['_settings']['ParameterValues']['I_cut'][0])
        tmp_Pmax = abs(res[n]['_settings']['ParameterValues']['U_max'][0]*res[n]['_settings']['ParameterValues']['I_dis_max'][0])
        # tmp_Pmin = -tmp_Pmin
        
        ref_n_max = 10 if len(res) <=10 else 25 if len(res)-1 <=25 else 50 # CAUTION HERE!
        res[n][prefix_str+'P_cell'] = [tmp_Pmin+n/ref_n_max*(tmp_Pmax-tmp_Pmin)]*len(res[n]['P_W'])

        res[n][prefix_str+'U_cell'] = res[n].pop('U_V')
        res[n][prefix_str+'I_cell'] = res[n].pop('I_A')
        res[n][prefix_str+'Crate_cell'] = res[n][prefix_str+'I_cell']/res[n]['_settings']['ParameterValues']['Q_n'][0]
        res[n][prefix_str+'Crate'] = res[n][prefix_str+'I_cell']/res[n]['_settings']['ParameterValues']['Q_n'][0]


        res[n][prefix_str+'T_cell'] = res[n].pop('T_2_C')

        # key_Qref = [key for key in res[n]['_settings']['ParameterValues'] if 'Q_SCT' in key][-1]
        key_Qref = 'Q_n'
        if mode == 'chg':
            res[n][prefix_str+'SoC_cell'] = res[n]['Qpos_As']/3600/res[n]['_settings']['ParameterValues'][key_Qref][0]
        else:
            res[n][prefix_str+'SoC_cell'] = 1 + res[n]['Qneg_As']/3600/res[n]['_settings']['ParameterValues'][key_Qref][0]

        # Limits   
        res[n]['U_max'] = res[n]['_settings']['ParameterValues']['U_max']
        res[n]['U_min'] = res[n]['_settings']['ParameterValues']['U_min']
        # res[n]['_settings']['ParameterValues']['I_dis_max'] = res[n]['_settings']['ParameterValues']['I_dis_max']
        res[n]['I_dis_max'] = res[n]['_settings']['ParameterValues']['I_dis_max']
        res[n]['I_chg_max'] = res[n]['_settings']['ParameterValues']['I_chg_max']
        res[n]['I_cut'] = res[n]['_settings']['ParameterValues']['I_cut']        
        
        res[n]['_settings']['ParameterValues']['Crate_max'] = res[n]['_settings']['ParameterValues']['I_dis_max']/res[n]['_settings']['ParameterValues']['Q_n'][0]
        res[n]['Crate_max'] = res[n]['_settings']['ParameterValues']['Crate_max']
        
        res[n]['_settings']['ParameterValues']['T_max'] = res[n]['_settings']['ParameterValues']['T_dis_max']
        res[n]['T_max'] = res[n]['_settings']['ParameterValues']['T_max']
        
        res[n]['_settings']['ParameterValues']['T_min'] = res[n]['_settings']['ParameterValues']['T_dis_min']
        res[n]['T_min'] = res[n]['_settings']['ParameterValues']['T_min']
                
        # change sign
        for key in [prefix_str+'E_cell',
                    prefix_str+'P_cell',
                    prefix_str+'I_cell',
                    'I_dis_max',
                    prefix_str+'Crate_cell',
                    'Crate_max',
                    ]:
            if mode == 'chg':
                res[n][key] = -np.abs(res[n][key])
            else:
                res[n][key] = np.abs(res[n][key])


        # plot_limits
        lim_scale = 1.1
        res[n]['lim_E_dis_max'] = np.ceil(res[n]['_settings']['ParameterValues']['Q_n'][0]*res[n]['_settings']['ParameterValues']['U_max'][0]*lim_scale)
        res[n]['lim_P_dis_max'] = np.ceil(np.abs(res[n]['_settings']['ParameterValues']['I_dis_max'][0]*res[n]['_settings']['ParameterValues']['U_max'][0])*lim_scale)
            
    return res


#%%
def set_app(no):
    """
    Selects and returns predefined application settings for power grid supply
    scenarios and energy storage systems. Each setting includes details such
    as name, source, required active power (Pac_req), required energy capacity
    (Eac_req), and the calculated energy-to-power ratio (E/P_req).

    Parameters:
    - no (int): Identifier number for selecting a specific set of application
      settings.

    Returns:
    - dict: A dictionary containing the selected application settings. Keys
      include 'name', 'source', 'Pac_req', 'Eac_req', and 'E/P_req'.

    Note:
    - Predefined settings for various scenarios are identified by unique
      integers 'no'.
    - Includes settings for different power grid supply scenarios and energy
      storage systems, with references to external sources.
    - Calculates the 'E/P_req' ratio for the selected setting.
    - If an unrecognized 'no' is provided, no settings are returned.
    """

    ## ESS Examples used in Publication (Case Studies)
    # DSS Explanation
    if (no==0):
        d = {'name':             'Power Grid Supply #1',
             'source':           '',
             'Pac_req':          820000,
             'Eac_req':          330000}
    # Case Study: DSS available
    if (no==1):
        d = {'name':             'Power Grid Supply #2',
             'source':           '',
             'Pac_req':          820000,
             'Eac_req':          75000}
    # Case Study: no DSS available --> Limit Variation
    if (no==2):
        d = {'name':             'Power Grid Supply #3',
             'source':           '',
             'Pac_req':          820000,
             'Eac_req':          400000}
    # Case Study: no DSS ever available
    if (no==3):
        d = {'name':             'Power Grid Supply #4',
             'source':           '',
             'Pac_req':          820000,
             'Eac_req':          900000}
    
    ## ESS Examples extracted from Literature
    if (no==100):
        d = {'name':             '250 kW/500 kWh Li-ion battery integrated with the grid and solar farm - Qatar',
              'source':          'https://doi.org/10.1016/j.jpowsour.2017.10.063',
              'Pac_req':         250000,
              'Eac_req':         500000}
    if (no==101):
        d = {'name':             '50 kW / 60 kWh Energy Storage System - BYD',
              'source':          'https://www.energystorageexchange.org',
              'Pac_req':         50000,
              'Eac_req':         60000}
    if (no==102):
        d = {'name':             'Tehachapi Wind Energy Storage Project - Southern California Edison',
              'source':          'https://www.energystorageexchange.org',
              'Pac_req':         8000000,
              'Eac_req':         32000000}
    if (no==103):
        d = {'name':             'Falköping Substation Smart Grid - ABB',
              'source':          'https://www.energystorageexchange.org',
              'Pac_req':         75000,
              'Eac_req':         75000}
    if (no==104):
        d = {'name':             'SDG&E Century Park DESS - Greensmith Energy',
             'source':           'https://www.energystorageexchange.org',
             'Pac_req':          50000,
             'Eac_req':          82000}

    # calculate E/P ratio
    d['E/P_req'] = d['Eac_req']/d['Pac_req']
    return d


#%%
def set_pec(no):
    """
    Selects and returns predefined settings for various models of SMA Solar
    Power Electronics Converters (PECs). The function provides detailed
    specifications including nominal and maximum power capacities, efficiency,
    current and voltage limits for different PEC models.
    
    Parameters:
    - no (int): An identifier number used to select a specific PEC model.
    
    Returns:
    - dict: A dictionary containing the selected PEC's settings. Includes
      'name', 'Pac_nom', 'Pac_max', 'Pdc_max', 'eta_max', 'Idc_max', 'Uac_nom',
      'Uac_max', 'Uac_min', 'Udc_nom', 'Udc_max', 'Udc_min', and other 
      operational parameters.
    
    Note:
    - Predefined settings for various SMA Solar PEC models, each identified by
      a unique integer 'no'.
    - Detailed specifications for each model are provided, suitable for 
      integration into energy system designs.
    - Additional operational parameters like 'n' and 'eta_op' are also initialized.
    - If an unrecognized 'no' is provided, the function will not return any settings.
    """

    # Source: SMA Solar Technology AG. SUNNY CENTRAL STORAGE 500 / 630 / 720 / 760 /800 / 850 / 900 / 1000 - Inverter for Large-Scale Battery Storage Systems. Niestetal, Germany: 2020.
    if (no==0):
        d = {'name':            'SMA Solar SCS500',
             'Pac_nom':         500000,
             'Pac_max':         550000,
             'Pdc_max':         560000,
             'eta_max':         0.986,
             'Idc_max':         1400,
             'Uac_nom':         270,
             'Uac_max':         310,
             'Uac_min':         243,
             'Udc_nom':         449,
             'Udc_max':         850,
             'Udc_min':         430}
    if (no==1):
        d = {'name':            'SMA Solar SCS630',
             'Pac_nom':         630000,
             'Pac_max':         700000,
             'Pdc_max':         713000,
             'eta_max':         0.987,
             'Idc_max':         1400,
             'Uac_nom':         315,
             'Uac_max':         362,
             'Uac_min':         284,
             'Udc_nom':         529,
             'Udc_max':         850,
             'Udc_min':         500}
    if (no==2):
        d = {'name':            'SMA Solar SCS720',
             'Pac_nom':         720000,
             'Pac_max':         792000,
             'Pdc_max':         808000,
             'eta_max':         0.986,
             'Idc_max':         1400,
             'Uac_nom':         324,
             'Uac_max':         372,
             'Uac_min':         292,
             'Udc_nom':         577,
             'Udc_max':         850,
             'Udc_min':         480}
    if (no==3):
        d = {'name':            'SMA Solar SCS760',
             'Pac_nom':         760000,
             'Pac_max':         836000,
             'Pdc_max':         853000,
             'eta_max':         0.986,
             'Idc_max':         1400,
             'Uac_nom':         342,
             'Uac_max':         393,
             'Uac_min':         308,
             'Udc_nom':         609,
             'Udc_max':         850,
             'Udc_min':         505}
    if (no==4):
        d = {'name':            'SMA Solar SCS800',
             'Pac_nom':         800000,
             'Pac_max':         880000,
             'Pdc_max':         898000,
             'eta_max':         0.986,
             'Idc_max':         1400,
             'Uac_nom':         360,
             'Uac_max':         414,
             'Uac_min':         324,
             'Udc_nom':         641,
             'Udc_max':         950,
             'Udc_min':         530}
    if (no==5):
        d = {'name':            'SMA Solar SCS850',
             'Pac_nom':         850000,
             'Pac_max':         935000,
             'Pdc_max':         954000,
             'eta_max':         0.986,
             'Idc_max':         1400,
             'Uac_nom':         386,
             'Uac_max':         443,
             'Uac_min':         348,
             'Udc_nom':         681,
             'Udc_max':         950,
             'Udc_min':         568}
    if (no==6):
        d = {'name':            'SMA Solar SCS900',
             'Pac_nom':         900000,
             'Pac_max':         990000,
             'Pdc_max':         1010000,
             'eta_max':         0.986,
             'Idc_max':         1400,
             'Uac_nom':         405,
             'Uac_max':         465,
             'Uac_min':         365,
             'Udc_nom':         722,
             'Udc_max':         950,
             'Udc_min':         596}
    if (no==7):
        d = {'name':            'SMA Solar SCS1000',
             'Pac_nom':         1000000,
             'Pac_max':         1100000,
             'Pdc_max':         1122000,
             'eta_max':         0.987,
             'Idc_max':         1400,
             'Uac_nom':         405,
             'Uac_max':         465,
             'Uac_min':         365,
             'Udc_nom':         688,
             'Udc_max':         900,
             'Udc_min':         596}

    
    # Parameter Initialization (n > 1 recommended to avoid errors when creating plot objects for the first time)
    d['n'] = 2
    
    d['eta_op'] = d['eta_max']
    d.update(_fcn_eta(d['name']))
    
    return d


#%%
def _fcn_eta(name_pec):
    """
    Retrieves and processes the efficiency curve for a specified Power Converter 
    model, particularly for SMA Solar models. The function reads efficiency data 
    from an external file and creates interpolation functions based on the power 
    ratio and efficiency values. It supports different voltage conditions (Umin 
    and Umax).
    
    Parameters:
    - name_pec (str): The name of the Power Converter model.
    
    Returns:
    - dict: A dictionary containing interpolation functions for the efficiency 
      curve under different voltage conditions. Keys include 'fcn_eta(Umin)' and 
      'fcn_eta(Umax)'.
    
    Raises:
    - ValueError: If no efficiency curve is available for the selected Power 
      Converter model.
    
    Note:
    - Currently, the function is tailored for SMA Solar inverter models.
    - It relies on external data files containing the efficiency curves, expected
      to be in a specific format.
    - Interpolation functions are created to estimate efficiency for given power
      ratios, aiding in the accurate modeling of converter behavior.
    """
    
    d = {}

    if 'SMA Solar' in name_pec:
        # Source: SMA Solar Technology AG. SUNNY TRIPOWER STORAGE 60 - Highest power density for flexible applications. Niestetal, Germany: 2022.
        fn = 'SMA_STPS60-DS.xls'
    else:
        raise ValueError('No Efficiency Curve availble for selected Power Converter: "{}".'.format(name_pec))
    
    tmp_xls = xls2dict_digidata.read_xls(r'{}/data/pec/{}'.format('.',fn))

    tmp_d = tmp_xls['SMA_STPS60-DS']['eta(ratioP)']
    tmp_x = 'ratioP_-'
    tmp_y = 'eta_-'

    for n,key in enumerate(tmp_d):
        if n==1:
            continue
        
        notnan = ~np.isnan(tmp_d[key][tmp_x])
        if n==0:
            d['fcn_eta(Umin)'] = interpolate.interp1d(np.array(tmp_d[key][tmp_x])[notnan],np.array(tmp_d[key][tmp_y])[notnan])
        else:
            d['fcn_eta(Umax)'] = interpolate.interp1d(np.array(tmp_d[key][tmp_x])[notnan],np.array(tmp_d[key][tmp_y])[notnan])
    
    return d


#%%
def _set_by_arg(fcn,key,arg):
    """
    Iterates through settings generated by a function (`fcn`), searching for a 
    key-value pair match. When it finds a dictionary where the specified key matches 
    the given argument, it returns that dictionary of settings.

    Parameters:
    - fcn (function): A function returning a dictionary of settings when called 
      with an integer.
    - key (str): The key to search within the dictionaries returned by `fcn`.
    - arg: The value to match with the specified key.

    Returns:
    - dict: Dictionary with settings where the specified key matches `arg`. If no 
      match is found, returns an empty dictionary.

    Note:
    - Calls `fcn` iteratively with integer arguments from 0 to 99.
    - Checks each dictionary returned by `fcn` for the presence of `key` and 
      compares its value to `arg`.
    - Returns the first matching dictionary or continues searching until 100 
      iterations or an exception occurs.
    - Useful for searching through predefined lists of settings or configurations.
    """
    tmp_d = {}
    for i in range(100):
        try:
            if key in fcn(i).keys():
                if fcn(i)[key] == arg:
                    tmp_d.update(fcn(i))
            else:
                continue
        except:
            break
    return tmp_d


# #%%
# def load_DriveCycle(fname):
#     """
#     Loads and parses drive cycle data from a specified file into a structured 
#     dictionary. The file is expected to follow a specific format, optionally 
#     including headers to denote data columns.

#     Parameters:
#     - fname (str): File name or path containing drive cycle data.

#     Returns:
#     - dict: Dictionary containing parsed drive cycle data. Keys correspond to 
#       different sections or categories in the file, each with a sub-dictionary 
#       of data series or arrays based on the file columns.

#     Note:
#     - Expected file format: data separated by spaces, with optional headers 
#       prefixed by '#'.
#     - If headers are present, indicated by '#', they are used as keys in the 
#       sub-dictionaries; otherwise, numerical indices are used.
#     - Utilizes Python's `eval` for numerical data parsing, caution advised 
#       with untrusted files.
#     """

#     tmp_d = {}
#     with open(fname, mode='r', closefd=True) as f:
#         content = f.readlines()
#         for x in content:
#             row = x.split()
#             if len(row) == 0 or row[0] == '#1':
#                 if 'tmp_header' in locals():
#                     del tmp_header
#                 continue
#             elif row[0] == 'double':
#                 tmp_key = '{}'.format(row[1])
#                 tmp_d[tmp_key] = {}
                
#                 # if tableName given
#                 if '#' in x:
#                     tmp_header = x.split('#')[-1]
#                     tmp_header = tmp_header.split('\n')[0]
#                     tmp_header = tmp_header.split(',')                    
#                     for n,ele in enumerate(tmp_header):
#                         tmp_header[n] = ele.strip()
#                         tmp_d[tmp_key][ele.strip()] = []
#             else:
#                 tmp_list = list(eval(row[0]))
#                 for app_load in range(len(tmp_list)):
#                     if 'tmp_header' in locals():
#                         tmp_d[tmp_key][tmp_header[app_load]].append(tmp_list[app_load])
#                     else:
#                         # if tableName not given
#                         if app_load not in tmp_d[tmp_key].keys():
#                             tmp_d[tmp_key][app_load] = []
#                         tmp_d[tmp_key][app_load].append(tmp_list[app_load])
#     return tmp_d