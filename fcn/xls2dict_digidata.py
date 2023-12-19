"""
xls2dict_digidata Module

Author:     Sven Wiegelmann,
            Leibniz University Hannover,
            Institute of Electric Power Systems,
            Electric Energy Storage Systems Section,
            Appelstraße 9A,
            30167 Hannover,
            Germany
Version:    12.12.2023

Overview:
This module, xls2dict_digidata, specializes in converting Excel (.xls) files into 
structured Python dictionaries. It is adept at handling diverse data layouts within 
Excel sheets, providing a robust tool for Python-based data processing and analysis.

Features:
- Efficient transformation of Excel data into accessible dictionary format.
- Capability to handle complex data structures within Excel sheets.
- Selective sheet processing, enabling targeted data extraction.
- Support for nested data arrangements and multiple data types.

Usage:
This module is ideal for users requiring Excel data manipulation in Python, especially 
in fields like data analysis, scientific research, and engineering. It simplifies the 
integration of Excel-based data into Python workflows.

Execution:
Execute in a Python environment with 'xlrd' package installed. The module is designed 
to be imported and used in conjunction with other data processing or analysis scripts.

Notes:
- Specifically developed for .xls format; adaptations may be required for other Excel formats.
- Structured to be versatile and adaptable across various data extraction needs.
- Ensure understanding of the Excel data structure for optimal use of the module's capabilities.
"""

import xlrd
import os

#%%
def read_xls(workbook_url, NN_Sheet=[]):
    """
    Converts data from an Excel (.xls) file into a structured dictionary format. 
    This function reads each sheet within the workbook, transforming the data 
    into nested dictionaries for easy accessibility and manipulation in Python.

    Parameters:
    - workbook_url (str): Path to the .xls file to be read.
    - NN_Sheet (list, optional): List of sheet names to be excluded from processing.
      Defaults to an empty list, meaning all sheets are processed.

    Returns:
    - dict: A dictionary where each key represents a sheet name (with spaces removed),
      and its value is another dictionary representing the sheet's content. The nested
      dictionaries are structured based on the columns (axes, lines, headers) of the
      Excel sheet, making the data easy to navigate and use in Python.

    Note:
    - The function assumes a specific layout of the Excel sheets, with the first
      three rows containing the axes, lines, and headers, respectively.
    - Empty cells are converted to NaN (Not a Number) to maintain consistency in data type.
    - If a sheet does not adhere to the expected layout (i.e., missing single column headers),
      the function will print an error and stop processing that sheet.
    - Applicable in data analysis, research, or any scenario requiring Excel data
      integration into Python environments.
    """
    workbook_dict = {}
    book = xlrd.open_workbook(workbook_url, logfile=open(os.devnull, 'w'))
    sheets = book.sheets()
    for sheet in sheets:
        # continue if sheet-data is not needed
        if sheet.name in NN_Sheet:
            continue
        else:
            num_axis = 0
            num_line = 1
            num_head = 2
            skip = 0
            
        # read axis column
        axes = sheet.row_values(num_axis)
        if all(ele == '' for ele in axes):
            axes[0] = 'NaN'
        # read line column
        lines = sheet.row_values(num_line)
        if all(ele == '' for ele in lines):
            lines[0] = 'NaN'
        # read header column
        headers = sheet.row_values(num_head)
        if '' in headers:
            print('E: Sheet "{}" not readable. Single Column Headers needed!'.format(sheet.name))
            break
        
        # create array
        rows = []
        # read in every row and save into rows
        for row_index in range(num_head+1+skip, sheet.nrows):
            row = [float("nan") if i == '' else i for i in sheet.row_values(row_index)]
            rows.append(row)
    
        # create dictionary for sheet
        sheet_dict = {}
        tmp_ax = ''
        tmp_line = ''

        i = 0         
        # iterate through headers of sheet
        for head in headers:
            if axes[i] != '':
                tmp_ax = axes[i]
            if lines[i] != '':
                tmp_line = lines[i]
            
            if axes[i] == tmp_ax and not check_dict(tmp_ax.replace(' ',''),sheet_dict,[sheet.name]):
                sheet_dict.update({tmp_ax.replace(' ',''): {}})
            if lines[i] == tmp_line and not check_dict(tmp_line.replace(' ',''),sheet_dict[tmp_ax.replace(' ','')],[sheet.name,tmp_ax]):
                sheet_dict[tmp_ax.replace(' ','')].update({tmp_line.replace(' ','') : {}})
            if not check_dict(headers[i].replace(' ',''),sheet_dict[tmp_ax.replace(' ','')][tmp_line.replace(' ','')],[sheet.name,tmp_ax,tmp_line]):
                sheet_dict[tmp_ax.replace(' ','')][tmp_line.replace(' ','')].update({headers[i].replace(' ',''): [rows[k][i] for k in range(len(rows))]})
            i+=1
    
        # save dictionary of sheet into workbook and delete spaces in key
        workbook_dict[sheet.name.replace(' ','')] = sheet_dict
        
    return workbook_dict


#%%
def check_dict(tmp_key, tmp_dict, tmp_sheet):
    """
    Checks if a given key exists in a specified dictionary. This function is 
    utilized during the processing of Excel data to ensure uniqueness of keys 
    when constructing nested dictionaries.
    
    Parameters:
    - tmp_key (str): The key to be checked in the dictionary.
    - tmp_dict (dict): The dictionary in which the key's existence is checked.
    - tmp_sheet (list): List containing the context information (like sheet name) 
      for the key, used for logging purposes.
    
    Returns:
    - int: Returns 1 if the key exists in the dictionary, and 0 otherwise.
    
    Note:
    - This function aids in maintaining data integrity and avoiding key collisions 
      in dictionaries.
    - It's primarily used in conjunction with the read_xls function for handling 
      Excel data.
    - When a key collision occurs, the function prints a warning with the problematic 
      key and its corresponding context information for debugging.
    """

    if tmp_key in tmp_dict:
        i = 1
        print('W: Key "{}" already exists in "{}".'.format(tmp_key,'/'.join(tmp_sheet)))
    else:
        i = 0
    return i