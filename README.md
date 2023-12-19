# System Design Method (SDM) Tool

[Graphical Abstract: XXX]

## Description
The System Design Method (SDM) tool is a Python-based interactive suite designed
for the analysis and design of energy storage systems. It integrates various
functionalities to facilitate dynamic visualization, manipulation, and
understanding of energy system parameters, offering an immersive user experience
for system designers, researchers, and educators in the field of energy
engineering.

## Associated Work

The associated research paper is available at
<still in peer review>

## Features
- Interactive Extended Ragone Plot (ERP) Generation: Visualize and analyze the
  interplay between energy and power parameters in storage systems.
- Constraint and Performance Visualization: Dynamically display system
  constraints and performance metrics.
- Parameter Adjustment via Interactive Sliders: Real-time adjustment of system
  parameters like power, energy, and efficiency for on-the-fly analysis.
- Integrated Case Studies: Practical examples and scenarios to illustrate the
  application of the tool in real-world system designs.
- Graphical Synthesis: Aids in decision-making by visually synthesizing system
  design parameters.

## Installation
1. Clone the Repository: 
   Use the command 'git clone' followed by your repository URL.
2. Install Required Dependencies:
   Ensure you have Python installed, along with libraries like matplotlib, numpy,
   and xlrd. Use 'pip install' to install these libraries.

## Usage
This tool is intended for use by individuals with a basic understanding of Python
and familiarity with energy storage concepts. 
- Run the main script in a Python environment: Use 'python SDM_main.py' to
  execute the script.
- Interact with the tool through the provided graphical interface and sliders to
  explore different system designs.

## Modules
- SDM_main.py: The core interactive module integrating SDM functionalities.
- SDM_fcn.py: Contains functions and methods for system analysis and design.
- SDM_data.py: Manages data loading and preprocessing for the tool.
- xls2dict_digidata.py: Converts Excel data into Python dictionaries, aiding in
  data management.

## Documentation
For detailed information and usage instructions, refer to the docstrings within
each submodule and function.

## Contributing
Contributions to the SDM tool are welcome. If you have suggestions for
improvements or want to contribute code, please feel free to create an issue
or submit a pull request.

## License
This software is licensed under GPLv3, excluding later versions.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

For details see [XXX].

GPLv3 explicitely allows a commercial usage without any royalty or further
implications. However, any contact to discuss possible cooperations is
appreciated.

## Author
sdm - System Design Method\
Copyright (C) 2023\
Sven Wiegelmann\
wiegelmann@ifes.uni-hannover.de

Leibniz Universität Hannover\
Institut für Elektrische Energiesysteme\
Fachgebiet für Elektrische Energiespeichersysteme

Leibniz University Hannover\
Institute of Electric Power Systems\
Electric Energy Storage Systems Section

https://www.ifes.uni-hannover.de/ees.html
