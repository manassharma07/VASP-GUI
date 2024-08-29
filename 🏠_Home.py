import streamlit as st
import platform

# Set page config
st.set_page_config(page_title='GUI for VASP', layout='wide', page_icon="⚛️",
menu_items={
         'About': "A web app to help you create input files for running KS-DFT calculations with the [VASP](https://www.vasp.at/)"
     })

# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write(' The web app is made By [Manas Sharma](https://manas.bragitoff.com)')
st.sidebar.write(' in the groups of [Prof. Ananth Govind Rajan](https://www.agrgroup.org) and [Prof. Sudeep Punnathanam](https://chemeng.iisc.ac.in/sudeep/).')
st.sidebar.write('### *Powered by*')
st.sidebar.write('* [Py3Dmol](https://3dmol.csb.pitt.edu/) for Chemical System Visualizations')
st.sidebar.write('* [Streamlit](https://streamlit.io/) for making of the Web App')
st.sidebar.write('* [PyMatgen](https://pymatgen.org/) for Periodic Structure Representations')
st.sidebar.write('* [PubChempy](https://pypi.org/project/PubChemPy/1.0/) for Accessing the PubChem Database')
st.sidebar.write('* [MP-API](https://pypi.org/project/mp-api/) for Accessing the Materials Project Database')
st.sidebar.write('* [ASE](https://wiki.fysik.dtu.dk/ase/) for File Format Conversions')
st.sidebar.write('### *Useful links*')
st.sidebar.write('[Web App Source Code](https://github.com/manassharma07/VASP-GUI)')
st.sidebar.write('[VASP Wiki](https://www.vasp.at/wiki/index.php/The_VASP_Manual)')
st.sidebar.write('[VASP Official Website](https://www.vasp.at/)')
st.sidebar.write('[VASP Forum](https://www.vasp.at/forum/)')


# Main app
st.header('VASP GUI')
st.write('#### A web app to help you create input files for running KS-DFT calculations with the [VASP](https://www.vasp.at/)')

st.write('''The Vienna Ab initio Simulation Package (VASP) is a computer program for atomic scale materials modelling, e.g. electronic structure calculations and quantum-mechanical molecular dynamics, from first principles.''')
# "to create input (control) files from a given [CIF](https://en.wikipedia.org/wiki/Crystallographic_Information_File), [XYZ](https://en.wikipedia.org/wiki/XYZ_file_format), etc"
# "You can also parse the output of the RIPER program for further analysis."
text_intro = """
This app allows you to create `VASP` input files for a material in the [materialsproject.org](https://next-gen.materialsproject.org/) database or a molecule in the [PubChem database](https://pubchem.ncbi.nlm.nih.gov/). 
Additionally it allows to convert `CIF`, `XYZ`, `PWSCF` and other popular formats to `POSCAR` format.

You can select from a variety of options from the menu in the left sidebar."""
st.write(text_intro)

st.write('To run a DFT calculation, simply invoke `vasp` as')
st.code('mpirun -n N ./vasp_std', language='shell')
st.write('to run with `N` processors.')
