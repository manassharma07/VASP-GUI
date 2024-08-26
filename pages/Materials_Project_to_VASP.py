import streamlit as st
from mp_api.client import MPRester
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import py3Dmol
import pandas as pd
import streamlit.components.v1 as components
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from io import StringIO
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core import Structure

# Set page config
st.set_page_config(page_title='Materials Project ➡️ VASP', layout='wide', page_icon="⚛️",
menu_items={
         'About': "A web app to help you with DFT related calculations using [VASP](https://www.vasp.at/)"
     })

# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write('Originally Made By [Manas Sharma](https://manas.bragitoff.com)')
st.sidebar.write('In the group of [Prof. Dr. Marek Sierka](https://cmsg.uni-jena.de)')
st.sidebar.write('### *Powered by*')
st.sidebar.write('* [Py3Dmol](https://3dmol.csb.pitt.edu/) for Chemical System Visualizations')
st.sidebar.write('* [Streamlit](https://streamlit.io/) for making of the Web App')
st.sidebar.write('* [PyMatgen](https://pymatgen.org/) for Periodic Structure Representations')
st.sidebar.write('* [PubChempy](https://pypi.org/project/PubChemPy/1.0/) for Accessing the PubChem Database')
st.sidebar.write('* [MP-API](https://pypi.org/project/mp-api/) for Accessing the Materials Project Database')
st.sidebar.write('* [ASE](https://wiki.fysik.dtu.dk/ase/) for File Format Conversions')
st.sidebar.write('### *Contributors*')
st.sidebar.write('[Ya-Fan Chen ](https://github.com/Lexachoc)')
st.sidebar.write('### *Source Code*')
st.sidebar.write('[GitHub Repository](https://github.com/manassharma07/RIPER-Tools-for-TURBOMOLE)')

# Set your Materials Project API key
api_key = st.secrets["MP_API"]

# Function to convert a structure to POSCAR format and save to file
def convert_to_poscar_pymatgen(structure, filename):
    poscar = Poscar(structure)
    poscar.write_file(filename)
    with open(filename, 'r') as file:
        poscar_content = file.read()
    return poscar_content

# Function to generate VASP input files (KPOINTS, INCAR)
def generate_vasp_input_files(structure):
    # Generate POSCAR
    poscar_content = convert_to_poscar_pymatgen(structure, 'POSCAR')

    # Generate KPOINTS
    kpoints_content = """Automatic generation
0
Monkhorst-pack
 4 4 4
 0 0 0
    """

    # Generate INCAR
    incar_content = """
    # General settings
    ENCUT = 520       # Energy cutoff for plane waves
    PREC = Accurate   # Precision mode
    EDIFF = 1E-6      # Energy convergence criterion
    ISMEAR = 0        # Gaussian smearing
    SIGMA = 0.05      # Smearing width
    # Electronic relaxation
    ALGO = Fast       # Algorithm (Fast or Normal)
    NELM = 100        # Maximum number of electronic steps
    # Ionic relaxation
    IBRION = 2        # Algorithm for ionic relaxation (2 = conjugate gradient)
    NSW = 50          # Maximum number of ionic steps
    ISIF = 3          # Stress and relaxation (3 = relax cell shape and volume)
    # Output settings
    LWAVE = .FALSE.   # Write WAVECAR file
    LCHARG = .FALSE.  # Write CHGCAR file
    """

    return poscar_content, kpoints_content, incar_content

def parse_cif_ase(stringio):
    # Read CIF
    atoms = read(stringio, format="cif")

    # Convert ASE Atoms to pymatgen Structure
    structure = AseAtomsAdaptor().get_structure(atoms)
    return structure

# Function to visualize the structure using py3Dmol
def visualize_structure(structure, html_file_name='viz.html'):
    spin = st.checkbox('Spin', value=False, key='key'+html_file_name)
    view = py3Dmol.view(width=500, height=400)
    cif_for_visualization = structure.to(fmt="cif")
    view.addModel(cif_for_visualization, 'cif')
    view.setStyle({'sphere': {'colorscheme': 'Jmol', 'scale': 0.3},
                   'stick': {'colorscheme': 'Jmol', 'radius': 0.2}})
    view.addUnitCell()
    view.zoomTo()
    view.spin(spin)
    view.setClickable({'clickable':'true'})
    view.enableContextMenu({'contextMenuEnabled':'true'})
    view.show()
    view.render()
    t = view.js()
    f = open(html_file_name, 'w')
    f.write(t.startjs)
    f.write(t.endjs)
    f.close()

    HtmlFile = open(html_file_name, 'r', encoding='utf-8')
    source_code = HtmlFile.read() 
    components.html(source_code, height=300, width=500)

# Streamlit app
st.write('# Materials Project ➡️ VASP')
st.write("#### Generate VASP input files from structures obtained from the Materials Project.")

# Select Material ID
material_id = st.text_input("Enter the Material ID from Materials Project", value="mp-149")

# Fetch the structure from Materials Project
if material_id:
    with MPRester(api_key) as mpr:
        structure = mpr.get_structure_by_material_id(material_id)
    
    if structure:
        st.subheader("Structure Information")
        st.write("Formula: ", structure.composition.reduced_formula)
        visualize_structure(structure, 'structure_visualization.html')

        poscar_content, kpoints_content, incar_content = generate_vasp_input_files(structure)

        # Display POSCAR and KPOINTS in editable text boxes
        col1, col2 = st.columns(2)

        with col1:
            st.write("### POSCAR")
            poscar_editable = st.text_area("POSCAR Content", poscar_content, height=300)
            st.download_button(
                label="Download POSCAR",
                data=poscar_editable,
                file_name='POSCAR',
                mime='text/plain',
            )

        with col2:
            st.write("### KPOINTS")
            kpoints_editable = st.text_area("KPOINTS Content", kpoints_content, height=300)
            st.download_button(
                label="Download KPOINTS",
                data=kpoints_editable,
                file_name='KPOINTS',
                mime='text/plain',
            )

        # Display INCAR file
        st.subheader("Sample INCAR")
        incar_editable = st.text_area("INCAR Content", incar_content, height=300)
        st.download_button(
            label="Download INCAR",
            data=incar_editable,
            file_name='INCAR',
            mime='text/plain',
        )
    else:
        st.warning("Could not retrieve structure. Please check the Material ID.")
else:
    st.warning("Please enter a Material ID.")

