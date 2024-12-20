import streamlit as st
from mp_api.client import MPRester
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.cif import CifParser
from pymatgen.io.xyz import XYZ
from pymatgen.io.vasp.inputs import Poscar
from pymatgen.io.pwscf import PWInput
from pymatgen.core import Structure, Lattice
from ase import Atoms
import numpy as np

import py3Dmol
import pandas as pd
import streamlit.components.v1 as components
from io import StringIO
from ase.io.espresso import read_espresso_in
from ase.io.extxyz import read_extxyz
from ase.io.cif import read_cif
from ase.io import read
from ase.io.vasp import write_vasp
from pymatgen.io.ase import AseAtomsAdaptor
from ase.io.dmol import read_dmol_car

# Set page config
st.set_page_config(page_title='CIF/XYZ/CAR/PWSCF ➡️ VASP POSCAR', layout='wide', page_icon="⚛️",
                   menu_items={
                       'About': "A web app to help you convert various file formats to VASP POSCAR format."
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



# Function to convert a structure to CIF
def convert_to_cif(structure, filename):
    cif_writer = CifWriter(structure)
    cif_writer.write_file(filename)


# Function to convert a structure to POSCAR format and save to file
def convert_to_poscar_pymatgen(structure, filename):
    poscar = Poscar(structure)
    poscar.write_file(filename)
    with open(filename, 'r') as file:
        poscar_content = file.read()
    


# Function to convert a structure to POSCAR format and save to file
def convert_to_poscar_ase(structure, filename, direct=False):
    ase_atoms = structure.to_ase_atoms()
    with open(filename, 'w') as f:
        write_vasp(f, ase_atoms, direct=direct, sort=True)
    
    


# Function to visualize the structure using py3Dmol
def visualize_structure(structure, html_file_name='viz.html'):
    spin = st.checkbox('Spin', value=False, key='key' + html_file_name)
    view = py3Dmol.view(width=500, height=400)
    cif_for_visualization = structure.to(fmt="cif")
    view.addModel(cif_for_visualization, 'cif')
    view.setStyle({'sphere': {'colorscheme': 'Jmol', 'scale': 0.3},
                   'stick': {'colorscheme': 'Jmol', 'radius': 0.2}})
    view.addUnitCell()
    view.zoomTo()
    view.spin(spin)
    view.setClickable({'clickable': 'true'})
    view.enableContextMenu({'contextMenuEnabled': 'true'})
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
    HtmlFile.close()


# Function to visualize the structure using py3Dmol
def visualize_molecule(structure, html_file_name='viz.html'):
    spin = st.checkbox('Spin', value=False, key='key' + html_file_name)
    view = py3Dmol.view(width=500, height=400)
    xyz_for_visualization = structure.to(fmt="xyz")
    view.addModel(xyz_for_visualization, 'xyz')
    view.setStyle({'sphere': {'colorscheme': 'Jmol', 'scale': 0.3},
                   'stick': {'colorscheme': 'Jmol', 'radius': 0.2}})
    view.zoomTo()
    view.spin(spin)
    view.setClickable({'clickable': 'true'})
    view.enableContextMenu({'contextMenuEnabled': 'true'})
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
    HtmlFile.close()


# Function to display structure information
def display_structure_info(structure):
    st.subheader("Structure Information")
    st.write("Formula: ", structure.composition.reduced_formula)

    # Display lattice parameters
    a, b, c = structure.lattice.abc
    alpha, beta, gamma = structure.lattice.angles

    # Create a DataFrame for the lattice parameters and angles
    data = {
        "Lattice Parameters": [a, b, c, alpha, beta, gamma]
    }
    df_latt_params = pd.DataFrame(data, index=["a", "b", "c", "alpha", "beta", "gamma"])
    with st.expander("Lattice Parameters", expanded=False):
        st.table(df_latt_params)

    # Display lattice vectors
    lattice_vectors = structure.lattice.matrix
    df_vectors = pd.DataFrame(lattice_vectors, columns=["X", "Y", "Z"], index=["a", "b", "c"])
    with st.expander("Lattice Vectors", expanded=True):
        st.table(df_vectors)

    # Create a list of atomic coordinates
    with st.expander("Atomic Coordinates", expanded=False):
        coord_type = st.selectbox('Coordinate type', ['Cartesian', 'Fractional/Crystal'])
        if coord_type == 'Cartesian':
            atomic_coords = []
            for site in structure.sites:
                atomic_coords.append([site.species_string] + list(site.coords))
        else:
            atomic_coords = []
            for site in structure.sites:
                atomic_coords.append([site.species_string] + list(site.frac_coords))

        # Create a Pandas DataFrame from the atomic coordinates list
        df_coords = pd.DataFrame(atomic_coords, columns=["Element", "X", "Y", "Z"])

        # Display the atomic coordinates as a table
        st.table(df_coords)


def parse_cif_pymatgen(contents):
    # Parse the CIF file using pymatgen
    cif_parser = CifParser.from_str(contents)
    structure = cif_parser.get_structures(primitive=False)[0]  # Assuming there's only one structure in the CIF file
    return structure


# def parse_xyz(contents):
#     # Parse the XYZ file using pymatgen
#     xyz_parser = XYZ.from_str(contents)
#     structure = xyz_parser.molecule  # Assuming it's a molecule XYZ file
#     return structure

def parse_xyz(contents):
    # Step 1: Parse the XYZ file using pymatgen
    xyz_parser = XYZ.from_str(contents)
    molecule = xyz_parser.molecule  # Assuming it's a molecule XYZ file
    
    # Step 2: Convert pymatgen Molecule to ASE Atoms object
    ase_atoms = AseAtomsAdaptor.get_atoms(molecule)
    
    # Step 3: Calculate the new cell dimensions with a 15 Å buffer in each direction
    positions = ase_atoms.get_positions()
    min_pos = np.min(positions, axis=0)
    max_pos = np.max(positions, axis=0)
    buffer = 15.0
    cell_dimensions = (max_pos - min_pos) + 2 * buffer
    
    # Create a cubic cell with the calculated dimensions
    ase_atoms.set_cell(cell_dimensions)
    ase_atoms.center()  # Center the molecule in the cell
    ase_atoms.set_pbc(True)  # Set periodic boundary conditions
    
    # Step 4: Convert the ASE Atoms object back to pymatgen Structure
    structure = AseAtomsAdaptor.get_structure(ase_atoms)
    
    return structure


def parse_poscar(contents):
    # Parse the POSCAR file using pymatgen
    poscar_parser = Poscar.from_str(contents)
    structure = poscar_parser.structure
    return structure


def parse_quantum_espresso(contents):
    # Parse the POSCAR file using pymatgen
    qe_parser = PWInput.from_str(contents)
    structure = qe_parser.structure
    return structure


def parse_qe_ase(stringio):
    # Read QE input file
    atoms = read_espresso_in(stringio)

    # Convert ASE Atoms to pymatgen Structure
    structure = AseAtomsAdaptor().get_structure(atoms)
    return structure

def parse_extxyz_ase(stringio):
    # Read extended XYZ file
    atoms = read(stringio, format="extxyz")

    # Convert ASE Atoms to pymatgen Structure
    structure = AseAtomsAdaptor().get_structure(atoms)
    return structure


def parse_cif_ase(stringio):
    # Read CIF
    atoms = read(stringio, format="cif")

    # Convert ASE Atoms to pymatgen Structure
    structure = AseAtomsAdaptor().get_structure(atoms)
    return structure

def parse_car_ase(stringio):
    # Read CAR
    atoms = read_dmol_car(stringio)

    # Convert ASE Atoms to pymatgen Structure (determine if CAR file is 3D periodicity or not)
    if all(atoms.pbc):
        structure = AseAtomsAdaptor().get_structure(atoms)
    else:
        structure = AseAtomsAdaptor().get_molecule(atoms)
    return structure


# return filecontents
def read_file(filename):
    with open(filename, 'r') as file:
        return file.read()


# Streamlit app
st.write('# CIF/XYZ/CAR/PWSCF ➡️ VASP POSCAR')
st.write(
    "#### Convert various file formats to VASP POSCAR format.") 

st.write("Please select the file format")

# Select file format
file_format = st.selectbox("Select file format",
                           ("CIF", "XYZ", "CAR (Materials Studio)", "POSCAR", "Quantum ESPRESSO (PWSCF)", "Extended XYZ"))

if file_format == 'CIF':
    cif_parser_options = ['ASE', 'Pymatgen']
    cif_parser_selection = st.radio('Choose CIF parser', cif_parser_options)

uploaded_file = st.file_uploader("Choose a file", type=["cif", "xyz", "car", "POSCAR", "pwi", "extxyz"])
pasted_content = st.text_area('or paste the contents of the file here:', height=400)

contents = None

if uploaded_file is not None:
    contents = uploaded_file.read().decode('utf-8')
elif pasted_content:
    contents = pasted_content

if contents:

# if uploaded_file is not None:
#     # To read file as string:
#     contents = uploaded_file.read().decode('utf-8')

    if file_format == 'CIF' and cif_parser_selection == 'Pymatgen':
        structure = parse_cif_pymatgen(contents)
    elif file_format == 'CIF' and cif_parser_selection == 'ASE':
        stringio = StringIO(contents)
        structure = parse_cif_ase(stringio)
    elif file_format == 'XYZ':
        structure = parse_xyz(contents)
    elif file_format == 'CAR (Materials Studio)':
        stringio = StringIO(contents)
        structure = parse_car_ase(stringio)
    elif file_format == 'POSCAR':
        structure = parse_poscar(contents)
    elif file_format == 'Quantum ESPRESSO (PWSCF)':
        stringio = StringIO(contents)
        structure = parse_qe_ase(stringio)
    elif file_format == 'Extended XYZ':
        stringio = StringIO(contents)
        structure = parse_extxyz_ase(stringio)

    
    # Display structure information
    if isinstance(structure, Structure):  # type = Structure
        display_structure_info(structure)
    else:  # type = Molecule
        st.subheader("Atomic Coordinates")

    # Visualize the structure
    st.subheader("Structure Visualization")
    # if file_format == 'XYZ':
    #     visualize_molecule(structure, 'xyz_visualization.html')
    # else:
    #     visualize_structure(structure, 'structure_visualization.html')
    if isinstance(structure, Structure):  # type = Structure
        visualize_structure(structure, 'structure_visualization.html')
    else:
        visualize_molecule(structure, 'xyz_visualization.html')
    
    # Download CIF files
    if isinstance(structure, Structure):
        st.subheader("Download CIF Files")

        convert_to_cif(structure, "structure.cif")
        st.download_button('Download CIF', data=read_file("structure.cif"), file_name='structure.cif', key='cif_button')

    # Add a selection box for coordinate type
    coord_type = st.selectbox("Select coordinate type for POSCAR", ["Direct", "Cartesian"])
    if coord_type=='Direct':
        direct = True
    else:
        direct = False
    # Provide download link for POSCAR
    # convert_to_poscar_pymatgen(structure, 'POSCAR')
    convert_to_poscar_ase(structure, 'POSCAR', direct)
    poscar_content = read_file('POSCAR')

    # Generate a sample KPOINTS file
    kpoints_content = """Automatic generation
0
Monkhorst-pack
 4 4 4
 0 0 0
    """

    # Generate a sample INCAR file
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

    st.subheader("INPUT Files")

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
    st.warning("Please upload a file.")
