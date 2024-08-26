import streamlit as st
import pubchempy as pcp
from pymatgen.core.structure import Molecule
import py3Dmol
import pandas as pd
import streamlit.components.v1 as components
from pymatgen.io.ase import AseAtomsAdaptor
from ase.calculators.emt import EMT
from ase.calculators.singlepoint import SinglePointCalculator
from ase.optimize import BFGS
import requests
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.inputs import Poscar
from ase.io.vasp import write_vasp
import numpy as np
from pymatgen.core import Structure


# Set page config
st.set_page_config(page_title='PubChem ➡️ VASP', layout='wide', page_icon="⚛️",
menu_items={
         'About': "A web app to help you create VASP input files from PubChem entries."
     })

# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write(' The web app is made By [Manas Sharma](https://manas.bragitoff.com)')
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

# return filecontents
def read_file(filename):
    with open(filename, 'r') as file:
        return file.read()

# Function to convert a structure to POSCAR format and save to file
def convert_to_poscar_ase(structure, filename, direct=False):
    ase_atoms = structure.to_ase_atoms()
    with open(filename, 'w') as f:
        write_vasp(f, ase_atoms, direct=direct)

# Function to convert a structure to POSCAR format and save to file
def convert_to_poscar_pymatgen(structure, filename):
    poscar = Poscar(structure)
    poscar.write_file(filename)
    with open(filename, 'r') as file:
        poscar_content = file.read()
    return poscar_content

# Function to generate VASP input files (KPOINTS, INCAR)
def generate_vasp_input_files(structure, direct=False):
    # Generate POSCAR
    convert_to_poscar_ase(structure, 'POSCAR', direct=direct)
    poscar_content = read_file('POSCAR')
    # Generate KPOINTS
    kpoints_content = """Automatic generation
0
Gamma
 1 1 1
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
# CID to XYZ
def generate_xyz_coordinates(cid):
    compound = pcp.Compound.from_cid(cid, record_type='3d')
    atoms = compound.atoms
    coords = [(atom.x, atom.y, atom.z) for atom in atoms]

    num_atoms = len(atoms)
    xyz_text = f"{num_atoms}\n{compound.cid}\n"

    for atom, coord in zip(atoms, coords):
        atom_symbol = atom.element
        x, y, z = coord
        xyz_text += f"{atom_symbol} {x} {y} {z}\n"

    return xyz_text

def create_structure_with_vacuum(molecule, vacuum):
    # Get the max and min coordinates in each dimension
    coords = molecule.cart_coords
    min_coords = np.min(coords, axis=0)
    max_coords = np.max(coords, axis=0)
    
    # Calculate the molecule's dimensions
    dimensions = max_coords - min_coords
    
    # Calculate the new box dimensions with vacuum
    new_dimensions = dimensions + 2 * vacuum
    
    # Create a new lattice
    lattice = np.diag(new_dimensions)
    
    # Center the molecule in the new box
    center = (max_coords + min_coords) / 2
    new_center = new_dimensions / 2
    translation = new_center - center
    
    # Translate the molecule
    new_coords = coords + translation
    
    # Create a new structure
    structure = Structure(lattice, molecule.species, new_coords, coords_are_cartesian=True)
    
    return structure

# Function to visualize the structure using py3Dmol
def visualize_structure(structure, html_file_name='viz.html'):
    spin = st.checkbox('Spin', value = False, key='key'+html_file_name)
    view = py3Dmol.view(width=500, height=400)
    cif_for_visualization = structure.to(fmt="cif")
    view.addModel(cif_for_visualization, 'cif')
    # view.setStyle({'stick': {'radius': 0.2}})
    view.setStyle({'sphere':{'colorscheme':'Jmol','scale':0.3},
                    'stick':{'colorscheme':'Jmol', 'radius':0.2}})
    view.addUnitCell()
    view.zoomTo()
    view.spin(spin)
    view.setClickable({'clickable':'true'});
    view.enableContextMenu({'contextMenuEnabled':'true'})
    view.show()
    view.render()
    # view.png()
    t = view.js()
    f = open(html_file_name, 'w')
    f.write(t.startjs)
    f.write(t.endjs)
    f.close()

    HtmlFile = open(html_file_name, 'r', encoding='utf-8')
    source_code = HtmlFile.read() 
    components.html(source_code, height = 300, width=500)
    HtmlFile.close()
    

def search_pubchem(query):
    compounds = pcp.get_compounds(query, namespace='name', as_dataframe=False)
    return compounds


def get_molecule(cid):
    xyz_str = generate_xyz_coordinates(cid)
    return Molecule.from_str(xyz_str, fmt='xyz')


def format_xyz(molecule):
    xyz = molecule.to(fmt="xyz")
    return xyz



st.title("PubChem ➡️ `VASP`")

search_query = st.text_input("Enter molecule name or formula to search in the PubChem database", placeholder='Water / H2O')
molecule_df = None
compounds = None

if not search_query=="":
    compounds = search_pubchem(search_query)
    # st.dataframe(compounds) # if using the pubchempy dataframe

if compounds:
    molecule_df = pd.DataFrame(
        [(compound.cid, compound.iupac_name, compound.molecular_formula, compound.molecular_weight, compound.isomeric_smiles)
            for compound in compounds],
        columns=["CID", "Name", "Formula", "Weight", "Isomeric SMILES"]
    )

    st.success(f"{len(molecule_df)} molecule(s) found!")
    st.dataframe(molecule_df)

if compounds is not None:
    # selected_cid = st.selectbox("Select a molecule", [cid for cid in compounds.index]) # if using the pubchempy dataframe
    selected_cid = st.selectbox("Select a molecule", molecule_df["CID"])
    # st.write(generate_xyz_coordinates(selected_cid))
    selected_molecule = get_molecule(selected_cid)

    st.subheader("3D Atomic Coordinates")
    # Create a dataframe with atomic symbols and atomic coordinates
    st.dataframe(pd.DataFrame({"Atomic Symbol": selected_molecule.species, "X": selected_molecule.cart_coords[:, 0], 
                               "Y": selected_molecule.cart_coords[:, 1], "Z": selected_molecule.cart_coords[:, 2]}))


    # opt_geom = st.checkbox(label= 'Optimize Geometry via ASE EMT calculator (Beta - Does not work well yet)', value=False)
    # if opt_geom:
    #     ase_atoms = AseAtomsAdaptor().get_atoms(selected_molecule)
    #     # calc = SinglePointCalculator(ase_atoms, EMT())
    #     # Set up the calculator for energy and forces using EMT (or other forcefield)
    #     calc = EMT()
    #     ase_atoms.set_calculator(calc)
    #     # Set up the optimizer (BFGS in this example)
    #     optimizer = BFGS(ase_atoms)

    #     # Run the optimization
    #     optimizer.run(fmax=0.05, steps=10)  # Adjust fmax value as needed

    #     # Get the optimized structure as Pymatgen structure
    #     selected_molecule = AseAtomsAdaptor().get_molecule(ase_atoms)
    



    st.download_button(
        "Download XYZ",
        data=format_xyz(selected_molecule),
        file_name="molecule.xyz",
        mime="chemical/x-xyz"
    )

    # # Create a structure with vacuum
    # structure = selected_molecule.get_boxed_structure(15, 15, 15)  # 15 Å vacuum on each side
    vacuum_size = st.slider("Select vacuum size (Å)", min_value=10, max_value=35, value=15, step=1)
    # Create a structure with vacuum
    structure = create_structure_with_vacuum(selected_molecule, vacuum=vacuum_size)  # 15 Å vacuum on each side
    # Visualization
    visualize_structure(structure)

    # Generate VASP input files
    poscar_content, kpoints_content, incar_content = generate_vasp_input_files(structure)

    st.subheader("VASP Input Files")
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