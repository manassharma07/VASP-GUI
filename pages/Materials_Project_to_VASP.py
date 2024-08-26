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
from ase.io.vasp import write_vasp

# Set page config
st.set_page_config(page_title='Materials Project ➡️ VASP', layout='wide', page_icon="⚛️",
menu_items={
         'About': "A web app to help you create VASP input files from materialsproject entries."
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

# Set your Materials Project API key
api_key = st.secrets["MP_API"]

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
    poscar_content = convert_to_poscar_ase(structure, 'POSCAR', direct=False)

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

# Function to convert atomic coordinates to Bohr units
def convert_to_bohr(structure):
    coords = [(site.coords[0], site.coords[1], site.coords[2], site.species_string) for site in structure.sites]
    return [(x * 1.88972612456506, y * 1.88972612456506, z * 1.88972612456506, element.lower()) for x, y, z, element in coords]

# Function to generate coordinate text
def generate_coord_text(coords_bohr):
    coord_text = "$coord\n"
    for coord in coords_bohr:
        coord_text += f"    {coord[0]:.8f}   {coord[1]:.8f}   {coord[2]:.8f}    {coord[3]}\n"
    coord_text += "$end"
    return coord_text

# Function to generate lattice parameter text
def generate_lattice_text(structure):
    lattice_params = structure.lattice.abc
    angles = structure.lattice.angles
    lattice_text = "$cell angs\n"
    lattice_text += f"  {lattice_params[0]:.8f}   {lattice_params[1]:.8f}   {lattice_params[2]:.8f}   {angles[0]}   {angles[1]}   {angles[2]}\n"
    lattice_text += "$periodic 3\n"
    lattice_text += "$kpoints\n"
    lattice_text += "    nkpoints 1 1 1 # Gamma point calculation"
    return lattice_text

# Function to convert a structure to CIF
def convert_to_cif(structure, filename):
    cif_writer = CifWriter(structure)
    cif_writer.write_file(filename)

def parse_cif_ase(stringio):
    # Read CIF
    atoms = read(stringio, format="cif")

    # Convert ASE Atoms to pymatgen Structure
    structure = AseAtomsAdaptor().get_structure(atoms)
    return structure

def convert_pymatgen_to_ase_to_pymatgen(structure):
    convert_to_cif(structure, "temp.cif")
    file = open("temp.cif", 'r')
    return parse_cif_ase(file)


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

# return filecontents
def read_file(filename):
    with open(filename, 'r') as file:
        return file.read()

# Function to display structure information
def display_structure_info(structure):
    st.subheader("Structure Information")
    st.write("Formula: ", structure.composition.reduced_formula)

    #Display lattice parameters
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
        # st.write("Lattice Vectors:")
        st.table(df_vectors)

    # Create a list of atomic coordinates
    with st.expander("Atomic Coordinates", expanded=False):
        coord_type = st.selectbox('Coordinate type', ['Cartesian', 'Fractional/Crystal'])
        if coord_type=='Cartesian':
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
        # st.write("Atomic Coordinates:")
        st.table(df_coords)



# Create an instance of MPRester
mpr = MPRester(api_key)

# Streamlit app
docs = None
st.write('# Materials Project ➡️ VASP')
st.write("#### Get POSCAR, KPOINTS, and INCAR files from Materials Project Database")

# Search by material_id, formula or elements?
# input_type = st.selectbox("Search by:", ['Formula','Material ID','Elements'])
input_type = st.selectbox("Search by:", ['Formula','Material ID'])


# Search for materials
if input_type=='Formula':
    formula = st.text_input("Enter formula:", placeholder='NaCl')
elif input_type=='Material ID':
    material_id = st.text_input("Enter material id:", placeholder='mp-22862')
elif input_type=='Elements':
    _elements = st.text_input("Enter elements:", placeholder='Na, Cl')
    elements = [item.strip() for item in _elements.split(',')]
if input_type=='Formula':
    if not formula=="":
        with st.spinner("Searching..."):
            docs = mpr.summary.search(formula=[formula], fields=["structure", "band_gap", "material_id", "is_stable", "is_metal", "symmetry", "formula_pretty"])
elif input_type=='Material ID':
    if not material_id=="":
        with st.spinner("Searching..."):
            docs = mpr.summary.search(material_ids=[material_id], fields=["structure", "band_gap", "material_id", "is_stable", "is_metal", "symmetry", "formula_pretty"])
elif input_type=='Elements':
    if not _elements=="":
        with st.spinner("Searching..."):
            docs = mpr.summary.search(elements=elements, fields=["structure", "band_gap", "material_id", "is_stable", "is_metal", "symmetry", "formula_pretty"])

if docs is not None:
    if len(docs) > 0:
        st.success(f"Matching materials found: {len(docs)}")

        # Display materials as a table
        table_data = []
        for doc in docs:
            table_data.append({
                "Material ID": doc.material_id,
                "Band Gap": doc.band_gap,
                "Crystal System": str(doc.symmetry.crystal_system),
                "Symbol": str(doc.symmetry.symbol),
                "Group": str(doc.symmetry.point_group),
                "Formula": doc.formula_pretty,
                "Is Stable": doc.is_stable,
                "Is Metal": doc.is_metal
            })
        st.table(table_data)
    else:
        st.warning("No materials found for the given formula.")

if docs is not None:

    # Select a material
    selected_material = st.selectbox("Select a material:", [doc.material_id for doc in docs])
    selected_doc = next((doc for doc in docs if doc.material_id == selected_material), None)
    structure = selected_doc.structure
    # Get conventional structure
    analyzer = SpacegroupAnalyzer(structure)
    conventional_structure = analyzer.get_conventional_standard_structure()
    # The structure returned by pymatgen doesn't necessarily have the lattice vector a parallel to the x-axis. 
    # See for example Ru (Hexagonal) mp-33
    # A trick could be to convert it to CIF and then use ASE to read that CIF (as ASE seems to be following the convention)
    # Then convert the ASE structure to pymatgen
    conventional_structure = convert_pymatgen_to_ase_to_pymatgen(conventional_structure)


    # Get primitive structure
    primitive_structure = analyzer.get_primitive_standard_structure()
    # The structure returned by pymatgen doesn't necessarily have the lattice vector a parallel to the x-axis. 
    # See for example Ru (Hexagonal) mp-33
    # A trick could be to convert it to CIF and then use ASE to read that CIF (as ASE seems to be following the convention)
    # Then convert the ASE structure to pymatgen
    primitive_structure = convert_pymatgen_to_ase_to_pymatgen(primitive_structure)

    # Choose between primitive or conventional
    selected_structure_type = st.selectbox("Unit cell type:", ['Primitive Cell','Conventional Unit Cell'])

    # Display structure information
    if selected_structure_type=='Primitive Cell':
        display_structure_info(primitive_structure)
    else:
        display_structure_info(conventional_structure)

    # Visualize the structure
    if selected_structure_type=='Primitive Cell':
        visualize_structure(primitive_structure, "viz1.html")
    else:
        visualize_structure(conventional_structure, "viz1.html")

    # Download CIF files
    st.subheader("Download CIF Files")
    col1, col2 = st.columns(2)
    if primitive_structure is not None:
        convert_to_cif(primitive_structure, "primitive_unit_cell.cif")
        col1.download_button('Download Primitive Unit Cell CIF', data=read_file("primitive_unit_cell.cif"), file_name='primitive_unit_cell.cif', key='primitive_cif_button')

    if conventional_structure is not None:
        convert_to_cif(conventional_structure, "conventional_unit_cell.cif")
        col2.download_button('Download Conventional Unit Cell CIF', data=read_file("conventional_unit_cell.cif"), file_name='conventional_unit_cell.cif', key='conventional_cif_button')

    # Get VASP POSCAR KPOINTS and INCAR file contents
    st.subheader("VASP Files")
    # Add a selection box for coordinate type
    coord_type = st.selectbox("Select coordinate type for POSCAR", ["Direct", "Cartesian"])
    if coord_type=='Direct':
        direct = True
    else:
        direct = False
    if selected_structure_type=='Primitive Cell':
        poscar_content, kpoints_content, incar_content = generate_vasp_input_files(primitive_structure, direct=direct)
    else:
        poscar_content, kpoints_content, incar_content = generate_vasp_input_files(conventional_structure, direct=direct)

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


    # Create supercells
    # with st.expander("Model Supercell", expanded=False):
    st.subheader('Model Supercell')
    # Create three columns for inputs
    col1, col2, col3 = st.columns(3)

    # Add input fields for nx, ny, and nz in separate columns
    nx = col1.number_input("Enter nx", min_value=1, step=1)
    ny = col2.number_input("Enter ny", min_value=1, step=1)
    nz = col3.number_input("Enter nz", min_value=1, step=1)
    # st.write(primitive_structure)
    # st.write(conventional_structure)
    if selected_structure_type=='Primitive Cell':
        supercell_structure = primitive_structure.copy()
        supercell_structure.make_supercell([int(nx), int(ny), int(nz)])
    else:
        supercell_structure = conventional_structure.copy()
        supercell_structure.make_supercell([int(nx), int(ny), int(nz)])
    # Get the number of atoms
    num_atoms_supercell = supercell_structure.num_sites
    if num_atoms_supercell<500:
        visualize_structure(supercell_structure, 'viz2.html')
    else:
        st.warning("We can't visualize your supercell as it contains more than 500 atoms which is a bit too much for a free web app.\n But don't worry, RIPER can still do the calculations with ease (provided you have the required resources).")

    # Get VASP POSCAR KPOINTS and INCAR file contents
    st.subheader("VASP Files for the Supercell")
    poscar_content, kpoints_content, incar_content = generate_vasp_input_files(supercell_structure)

    # Display POSCAR and KPOINTS in editable text boxes
    col1, col2 = st.columns(2)

    with col1:
        st.write("### POSCAR")
        poscar_editable2 = st.text_area("POSCAR Content", poscar_content, height=300)
        st.download_button(
            label="Download POSCAR",
            data=poscar_editable2,
            file_name='POSCAR',
            mime='text/plain',
        )

    with col2:
        st.write("### KPOINTS")
        kpoints_editable2 = st.text_area("KPOINTS Content", kpoints_content, height=300)
        st.download_button(
            label="Download KPOINTS",
            data=kpoints_editable2,
            file_name='KPOINTS',
            mime='text/plain',
        )

    # Display INCAR file
    st.subheader("Sample INCAR")
    incar_editable2 = st.text_area("INCAR Content", incar_content, height=300)
    st.download_button(
        label="Download INCAR",
        data=incar_editable2,
        file_name='INCAR',
        mime='text/plain',
    )




# import streamlit as st
# from mp_api.client import MPRester
# from pymatgen.io.cif import CifWriter
# from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
# import py3Dmol
# import pandas as pd
# import streamlit.components.v1 as components
# from ase.io import read
# from pymatgen.io.ase import AseAtomsAdaptor
# from io import StringIO
# from pymatgen.io.vasp.inputs import Poscar
# from pymatgen.core import Structure

# # Set page config
# st.set_page_config(page_title='Materials Project ➡️ VASP', layout='wide', page_icon="⚛️",
# menu_items={
#          'About': "A web app to help you create VASP input files from materialsproject entries."
#      })

# # Sidebar stuff
# st.sidebar.write('# About')
# st.sidebar.write(' The web app is made By [Manas Sharma](https://manas.bragitoff.com)')
# st.sidebar.write('### *Powered by*')
# st.sidebar.write('* [Py3Dmol](https://3dmol.csb.pitt.edu/) for Chemical System Visualizations')
# st.sidebar.write('* [Streamlit](https://streamlit.io/) for making of the Web App')
# st.sidebar.write('* [PyMatgen](https://pymatgen.org/) for Periodic Structure Representations')
# st.sidebar.write('* [PubChempy](https://pypi.org/project/PubChemPy/1.0/) for Accessing the PubChem Database')
# st.sidebar.write('* [MP-API](https://pypi.org/project/mp-api/) for Accessing the Materials Project Database')
# st.sidebar.write('* [ASE](https://wiki.fysik.dtu.dk/ase/) for File Format Conversions')
# st.sidebar.write('### *Useful links*')
# st.sidebar.write('[Web App Source Code](https://github.com/manassharma07/VASP-GUI)')
# st.sidebar.write('[VASP Wiki](https://www.vasp.at/wiki/index.php/The_VASP_Manual)')
# st.sidebar.write('[VASP Official Website](https://www.vasp.at/)')
# st.sidebar.write('[VASP Forum](https://www.vasp.at/forum/)')

# # Set your Materials Project API key
# api_key = st.secrets["MP_API"]



# def parse_cif_ase(stringio):
#     # Read CIF
#     atoms = read(stringio, format="cif")

#     # Convert ASE Atoms to pymatgen Structure
#     structure = AseAtomsAdaptor().get_structure(atoms)
#     return structure

# # Function to visualize the structure using py3Dmol
# def visualize_structure(structure, html_file_name='viz.html'):
#     spin = st.checkbox('Spin', value=False, key='key'+html_file_name)
#     view = py3Dmol.view(width=500, height=400)
#     cif_for_visualization = structure.to(fmt="cif")
#     view.addModel(cif_for_visualization, 'cif')
#     view.setStyle({'sphere': {'colorscheme': 'Jmol', 'scale': 0.3},
#                    'stick': {'colorscheme': 'Jmol', 'radius': 0.2}})
#     view.addUnitCell()
#     view.zoomTo()
#     view.spin(spin)
#     view.setClickable({'clickable':'true'})
#     view.enableContextMenu({'contextMenuEnabled':'true'})
#     view.show()
#     view.render()
#     t = view.js()
#     f = open(html_file_name, 'w')
#     f.write(t.startjs)
#     f.write(t.endjs)
#     f.close()

#     HtmlFile = open(html_file_name, 'r', encoding='utf-8')
#     source_code = HtmlFile.read() 
#     components.html(source_code, height=300, width=500)

# # Streamlit app
# st.write('# Materials Project ➡️ VASP')
# st.write("#### Generate VASP input files from structures obtained from the Materials Project.")

# # Select Material ID
# material_id = st.text_input("Enter the Material ID from Materials Project", value="mp-149")

# # Fetch the structure from Materials Project
# if material_id:
#     with MPRester(api_key) as mpr:
#         structure = mpr.get_structure_by_material_id(material_id)
    
#     if structure:
#         st.subheader("Structure Information")
#         st.write("Formula: ", structure.composition.reduced_formula)
#         visualize_structure(structure, 'structure_visualization.html')

#         poscar_content, kpoints_content, incar_content = generate_vasp_input_files(structure)

#         # Display POSCAR and KPOINTS in editable text boxes
#         col1, col2 = st.columns(2)

#         with col1:
#             st.write("### POSCAR")
#             poscar_editable = st.text_area("POSCAR Content", poscar_content, height=300)
#             st.download_button(
#                 label="Download POSCAR",
#                 data=poscar_editable,
#                 file_name='POSCAR',
#                 mime='text/plain',
#             )

#         with col2:
#             st.write("### KPOINTS")
#             kpoints_editable = st.text_area("KPOINTS Content", kpoints_content, height=300)
#             st.download_button(
#                 label="Download KPOINTS",
#                 data=kpoints_editable,
#                 file_name='KPOINTS',
#                 mime='text/plain',
#             )

#         # Display INCAR file
#         st.subheader("Sample INCAR")
#         incar_editable = st.text_area("INCAR Content", incar_content, height=300)
#         st.download_button(
#             label="Download INCAR",
#             data=incar_editable,
#             file_name='INCAR',
#             mime='text/plain',
#         )
#     else:
#         st.warning("Could not retrieve structure. Please check the Material ID.")
# else:
#     st.warning("Please enter a Material ID.")

