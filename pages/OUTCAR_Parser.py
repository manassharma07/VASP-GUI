import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
from io import StringIO
from ase.io import read, write
from ase.visualize import view
from ase import Atoms
import py3Dmol
import streamlit.components.v1 as components
import base64
import tempfile
from pymatgen.io.ase import AseAtomsAdaptor

# Function to parse energies, structures, and forces from OUTCAR
def parse_outcar(outcar_text):
    # Create a temporary file
    with open('temp_outcar', 'w') as temp_outcar:
        # Write the OUTCAR content to the temporary file
        temp_outcar.write(outcar_text)
        temp_outcar.flush()  # Ensure content is written
    structures = read('temp_outcar', format='vasp-out', index=':')
    energies = []
    forces = []

    lines = outcar_text.split('\n')
    for line in lines:
        if "free  energy   TOTEN" in line:
            energies.append(float(line.split()[-2]))

    for structure in structures:
        forces.append(structure.get_forces())
    
    return structures, energies, forces

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

# Function to visualize the structure using ASE and py3Dmol
def visualize_structure(structure):
    xyz_str = structure.write("xyz", format="xyz")
    xyz_str = xyz_str.replace("\n", "\\n").replace("\r", "")

    view = py3Dmol.view(width=500, height=400)
    view.addModel(xyz_str, "xyz")
    view.setStyle({'sphere': {'colorscheme': 'Jmol', 'scale': 0.3},
                   'stick': {'colorscheme': 'Jmol', 'radius': 0.2}})
    view.addUnitCell()
    view.zoomTo()
    view.spin(False)
    view.show()

    t = view.js()
    components.html(t, height=400, width=500)

# Function to allow CIF download
def cif_download_link(structure):
    cif_str = StringIO()
    write(cif_str, structure, format='cif')
    b64 = base64.b64encode(cif_str.getvalue().encode()).decode()
    href = f'<a href="data:text/cif;base64,{b64}" download="structure.cif">Download CIF</a>'
    return href

st.title("VASP OUTCAR Parser")
st.write("This tool parses a VASP OUTCAR file and extracts relevant information.")

st.write('You can either paste the OUTCAR file contents below or upload the source file')
contents = st.text_area(label='Enter the contents of the OUTCAR file here', value='', placeholder='Put your text here',
                        height=400, key='input_text_area')
# Create a file uploader widget
file = st.file_uploader("or Upload the file")

if file is not None:
    # If a file is uploaded, read its contents
    # contents = file.read()
    # To read file as bytes:
    bytes_data = file.getvalue()

    # To convert to a string based IO:
    stringio = StringIO(file.getvalue().decode("utf-8"))

    # To read file as string:
    contents = stringio.read()
    # st.write(contensts)

if contents != '':
    # contents = file.getvalue().decode("utf-8")
    
    # Parse structures, energies, and forces
    structures, energies, forces = parse_outcar(contents)
    
    # Display the number of structures found
    num_structures = len(structures)
    st.subheader(f"Number of Structures Found: {num_structures}")
    
    if num_structures > 1:
        selected_structure_index = st.selectbox("Select structure to display:", range(num_structures))
    else:
        selected_structure_index = 0
    
    structure = structures[selected_structure_index]
    
    # Display parsed data
    st.subheader("Parsed Energies")
    energy_data = {
        "SCF Iteration": list(range(1, len(energies) + 1)),
        "Total Energy": energies,
    }
    df_energy = pd.DataFrame(energy_data)
    
    st.dataframe(df_energy)

    tab1, tab2 = st.tabs(["Plotly", "Matplotlib"])
    with tab1:
        fig = px.scatter(df_energy, x="SCF Iteration", y="Total Energy", title="Total Energy vs SCF Iteration")
        fig.update_traces(mode='lines+markers', marker={'size': 8})
        st.plotly_chart(fig)
    with tab2:
        plt.figure(figsize=(10, 6))
        plt.plot(df_energy["SCF Iteration"], df_energy["Total Energy"], marker='o', linestyle='-', color='b', label='Total Energy')
        plt.xlabel("SCF Iteration")
        plt.ylabel("Energy (eV)")
        plt.title("Total Energy vs SCF Iteration")
        plt.legend()
        st.pyplot(plt)
    
    # Display structure information and visualization
    display_structure_info(AseAtomsAdaptor.get_structure(structure))
    visualize_structure(structure)
    
    # CIF download link
    st.markdown(cif_download_link(structure), unsafe_allow_html=True)
    
    # Display forces
    st.subheader("Atomic Forces")
    forces_df = pd.DataFrame(forces[selected_structure_index], columns=["Fx", "Fy", "Fz"])
    st.dataframe(forces_df)
