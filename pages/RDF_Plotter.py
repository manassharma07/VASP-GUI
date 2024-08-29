import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from ase import Atoms
from ase.io import read
from scipy.spatial.distance import pdist
from io import StringIO
import io

# Set page config
st.set_page_config(page_title='RDF Plot from Trajectory', layout='wide', page_icon="⚛️",
menu_items={
         'About': "A web app to help you plot the RDF from XYZ trjectory file."
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

def calculate_rdf(trajectory, r_range, dr, elements=None):
    rdf = np.zeros(len(r_range) - 1)
    n_frames = len(trajectory)
    
    for atoms in trajectory:
        if elements:
            indices = [atom.index for atom in atoms if atom.symbol in elements]
            positions = atoms.positions[indices]
        else:
            positions = atoms.positions
        
        dists = pdist(positions)
        hist, _ = np.histogram(dists, bins=r_range)
        rdf += hist
    
    rdf = rdf / n_frames
    
    # Normalize RDF
    volume = atoms.get_volume()
    n_atoms = len(positions)
    norm = 4 * np.pi * r_range[1:]**2 * dr * n_atoms * (n_atoms - 1) / (2 * volume)
    rdf = rdf / norm
    
    return rdf

st.title("Radial Distribution Function Plotter")

# Input method selection
input_method = st.radio("Choose input method:", ("Upload file", "Paste content"))

if input_method == "Upload file":
    uploaded_file = st.file_uploader("Upload trajectory file (XYZ format)", type="xyz")
    if uploaded_file is not None:
        # trajectory = read(uploaded_file, index=":")
        # Read the file contents
        file_contents = uploaded_file.read().decode("utf-8")
    
        # Create a StringIO object
        string_io = io.StringIO(file_contents)
    
        # Read the trajectory
        trajectory = read(string_io, index=":", format="extxyz")
else:
    trajectory_content = st.text_area("Paste trajectory content (XYZ format)")
    # if trajectory_content:
    #     trajectory = read(StringIO(trajectory_content), format="xyz", index=":")

# Parameters
r_max = st.slider("Maximum distance (r_max)", min_value=1.0, max_value=10.0, value=5.0, step=0.1)
dr = st.slider("Distance resolution (dr)", min_value=0.01, max_value=0.5, value=0.1, step=0.01)

# Plot settings
plot_title = st.text_input("Plot title", "Radial Distribution Function")
x_axis_title = st.text_input("X-axis title", "r (Å)")
y_axis_title = st.text_input("Y-axis title", "g(r)")

# Plot type selection
plot_type = st.radio("Choose plot type:", ("Matplotlib", "Plotly"))

if 'trajectory' in locals():
    r_range = np.arange(0, r_max + dr, dr)
    rdf = calculate_rdf(trajectory, r_range, dr)

    if plot_type == "Matplotlib":
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(r_range[1:], rdf)
        ax.set_xlabel(x_axis_title)
        ax.set_ylabel(y_axis_title)
        ax.set_title(plot_title)
        ax.grid(True)
        st.pyplot(fig)
    else:
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=r_range[1:], y=rdf, mode='lines'))
        fig.update_layout(
            title=plot_title,
            xaxis_title=x_axis_title,
            yaxis_title=y_axis_title,
        )
        st.plotly_chart(fig)

    # Download button for calculated RDF data
    rdf_data = np.column_stack((r_range[1:], rdf))
    csv_buffer = StringIO()
    np.savetxt(csv_buffer, rdf_data, delimiter=",", header="r,g(r)", comments="")
    csv_str = csv_buffer.getvalue()
    
    st.download_button(
        label="Download RDF data as CSV",
        data=csv_str,
        file_name="rdf_data.csv",
        mime="text/csv",
    )
else:
    st.write("Please upload a trajectory file or paste trajectory content to calculate the RDF.")