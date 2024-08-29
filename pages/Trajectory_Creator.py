import streamlit as st
import tempfile
import os
from ase.io import read, write
import io

# Set page config
st.set_page_config(page_title='Create Trajectory (extxyz) file from OUTCAR', layout='wide', page_icon="⚛️",
menu_items={
         'About': "A web app to help you convert the VASP OUTCAR file to a trajectory file in extxyz format."
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

def convert_outcar_to_extxyz(outcar_content):
    # Create a temporary file to store the OUTCAR content
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.OUTCAR', delete=False) as temp_outcar:
        temp_outcar.write(outcar_content)
        temp_outcar_path = temp_outcar.name

    # Read the OUTCAR file
    input_data = read(temp_outcar_path, index=':', format='vasp-out')

    # Create a StringIO object to store the extxyz content
    extxyz_buffer = io.StringIO()

    # Write the data to the StringIO object in extxyz format
    write(extxyz_buffer, input_data, format='extxyz')

    # Get the content as a string
    extxyz_content = extxyz_buffer.getvalue()

    # Clean up the temporary file
    os.unlink(temp_outcar_path)

    return extxyz_content

st.title("OUTCAR to `trajectory` Converter")

uploaded_file = st.file_uploader("Upload OUTCAR file")

if uploaded_file is not None:
    outcar_content = uploaded_file.getvalue().decode("utf-8")
    
    if st.button("Convert to extxyz"):
        try:
            extxyz_content = convert_outcar_to_extxyz(outcar_content)
            st.success("Conversion successful!")
            
            # Provide download link
            st.download_button(
                label="Download extxyz file",
                data=extxyz_content,
                file_name="trajectory.xyz",
                mime="chemical/x-xyz"
            )
        except Exception as e:
            st.error(f"An error occurred during conversion: {str(e)}")
else:
    st.write("Please upload an OUTCAR file to begin.")