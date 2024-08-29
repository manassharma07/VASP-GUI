import streamlit as st
import matplotlib.pyplot as plt
import plotly.express as px
import re
import pandas as pd

# Set page config
st.set_page_config(page_title='NVT MD OSZICAR Analyzer', layout='wide', page_icon="⚛️",
menu_items={
         'About': "A web app to help you analyze your MD runs uing VASP."
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

# Function to parse the OSZICAR file and extract relevant data
def parse_oszicar(file_content):
    temperatures = []
    energies = []
    free_energies = []
    e0s = []
    eks = []
    sps = []
    sks = []

    # Regular expression pattern to match the lines containing T, E, F, E0, EK, SP, SK
    pattern = re.compile(r"T=\s*([\d\.\-]+)\s*E=\s*([\d\.\-E\+]+)\s*F=\s*([\d\.\-E\+]+)\s*E0=\s*([\d\.\-E\+]+)\s*EK=\s*([\d\.\-E\+]+)\s*SP=\s*([\d\.\-E\+]+)\s*SK=\s*([\d\.\-E\+]+)")

    for line in file_content:
        match = pattern.search(line)
        if match:
            temperatures.append(float(match.group(1)))
            energies.append(float(match.group(2)))
            free_energies.append(float(match.group(3)))
            e0s.append(float(match.group(4)))
            eks.append(float(match.group(5)))
            sps.append(float(match.group(6)))
            sks.append(float(match.group(7)))

    return temperatures, energies, free_energies, e0s, eks, sps, sks

# Streamlit app
st.title("NVT MD OSZICAR File Analyser")

uploaded_file = st.file_uploader("Upload OSZICAR file")

if uploaded_file is not None:
    # Decode the uploaded file content
    file_content = uploaded_file.read().decode("utf-8").splitlines()

    # Parse the file to extract data
    temperatures, energies, free_energies, e0s, eks, sps, sks = parse_oszicar(file_content)

    # Convert the extracted data to a DataFrame
    data = {
        "Temperature (T)": temperatures,
        "Total Energy (E)": energies,
        "Free Energy (F)": free_energies,
        "E0": e0s,
        "Kinetic Energy (EK)": eks,
        "SP": sps,
        "SK": sks
    }
    df = pd.DataFrame(data)

    # Display the DataFrame as a table
    st.subheader("Extracted Data Table")
    st.dataframe(df)

    # Allow the user to download the extracted data as CSV
    csv = df.to_csv(index=False)
    st.download_button(
        label="Download Extracted Data as CSV",
        data=csv,
        file_name='oszicar_data.csv',
        mime='text/csv',
    )

    # Plot the data
    st.subheader("Plots")
    
    # User selects the plotting library
    plot_type = st.radio(
        "Choose Plotting Library",
        ('matplotlib', 'plotly')
    )

    if plot_type == 'matplotlib':
        # Creating subplots
        fig, axs = plt.subplots(4, 2, figsize=(15, 20))
        axs = axs.flatten()

        axs[0].plot(temperatures, 'b-o')
        axs[0].set_title("Temperature (T)")

        axs[1].plot(energies, 'r-o')
        axs[1].set_title("Total Energy (E)")

        axs[2].plot(free_energies, 'g-o')
        axs[2].set_title("Free Energy (F)")

        axs[3].plot(e0s, 'c-o')
        axs[3].set_title("E0")

        axs[4].plot(eks, 'm-o')
        axs[4].set_title("Kinetic Energy (EK)")

        axs[5].plot(sps, 'y-o')
        axs[5].set_title("SP")

        axs[6].plot(sks, 'k-o')
        axs[6].set_title("SK")

        for ax in axs:
            ax.set_xlabel("Step")
            ax.set_ylabel("Value")
            ax.grid(True)

        fig.tight_layout()
        st.pyplot(fig)

    elif plot_type == 'plotly':
        fig = px.line(df, y=["Temperature (T)", "Total Energy (E)", "Free Energy (F)", 
                             "E0", "Kinetic Energy (EK)", "SP", "SK"], 
                      title="OSZICAR Data Over Time",
                      labels={"value": "Value", "index": "Step"})
        fig.update_layout(height=800)
        st.plotly_chart(fig)

    # Explanation of the terms
    st.subheader("Explanation of the Terms")
    explanation_data = {
        "Quantity": ["T", "E", "F", "E0", "EK", "SP", "SK"],
        "Explanation": [
            "The instantaneous temperature.",
            "The total energy including the potential energy F of the ionic degree of freedom, the potential energy SP and kinetic energy SK of the Nose Hoover thermostat, and the kinetic energy of the ionic motion EK. It is called ETOTAL in the OUTCAR file.",
            "The total free energy of the DFT calculation considering the artificial electronic temperature introduced by the smearing factor SIGMA. From the viewpoint of MD, this is the potential energy of the ionic degree of freedom. It is called TOTEN in the OUTCAR file.",
            "The internal energy of the ions.",
            "The kinetic energy of the ionic motion.",
            "The potential energy of the Nose Hoover thermostat.",
            "The kinetic energy of the Nose Hoover thermostat."
        ]
    }
    explanation_df = pd.DataFrame(explanation_data)
    st.table(explanation_df)
