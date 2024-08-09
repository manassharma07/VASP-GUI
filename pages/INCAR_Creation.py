import streamlit as st
import base64

# Initialize default INCAR settings
incar_defaults = {
    'SYSTEM': 'Calculation Name',  # Description of the system
    'ISTART': 0,                              # Job type (0 = new, 1 = orbitals from WAVECAR)
    'ICHARG': 2,                              # Charge type (2 = atom, 1 = from file)
    'ENCUT': 500.00,                          # Plane-wave cutoff energy
    'ALGO': 'Normal',                         # Algorithm for electronic optimization
    'NELM': 60,                               # Max number of electronic steps
    'EDIFF': 1E-06,                           # Stopping criterion for electronic minimization
    'ISMEAR': 0,                              # Smearing type (0 = Gaussian)
    'SIGMA': 0.05,                            # Smearing width
    'EDIFFG': -1E-02,                         # Stopping criterion for ionic relaxation
    'NSW': 20,                                # Number of steps for ionic relaxation
    'IBRION': 2,                              # Ionic relaxation method (2 = conjugate gradient)
    'POTIM': 0.5,                             # Step size for ionic motion
    'ISPIN': 2,                               # Spin polarization (1 = no, 2 = yes)
    'MAGMOM': '2*5.0',                        # Initial magnetic moment per atom
    'LORBIT': 11,                             # DOS and projections output
    'IVDW': 0                                 # van der Waals interaction correction
}

# Helper function to generate the INCAR text
def generate_incar(incar_settings):
    # Calculate the padding required for alignment
    max_key_length = max(len(key) for key in incar_settings.keys())
    max_value_length = max(len(str(value)) for value in incar_settings.values())

    incar_lines = []
    for key, value in incar_settings.items():
        comment = f" ! {key} description"
        if key == 'SYSTEM':
            comment = " ! Name of the calculation"
        elif key == 'ENCUT':
            comment = " ! Energy cutoff for plane-waves (eV)"
        elif key == 'PREC':
            comment = " ! Precision mode: (Low, Medium, Accurate)"
        elif key == 'EDIFF':
            comment = " ! Energy convergence criterion"
        elif key == 'IBRION':
            comment = " ! Algorithm for ionic relaxation"
        elif key == 'NSW':
            comment = " ! Number of ionic steps"
        elif key == 'ISMEAR':
            comment = " ! Smearing type: (-5 to 2, where 0 = Gaussian)"
        elif key == 'SIGMA':
            comment = " ! Smearing width (eV)"
        elif key == 'MAGMOM':
            comment = " ! Initial magnetic moment per atom"
        elif key == 'ISPIN':
            comment = " ! Spin polarization: 1 = no, 2 = yes"
        elif key == 'LORBIT':
            comment = " ! Density of states and projections"
        elif key == 'ALGO':
            comment = " ! Algorithm for electronic optimization"
        elif key == 'NELM':
            comment = " ! Maximum number of electronic steps"
        elif key == 'EDIFFG':
            comment = " ! Stopping criterion for ionic relaxation"
        elif key == 'POTIM':
            comment = " ! Step size for ionic motion"
        elif key == 'ISTART':
            comment = " ! Job start type (0 = new, 1 = from WAVECAR)"
        elif key == 'ICHARG':
            comment = " ! Charge initialization (2 = atom, 1 = from file)"
        elif key == 'IVDW':
            comment = " ! Van der Waals interaction correction"
        
        # Align comments by adding spaces between value and comment
        key_value_str = f"{key} = {value}"
        padding = " " * (max_key_length + max_value_length - len(key_value_str) + 4)  # Adjust padding to align comments
        incar_lines.append(f"{key_value_str}{padding}{comment}")
    
    return "\n".join(incar_lines)


# Helper function to generate the KPOINTS file text
def generate_kpoints(kmesh_type, kpoints, shift):
    if kmesh_type == 'Automatic':
        kpoints_text = "Automatic mesh\n0\nGamma\n"
        kpoints_text += f"{kpoints[0]} {kpoints[1]} {kpoints[2]}\n"
        kpoints_text += f"{shift[0]} {shift[1]} {shift[2]}"
    elif kmesh_type == 'Gamma':
        kpoints_text = "Gamma mesh\n0\nGamma\n"
        kpoints_text += f"{kpoints[0]} {kpoints[1]} {kpoints[2]}\n"
        kpoints_text += f"{shift[0]} {shift[1]} {shift[2]}"
    elif kmesh_type == 'Monkhorst-Pack':
        kpoints_text = "Monkhorst-Pack mesh\n0\nMonkhorst\n"
        kpoints_text += f"{kpoints[0]} {kpoints[1]} {kpoints[2]}\n"
        kpoints_text += f"{shift[0]} {shift[1]} {shift[2]}"
    else:
        kpoints_text = "KPOINTS file could not be generated."
    return kpoints_text

# Function to download files
def download_file(file_name, file_content):
    b64 = base64.b64encode(file_content.encode()).decode()
    href = f'<a href="data:file/txt;base64,{b64}" download="{file_name}">Download {file_name}</a>'
    st.markdown(href, unsafe_allow_html=True)

# Streamlit app layout
st.title("INCAR and KPOINTS File Generator")

# User inputs for INCAR settings
incar_settings = incar_defaults.copy()

incar_settings['SYSTEM'] = st.text_input("System Name (SYSTEM)", incar_defaults['SYSTEM'])

incar_settings['ISTART'] = st.selectbox("Job Start Type (ISTART)", [0, 1], index=[0, 1].index(incar_defaults['ISTART']))

incar_settings['ICHARG'] = st.selectbox("Charge Initialization (ICHARG)", [1, 2, 10], index=[1, 2, 10].index(incar_defaults['ICHARG']))

incar_settings['ENCUT'] = st.slider("Plane-Wave Cutoff Energy (ENCUT)", min_value=200.0, max_value=1000.0, value=incar_defaults['ENCUT'], step=10.0)

incar_settings['ALGO'] = st.selectbox("Electronic Optimization Algorithm (ALGO)", ['Normal', 'Fast', 'All'], index=['Normal', 'Fast', 'All'].index(incar_defaults['ALGO']))

incar_settings['NELM'] = st.number_input("Max Number of Electronic Steps (NELM)", min_value=0, max_value=200, value=incar_defaults['NELM'])

incar_settings['EDIFF'] = st.number_input("Energy Convergence Criterion (EDIFF)", min_value=1e-8, max_value=1e-3, value=incar_defaults['EDIFF'], format="%.1e")

incar_settings['ISMEAR'] = st.selectbox("Smearing Type (ISMEAR)", [-5, -1, 0, 1, 2], index=[-5, -1, 0, 1, 2].index(incar_defaults['ISMEAR']))

incar_settings['SIGMA'] = st.slider("Smearing Width (SIGMA)", min_value=0.01, max_value=0.5, value=incar_defaults['SIGMA'], step=0.01)

incar_settings['EDIFFG'] = st.number_input("Ionic Relaxation Criterion (EDIFFG)", min_value=-1e-1, max_value=-1e-7, value=incar_defaults['EDIFFG'], format="%.1e")

incar_settings['NSW'] = st.number_input("Number of Ionic Steps (NSW)", min_value=0, max_value=200, value=incar_defaults['NSW'])

incar_settings['IBRION'] = st.selectbox("Ionic Relaxation Method (IBRION)", [0, 1, 2, 3], index=[0, 1, 2, 3].index(incar_defaults['IBRION']))

incar_settings['POTIM'] = st.slider("Step Size for Ionic Motion (POTIM)", min_value=0.01, max_value=2.0, value=incar_defaults['POTIM'], step=0.01)

incar_settings['ISPIN'] = st.selectbox("Spin Polarization (ISPIN)", [1, 2], index=[1, 2].index(incar_defaults['ISPIN']))

incar_settings['MAGMOM'] = st.text_input("Initial Magnetic Moment (MAGMOM)", incar_defaults['MAGMOM'])

incar_settings['LORBIT'] = st.selectbox("Density of States (LORBIT)", [0, 10, 11], index=[0, 10, 11].index(incar_defaults['LORBIT']))

incar_settings['IVDW'] = st.selectbox("Van der Waals Correction (IVDW)", [0, 11, 12], index=[0, 11, 12].index(incar_defaults['IVDW']))

# Display the dynamically generated INCAR file
st.subheader("Generated INCAR File")
incar_text = generate_incar(incar_settings)
incar_display = st.text_area("INCAR", incar_text, height=300)

# KPOINTS settings
st.subheader("KPOINTS Settings")
kmesh_type = st.selectbox("K-Mesh Type", ['Automatic', 'Gamma', 'Monkhorst-Pack'])

kpoints = st.columns(3)
kpoint_values = [
    kpoints[0].number_input("Kx", min_value=1, max_value=10, value=3),
    kpoints[1].number_input("Ky", min_value=1, max_value=10, value=3),
    kpoints[2].number_input("Kz", min_value=1, max_value=10, value=3)
]

shift = st.columns(3)
shift_values = [
    shift[0].number_input("Shift x", min_value=0.0, max_value=1.0, value=0.0, step=0.1),
    shift[1].number_input("Shift y", min_value=0.0, max_value=1.0, value=0.0, step=0.1),
    shift[2].number_input("Shift z", min_value=0.0, max_value=1.0, value=0.0, step=0.1)
]

# Display the dynamically generated KPOINTS file
st.subheader("Generated KPOINTS File")
kpoints_text = generate_kpoints(kmesh_type, kpoint_values, shift_values)
kpoints_display = st.text_area("KPOINTS", kpoints_text, height=200)

# Download buttons
st.subheader("Download Files")
download_file("INCAR", incar_text)
download_file("KPOINTS", kpoints_text)

# Suggestions for additional INCAR options
st.subheader("Suggestions for Additional INCAR Options")
st.markdown("""
- **NELM**: Maximum number of electronic steps (Default: 60)
- **ALGO**: Algorithm for electronic minimization (e.g., Fast, Normal, All)
- **NCORE**: Number of cores to be used by VASP (parallelization)
- **LMAXMIX**: Maximum l-quantum number for non-spherical contributions (useful for d and f electrons)
- **LDAU**: DFT+U corrections for strongly correlated systems (requires additional settings)
- **GGA**: Type of exchange-correlation functional (e.g., PE, RPBE)
- **LWAVE**: Write WAVECAR file (default is .TRUE.)
- **LCHARG**: Write CHGCAR file (default is .TRUE.)
""")