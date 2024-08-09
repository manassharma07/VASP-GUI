import streamlit as st
import pandas as pd
from scipy import signal
import numpy as np
import io
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import re

# Set page config
st.set_page_config(page_title='GW/BSE Spectrum Plotter 📈', layout='wide', page_icon="⚛️",
menu_items={
         'About': "A web app to help you plot the GW/BSE spectrum using the dielectric function from VASP."
     })

# Sidebar stuff
st.sidebar.write('# About')
st.sidebar.write('Originally Made By [Your Name](https://your.website)')
st.sidebar.write('### *Powered by*')
st.sidebar.write('* [Streamlit](https://streamlit.io/) for making of the Web App')
st.sidebar.write('### *Source Code*')
st.sidebar.write('[GitHub Repository](https://github.com/your-repo-link)')

def parse_vasprun(content):
    energy, imaginary, real = [], [], []
    inside_dielectricfunction = False
    temp_real = []

    for line in content.splitlines():
        if "<dielectricfunction>" in line:
            inside_dielectricfunction = True
            continue

        if "</dielectricfunction>" in line:
            inside_dielectricfunction = False
            continue

        if inside_dielectricfunction and line.strip().startswith("<r>"):
            parts = re.split(r'\s+', line.strip())
            
            # Ensure there are at least 5 parts to avoid errors
            if len(parts) >= 5:
                energy.append(float(parts[1]))
                imaginary.append(float(parts[2]))
                temp_real.append(float(parts[3]))

    # The second half of temp_real contains the actual real part values
    n = len(energy) // 2
    real = temp_real[n:]  # Second half corresponds to the real part
    imaginary = imaginary[:n]  # First half corresponds to the imaginary part
    energy = energy[:n]

    df = pd.DataFrame({
        'Energy (eV)': energy,
        'Imaginary Part': imaginary,
        'Real Part': real
    })

    return df


st.title('GW/BSE Spectrum Plotter')
st.write('This utility helps you plot the GW/BSE spectrum by parsing the dielectric function data from `vasprun.xml` generated by VASP.')

st.sidebar.write('### Input Data')
uploaded_file = st.file_uploader('Upload the `vasprun.xml` file')
pasted_content = st.text_area('or paste the contents of the `vasprun.xml` file here:', height=400)

content = None

if uploaded_file is not None:
    content = uploaded_file.read().decode('utf-8')
elif pasted_content:
    content = pasted_content

if content:
    st.write('### Parsed data')
    df = parse_vasprun(content)
    st.dataframe(df)

    output = io.BytesIO()
    writer = pd.ExcelWriter(output, engine='xlsxwriter')
    df.to_excel(writer, sheet_name='Sheet1', index=False)
    writer.save()
    output.seek(0)
    st.download_button('Download data as Excel file', data=output, file_name='gw_bse_spectrum.xlsx')


    ## Plot the dielectric function
    st.write('### Plotting the GW/BSE Spectrum')

    plot_part = st.radio('Select part to plot', ('Imaginary Part', 'Real Part'), index=0)
    x_data = df['Energy (eV)']
    y_data = df['Imaginary Part'] if plot_part == 'Imaginary Part' else df['Real Part']

    if plot_part == 'Imaginary Part':
        # y_data = b
        y_label = 'Imaginary Dielectric Function '
    else:
        # y_data = c
        y_label = 'Real Dielectric Function '

    peaks_col1, peaks_col2 = st.columns(2)
    isFP = peaks_col1.checkbox('Find Peaks', value=False)
    prominence_val = peaks_col2.text_input(label='Prominence value', value='0', key='prominence')
    prominence_val = float(prominence_val)

    if isFP:
        peak_indices, dict_peak = signal.find_peaks(y_data, prominence=prominence_val)
        highest_peak_val = np.max(y_data[peak_indices])

    normalize_col1, normalize_col2 = st.columns(2)
    isNormalize = normalize_col1.checkbox('Normalize', value=True)
    norm_factor = normalize_col2.text_input(label='Value by which to normalize the highest peak', value='1', key='norm_factor')
    norm_factor = float(norm_factor)

    if isNormalize:
        highest_peak_val = np.max(y_data)
        if isFP:
            peak_indices, dict_peak = signal.find_peaks(y_data, prominence=prominence_val)
            highest_peak_val = np.max(y_data[peak_indices])
        y_data = (y_data / highest_peak_val) * norm_factor

    lim_col1, lim_col2, lim_col3, lim_col4 = st.columns(4)
    lxlim = lim_col1.text_input(label='Enter the lower limit for x-axis', value=str(min(x_data)), key='lxlim')
    lxlim = float(lxlim)
    uxlim = lim_col2.text_input(label='Enter the upper limit for x-axis', value=str(max(x_data)), key='uxlim')
    uxlim = float(uxlim)
    lylim = lim_col3.text_input(label='Enter the lower limit for y-axis', value=str(min(y_data) - 0.1 * max(y_data)), key='lylim')
    lylim = float(lylim)
    uylim = lim_col4.text_input(label='Enter the upper limit for y-axis', value=str(max(y_data) + 0.1 * max(y_data)), key='uylim')
    uylim = float(uylim)

    title_col1, title_col2, title_col3 = st.columns(3)
    xtitle = title_col1.text_input(label='Enter the label for x-axis', value='Energy (eV)', key='xtitle')
    ytitle = title_col2.text_input(label='Enter the label for y-axis', value=y_label, key='ytitle')
    plot_title = title_col3.text_input(label='Enter the title for the plot', value=f'GW/BSE Spectrum', key='plot_title')

    figsize_col1, figsize_col2 = st.columns(2)
    figsize_width = figsize_col1.text_input(label='Enter the plot figure width', value='10', key='figsize_width')
    figsize_width = float(figsize_width)
    figsize_height = figsize_col2.text_input(label='Enter the plot figure height', value='6', key='figsize_height')
    figsize_height = float(figsize_height)

    plot_type = st.selectbox('Select the plot style', ('outline', 'shaded', 'shaded+outline'))

    plot_col1, plot_col2, plot_col3, plot_col4 = st.columns(4)
    plot_color = plot_col1.color_picker('Pick a Color for the Spectrum plot', '#0023F9')
    transparency = plot_col2.checkbox('Make the plot transparent', value=True)
    line_width = plot_col3.slider('Line Width', 0.4, 10., 3.0, step=0.1)
    line_style = plot_col4.selectbox('Select the line style', ('solid', 'dashed', 'dotted', 'dashdot'))

    fig, ax = plt.subplots(figsize=[figsize_width, figsize_height])
    ax.tick_params(axis='both', which='major', labelsize=24)
    ax.tick_params(axis='both', which='minor', labelsize=24)

    if plot_type == 'outline':
        ax.plot(x_data, y_data, color=plot_color, lw=line_width, linestyle=line_style)
    if plot_type == 'shaded':
        ax.fill_between(x_data, lylim, y_data, facecolor=plot_color, alpha=0.5)
    if plot_type == 'shaded+outline':
        ax.plot(x_data, y_data, color=plot_color, lw=line_width, linestyle=line_style)
        ax.fill_between(x_data, lylim, y_data, facecolor=plot_color, alpha=0.5)

    ax.set_title(plot_title, fontsize=27)
    ax.set_xlabel(xtitle, fontsize=24)
    ax.set_ylabel(ytitle, fontsize=24)
    ax.set_xlim([lxlim, uxlim])
    ax.set_ylim([lylim, uylim])

    if isFP:
        for i in peak_indices:
            plt.annotate(f'({np.round(x_data[i], 4)}, {np.round(y_data[i], 4)})', [x_data[i], y_data[i]])

    plt.tight_layout()
    plt.savefig('spectrum.png', transparent=transparency)
    st.pyplot(fig)

    with open('spectrum.png', 'rb') as f:
        st.download_button('Download plot as PNG file', f, file_name='spectrum.png')
