import pandas as pd
import os
import plotly.express as px
import streamlit as st

import helper

file_path = "data/IF_data/"

patients = []
for filename in os.listdir(file_path + "IF1"):
    patients.append(filename[:4])

# TODO
# opis danych - skład komórek procentowy
if 'display_result' not in st.session_state:
    st.session_state.display_result = False
if 'reset' not in st.session_state:
    st.session_state.reset = False


def hide_callback():
    st.session_state.display_result = False
    st.session_state.reset = False


st.title('Tumor spatial data analysis')
patient = st.sidebar.selectbox('Select a patient', patients)

data_IF1 = helper.get_panel('IF1', patient)

st.markdown(
    "Main goal of this project is to locate Tertiary lymphoid structures (TLS) within analyzed tissue. "
    "The analysis was conducted based on data obtained from a biopsy of cancerous tissue. "
    "This tissue underwent spatial sequencing, which provided images that allowed for the location of cell nuclei. "
    "Following this, an examination for the presence of antibodies was carried out, leading to the "
    "determination of the phenotype of individual cells.\n\n"
    "Tertiary lymphoid structures (TLS) are ectopic lymphoid formations which are formed under long-lasting "
    "inflammatory conditions, including tumours. TLS are composed predominantly of B cells, T cells and dendritic cells, "
    "and display various levels of organisation, from locally concentrated aggregates of immune cells, "
    "through clearly defined B cell follicles to mature follicles containing germinal centres. Their presence has been "
    "strongly associated with improved survival and clinical outcome upon cancer immunotherapies for patients with "
    "solid tumours, indicating potential for TLS to be used as a prognostic and predictive factor\n\n"
    "***Marta Trüb, Alfred Zippelius***, *Tertiary Lymphoid Structures as a Predictive Biomarker of Response to Cancer Immunotherapies*"
    )

fig_hist = px.histogram(data_IF1, x='celltype',
                        color='celltype',
                        labels={'celltype': 'Cell type'},
                        title="Cell type composition in F1 panel")
st.plotly_chart(fig_hist)

cell_types = ['all'] + list(data_IF1.celltype.unique())
cell_types_opt = st.multiselect(
    'Choose cell types',
    cell_types,
    default=['CD15-Tumor', 'CD15+Tumor', 'Tcell', 'Bcell']  # set 'all' as the default selected option
)

# if 'all' is selected, use all cell types; otherwise, use the selected ones
if 'all' in cell_types_opt:
    data_IF1_filtered = data_IF1
else:
    data_IF1_filtered = data_IF1[data_IF1.celltype.isin(cell_types_opt)]

fig_IF1 = px.scatter(data_IF1_filtered, x='nucleus.x', y='nucleus.y',
                     color="celltype",
                     title=f'Patient {patient} IF1 panel',
                     hover_data=["tissue.type", "celltype", "nucleus.x", "nucleus.y"],
                     opacity=0.5,
                     )
fig_IF1.update_layout(autosize=False, width=800, height=600)
st.plotly_chart(fig_IF1)


if st.sidebar.button('Show F2 and F3 panels'):
    st.session_state.display_result = True

if st.session_state.display_result:
    button_b = st.sidebar.button('Hide', on_click=hide_callback)
    data_IF2 = helper.get_panel('IF2', patient)
    data_IF3 = helper.get_panel('IF3', patient)
    st.markdown('### IF2 panel')
    fig_IF2 = px.scatter(data_IF2, x='nucleus.x', y='nucleus.y',
                         color="celltype",
                         title=f'Patient {patient} IF2 panel',
                         hover_data=["tissue.type", "celltype", "nucleus.x", "nucleus.y"]
                         )
    fig_IF2.update_layout(autosize=False, width=800, height=600)
    st.plotly_chart(fig_IF2)

    st.markdown('### IF3 panel')
    fig_IF3 = px.scatter(data_IF3, x='nucleus.x', y='nucleus.y',
                         color="celltype",
                         title=f'Patient {patient} IF3 panel',
                         hover_data=["tissue.type", "celltype", "nucleus.x", "nucleus.y"]
                         )
    fig_IF3.update_layout(autosize=False, width=800, height=600)
    st.plotly_chart(fig_IF3)

# features = [
#       'cell.ID', 'nucleus.x', 'nucleus.y', 'CD15.score.normalized',
#       'CK.score.normalized', 'CD3.score.normalized', 'CD11c.score.normalized',
#       'CD20.score.normalized', 'CD163.score.normalized', 'tissue.type',
#       'phenotype'
#   ]


radius = st.sidebar.slider("Select radius", 0, 100, 30)
min_cell = st.sidebar.slider("Select minimal number of cell in cluster", 10, 30, 20)

plot_b, cc_b = helper.plot_cc(data_IF1, radius=radius, types=['Bcell'], min_cell=min_cell)
plot_t, cc_t = helper.plot_cc(data_IF1, radius=radius, types=['Tcell'], min_cell=min_cell)
plot_bt, cc_bt = helper.plot_cc(data_IF1, radius=radius, types=['Bcell', 'Tcell', 'BnT'], min_cell=min_cell)

plot_b.update_layout(autosize=False, width=800, height=600)
plot_t.update_layout(autosize=False, width=800, height=600)
plot_bt.update_layout(autosize=False, width=800, height=600)

st.plotly_chart(plot_b)
st.plotly_chart(plot_t)
st.plotly_chart(plot_bt)

data_no_tumor = data_IF1[~data_IF1.celltype.isin(['CD15-Tumor', 'CD15+Tumor'])]
data_tumor = data_IF1[data_IF1.celltype.isin(['CD15-Tumor', 'CD15+Tumor'])]
cell_types = list(data_no_tumor.celltype.unique())

tls_df = helper.find_surroundings(data_no_tumor, cc_bt, radius=radius)


# Create a figure
data_tumor['cluster'] = data_tumor.celltype
tls_df['cluster'] = tls_df.TLS

greys = ['#c5d8ef', '#e0d2d1']
with open('color_pallete.txt', 'r') as f:
    color_pallet = f.read().strip().lower().split('\n')
color_pallet = greys + color_pallet

df = pd.concat([data_tumor, tls_df.drop(columns=['TLS'])])
tls_plot = px.scatter(df, x='nucleus.x', y='nucleus.y',
                      color='cluster',
                      color_discrete_sequence=color_pallet,
                      title='Localization of TLS in Tumor Tissue',
                      )

st.plotly_chart(tls_plot)


comp_plot = helper.plot_composition(tls_df, types=cell_types)
comp_plot.update_layout(autosize=False, width=1000, height=600)
st.plotly_chart(comp_plot)
