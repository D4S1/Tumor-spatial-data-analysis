import pandas as pd
import os
import plotly.express as px
import streamlit as st

from helper import get_panel

file_path = "../data/IF_data/"

patients = []
for filename in os.listdir(file_path + "IF1"):
    patients.append(filename[:4])


st.title('Tumor spatial data analysis')
patient = st.sidebar.selectbox('Select a patient', patients)

data_IF1 = get_panel('IF1', patient)
data_IF2 = get_panel('IF2', patient)
data_IF3 = get_panel('IF3', patient)


st.markdown('### IF1 panel')

cell_types = ['all'] + list(data_IF1.celltype.unique())
cell_types_opt = st.multiselect(
    'Choose cell types',
    cell_types,
    default='all'  # set 'all' as the default selected option
)

# if 'all' is selected, use all cell types; otherwise, use the selected ones
if 'all' in cell_types_opt:
    data_IF1_filtered = data_IF1
else:
    data_IF1_filtered = data_IF1[data_IF1.celltype.isin(cell_types_opt)]

fig_IF1 = px.scatter(data_IF1_filtered, x='nucleus.x', y='nucleus.y', color="celltype", title=f'Patient {patient} IF1 panel')
fig_IF1.update_layout(autosize=False, width=800, height=600)
st.plotly_chart(fig_IF1)

st.markdown('### IF2 panel')

fig_IF2 = px.scatter(data_IF2, x='nucleus.x', y='nucleus.y', color="celltype", title=f'Patient {patient} IF2 panel')
fig_IF2.update_layout(autosize=False, width=800, height=600)
st.plotly_chart(fig_IF2)

st.markdown('### IF3 panel')

fig_IF3 = px.scatter(data_IF3, x='nucleus.x', y='nucleus.y', color="celltype", title=f'Patient {patient} IF3 panel')
fig_IF3.update_layout(autosize=False, width=800, height=600)
st.plotly_chart(fig_IF3)
# features = [
#       'cell.ID', 'nucleus.x', 'nucleus.y', 'CD15.score.normalized',
#       'CK.score.normalized', 'CD3.score.normalized', 'CD11c.score.normalized',
#       'CD20.score.normalized', 'CD163.score.normalized', 'tissue.type',
#       'phenotype'
#   ]




