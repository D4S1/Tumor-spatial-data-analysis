import pandas as pd
from scipy.sparse.csgraph import connected_components
from sklearn.neighbors import radius_neighbors_graph, KDTree
import plotly.express as px

FILE_PATH = "data/IF_data/"


def standardize_phenotype(phenotype):
    markers = set()
    start = 0
    for i in range(len(phenotype)):
        if phenotype[i] == '-' or phenotype[i] == '+':
            markers.add(phenotype[start:i + 1])
            start = i + 1
    return ''.join(sorted(markers))


def get_panel(panel, patient):
    mapping = pd.read_csv(f'data/{panel}_phen_to_cell_mapping.csv', sep=',', header=0)
    mapping['phenotype'] = mapping['phenotype'].apply(lambda x: standardize_phenotype(x))

    data = pd.read_csv(f'{FILE_PATH}{panel}/{patient}_{panel}.csv', sep=',', header=0, index_col=0)
    data['phenotype'] = data['phenotype'].apply(lambda x: standardize_phenotype(x))
    data = data.merge(mapping, on='phenotype', how='left')
    return data[['nucleus.x', 'nucleus.y', 'celltype', "tissue.type"]]


def plot_cc(data, radius, types, min_cell):
    with open('color_pallete.txt', 'r') as f:
        color_pallet = f.read().strip().lower().split('\n')

    data = data[data.celltype.isin(types)]

    graph = radius_neighbors_graph(data[['nucleus.x', 'nucleus.y']], radius=radius, include_self=True)
    n_cc, cc_labels = connected_components(graph)

    components = [[] for _ in range(n_cc)]
    clusters = []

    for node, label in enumerate(cc_labels):
        clusters.append(label)
        components[label].append(node)

    tls_i = 1
    valid_clusters = []
    for i, cc in enumerate(components):
        name = 'unclassified'
        if len(cc) >= min_cell:
            name = f'TLS no. {tls_i}'
            tls_i += 1
        valid_clusters.append(name)
    cluster_labels = [valid_clusters[c] for c in clusters]
    data['TLS'] = cluster_labels

    grey_idx = cluster_labels.index('unclassified')
    grey_idx = grey_idx if len(set(cluster_labels[:grey_idx])) == grey_idx else len(set(cluster_labels[:grey_idx]))
    color_pallet.insert(grey_idx, '#c5d8ef')

    plot = px.scatter(data,
                      x='nucleus.x', y='nucleus.y',
                      color=cluster_labels,
                      color_discrete_sequence=color_pallet,
                      hover_data=['celltype'],
                      title=f"{types} cells connected components"
                      )
    data.sort_values(by='TLS', inplace=True)
    data = data[data.TLS != 'unclassified']
    return plot, data


def plot_composition(data: pd.DataFrame, types: list):
    # data needs to have columns: cell type and tls

    res = {cl: [0]*len(types) for cl in data.TLS.unique()}

    for row, tls in enumerate(data.TLS):
        cell_type = data.iloc[row].celltype
        cell_type = types.index(cell_type)
        res[tls][cell_type] += 1

    composition_df = pd.DataFrame(columns=['cell_type', 'cluster', 'percent'])
    for tls, counts in res.items():
        n_cell = sum(counts)
        for i, cell in enumerate(types):
            composition_df.loc[len(composition_df)] = [cell, tls, round(counts[i]*100/n_cell, 2)]

    comp_plot = px.bar(composition_df, x='cluster', y='percent',
                       color='cell_type',
                       title='Cluster composition',
                       )

    return comp_plot


def find_surroundings(data, tls_df, radius):
    full_tls = pd.DataFrame(columns=['nucleus.x', 'nucleus.y', 'celltype', 'TLS'])

    kd_tree = KDTree(data[['nucleus.x', 'nucleus.y']])

    for tls in tls_df.TLS.unique():
        points = tls_df[tls_df.TLS == tls][['nucleus.x', 'nucleus.y']]

        surr_tmp = kd_tree.query_radius(points, radius)
        surroundings = set()
        for point in surr_tmp:
            surroundings.update(list(point))
        surroundings = list(surroundings)

        df = data.iloc[surroundings][['nucleus.x', 'nucleus.y', 'celltype']]
        df['TLS'] = [tls]*len(df)
        full_tls = pd.concat([full_tls, df])

    return full_tls

# zrobić wykres: komórki rakowe szare tls:kolorowe
# zrobić bar plot (composition plot)
