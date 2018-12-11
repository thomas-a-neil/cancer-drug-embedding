import json

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import umap

sns.set(style='white', rc={'figure.figsize': (14, 10)})
dims_epochs_prefix = 'dims_128_epochs_10'
prefix = 'deg_4_' + dims_epochs_prefix
graph2vec_embedding_path = 'data/embeddings/degree_4/nci_open_training_data_{}_lr_0.3_embeddings.txt'.format(dims_epochs_prefix)
nsc_ids_path = 'data/div5_nsc_ids.txt'
num_atoms_path = 'data/num_atoms_to_nsc.csv'



def get_div5_binary_labels(embedding_data_dict, div5_nsc_ids_df):
    div5_nsc_ids_set = set(div5_nsc_ids_df['nsc_id'].values)
    is_div5 = []
    for filepath in embedding_data_dict.keys():
        if get_nsc_id_from_path(filepath) in div5_nsc_ids_set:
            is_div5.append(1)
        else:
            is_div5.append(0)
    return np.array(is_div5)


def get_nsc_id_from_path(path):
    filename = path.split('/')[-1]
    no_extensions = filename.split('.')[0]
    return int(no_extensions)


# TODO increase DPI of saved images (280k points needs slightly higher quality)
def plot_with_colorbar(embedding, values_for_coloring, path, label=''):
    # cmap will generate a tuple of RGBA values for a given number in the range 0.0 to 1.0
    # (also 0 to 255 - not used in this example).
    # To map our z values cleanly to this range, we create a Normalize object.
    cmap = matplotlib.cm.get_cmap('viridis')
    normalize = matplotlib.colors.Normalize(vmin=min(values_for_coloring), vmax=max(values_for_coloring))
    colors = [cmap(normalize(value)) for value in values_for_coloring]

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.scatter(embedding[:, 0], embedding[:, 1], color=colors, s=0.02)

    # Optionally add a colorbar
    cax, _ = matplotlib.colorbar.make_axes(ax)
    cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize, label=label)
    plt.savefig(path)
    plt.close()


def plot_div5(embedding, is_div5, path):
    plt.scatter(embedding[:, 0], embedding[:, 1],
                c=[sns.color_palette()[x] for x in is_div5],
                s=[2.0*x + 0.02*(1 - x) for x in is_div5])
    plt.savefig(path)
    plt.close()


def filter_by_num_atoms(embedding_data_dict, num_atoms_df):
    to_embed_vectors = []
    for vector, num_atoms in zip(embedding_data_dict.values(), num_atoms_df['num_atoms'].values):
        if num_atoms < 50:
            to_embed_vectors.append(vector)
    graph_embedding_data = np.vstack((np.array(x) for x in to_embed_vectors))
    return graph_embedding_data


# TODO: save embedding (json, pickle) so that I can redo plots with the same embedding coordinates
# TODO: don't recompute embedding in this filtered step
def plot_filtered_umap_embedding(embedding_data_dict, num_atoms_df, path):
    filtered_num_atoms = [x for x in num_atoms_df['num_atoms'].values if x < 50]
    graph_embedding_data = filter_by_num_atoms(embedding_data_dict, num_atoms_df)
    umap_embedding = umap.UMAP().fit_transform(graph_embedding_data)
    plot_with_colorbar(umap_embedding, filtered_num_atoms, path, label='Number of Atoms')


def load_data():
    div5_nsc_ids_df = pd.read_csv(nsc_ids_path, header=None, names=['nsc_id'])

    with open(graph2vec_embedding_path, 'r') as f:
        embedding_data_dict = json.load(f)
    nsc_ids = [get_nsc_id_from_path(filepath) for filepath in embedding_data_dict.keys()]

    graph_embedding_data = np.vstack((np.array(x) for x in embedding_data_dict.values()))
    is_div5 = get_div5_binary_labels(embedding_data_dict, div5_nsc_ids_df)

    umap_embedding = umap.UMAP().fit_transform(graph_embedding_data)

    num_atoms_df = pd.read_csv(num_atoms_path, header=None, names=['num_atoms', 'nsc_id'])
    return embedding_data_dict, umap_embedding, nsc_ids, is_div5, num_atoms_df


def remake_plots():
    embedding_data_dict, umap_embedding, nsc_ids, is_div5, num_atoms_df = load_data()
    plot_div5(umap_embedding, is_div5, path='plots/{}_molecular_umap_with_div5.png'.format(prefix))
    # plot_with_colorbar(umap_embedding, num_atoms_df['num_atoms'].values, 'plots/{}_num_atoms_colorbar.png'.format(prefix), label='Number of Atoms')
    plot_with_colorbar(umap_embedding, nsc_ids, 'plots/{}_molecular_umap_nsc_colorbar.png'.format(prefix), label='NSC ID')
    plot_filtered_umap_embedding(embedding_data_dict, num_atoms_df, path='plots/{}_num_atoms_le_50.png'.format(prefix))
    # num_atoms stuff
    # plt.hist([x for x in num_atoms_df['num_atoms'].values if x < 50])
    # plt.savefig('plots/num_atoms_hist_ge_50.png')
    return embedding_data_dict, umap_embedding, nsc_ids, is_div5, num_atoms_df

# embedding_data_dict, umap_embedding, nsc_ids, is_div5, num_atoms_df = remake_plots()
# k = 50000
# write_bokeh_plot(embedding_data_dict, umap_embedding, nsc_ids, is_div5, num_atoms_df, k=k, output_filepath='{}_{}_hover_plot_is_div5.html'.format(k, prefix))
