import json
import pickle
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pandas as pd
import umap
from rdkit.Chem import PandasTools

sns.set(style='white', rc={'figure.figsize': (14, 10)})

all_nci_sdf_path = '/data/pass_data_800k.sdf'
all_nci_sdf_pandas_pickle_path = '/data/pass_data_800k.pkl'


feature_names = ['E_CAS',
                 'E_COMPLEXITY',
                 'E_DRUGLIKENESS',
                 'E_DRUGLIKENESS/2',
                 'E_LOGP',
                 'E_LOGP/2',
                 'E_LOGP/3',
                 'E_NHACCEPTORS',
                 'E_NHDONORS',
                 'E_NROTBONDS',
                 'E_NSC',
                 'E_PASS_DATA_PA',
                 'E_PASS_DATA_PI',
                 'E_SMILES',
                 'E_WEIGHT',
                 'ID',
                 'Molecule',
                 'SMILES']


def get_nsc_id_from_path(path):
    filename = path.split('/')[-1]
    no_extensions = filename.split('.')[0]
    return int(no_extensions)


# TODO increase DPI of saved images (280k points needs slightly higher quality)
def plot_with_colorbar(embedding, values_for_coloring, path, label='', vmin=None, vmax=None):
    # cmap will generate a tuple of RGBA values for a given number in the range 0.0 to 1.0
    # (also 0 to 255 - not used in this example).
    # To map our z values cleanly to this range, we create a Normalize object.
    cmap = matplotlib.cm.get_cmap('viridis')
    vmin = vmin if vmin is not None else min(values_for_coloring)
    vmax = vmax if vmax is not None else max(values_for_coloring)
    normalize = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    colors = [cmap(normalize(value)) for value in values_for_coloring]

    fig, ax = plt.subplots(figsize=(14, 10))
    ax.scatter(embedding['x'].values, embedding['y'].values, color=colors, s=0.02)

    # Optionally add a colorbar
    cax, _ = matplotlib.colorbar.make_axes(ax)
    cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize, label=label)
    plt.savefig(path)
    plt.close()


def get_dataframe_with_pass(umap_embedding, nsc_ids, pass_df):
    molecule_df = pd.DataFrame(umap_embedding, columns=('x', 'y'))
    molecule_df['NSC'] = nsc_ids
    pass_df['NSC'] = list(map(int, pass_df['E_NSC']))

    joined_df = molecule_df.merge(pass_df, on='NSC', how='left')
    return joined_df


# TODO: save embedding (json, pickle) so that I can redo plots with the same embedding coordinates
def load_data(graph2vec_embedding_path, umap_load_path=None, pass_df_load_path=None):
    with open(graph2vec_embedding_path, 'r') as f:
        embedding_data_dict = json.load(f)

    nsc_ids = [get_nsc_id_from_path(filepath) for filepath in embedding_data_dict.keys()]

    graph_embedding_data = np.vstack((np.array(x) for x in embedding_data_dict.values()))
    if umap_load_path:
        with open(umap_load_path, 'rb') as f:
            umap_embedding = pickle.load(f)
    else:
        umap_embedding = umap.UMAP().fit_transform(graph_embedding_data)

    if pass_df_load_path:
        with open(pass_df_load_path, 'rb') as f:
            pass_df = pickle.load(f)
    else:
        pass_df = PandasTools.LoadSDF(all_nci_sdf_path, smilesName='SMILES', molColName='Molecule')
        pass_df = split_out_pass_features(type_inference_pass_df(pass_df))

    joined_df = get_dataframe_with_pass(umap_embedding, nsc_ids, pass_df)
    return embedding_data_dict, umap_embedding, nsc_ids, pass_df, joined_df


all_feature_names = ['E_CAS',
                     'E_COMPLEXITY',
                     'E_DRUGLIKENESS',
                     'E_DRUGLIKENESS/2',
                     'E_LOGP',
                     'E_LOGP/2',
                     'E_LOGP/3',
                     'E_NHACCEPTORS',
                     'E_NHDONORS',
                     'E_NROTBONDS',
                     'E_NSC',
                     'E_PASS_DATA_PA',
                     'E_PASS_DATA_PI',
                     'E_SMILES',
                     'E_WEIGHT',
                     'ID',
                     'Molecule',
                     'SMILES']


def convert_logp_to_float(logp_string):
    if not type(logp_string) == str:
        if np.isnan(logp_string):
            return float('nan')
        else:
            print(logp_string)
            raise Exception
    return float(logp_string.split()[0])


plot_feature_name_converter = [float,
                               # 'i dont know druglikeness',
                               convert_logp_to_float,
                               int,
                               int,
                               int,
                               convert_logp_to_float]

# plot_feature_names = ['E_COMPLEXITY',
#                      # 'E_DRUGLIKENESS',
#                      'E_LOGP',
#                      'E_NHACCEPTORS',
#                      'E_NHDONORS',
#                      'E_NROTBONDS',
#                      'E_WEIGHT']

# I just want one.
plot_feature_names = ['E_COMPLEXITY',
                      'E_LOGP',
                      'E_WEIGHT']

# these bastards are long strings that we have to index into
pass_feature_names = ['E_PASS_DATA_PA',
                      'E_PASS_DATA_PI']

pass_feature_indexes = {'Antiinflammatory': 171,
                        'Carcinogenic': 494,
                        'Corticosteroid-like': 214,
                        'Immunomodulator': 48,
                        'Immunostimulant': 179,
                        'Immunosuppressant': 65,
                        'Mutagenic': 208}


embedding_paths = ['/data/graph2vec_embeddings/degree_4/nci_open_training_data_dims_128_epochs_10_lr_0.3_embeddings.txt',
                   '/data/graph2vec_embeddings/degree_4/nci_open_training_data_dims_128_epochs_100_lr_0.3_embeddings.txt',
                   '/data/graph2vec_embeddings/degree_3/nci_open_training_data_dims_256_epochs_10_lr_0.3_embeddings.txt',
                   '/data/graph2vec_embeddings/degree_3/nci_open_training_data_dims_256_epochs_100_lr_0.3_embeddings.txt',
                   ]

prefixes = ['deg_4_dims_128_epochs_10',
            'deg_4_dims_128_epochs_100',
            'deg_3_dims_256_epochs_10',
            'deg_3_dims_256_epochs_100']


def type_inference_pass_df(pass_df):
    for feature, converter in zip(plot_feature_names, plot_feature_name_converter):
        pass_df[feature] = list(map(converter, pass_df[feature]))
    return pass_df


def split_out_pass_features(pass_df):
    for feature_name in pass_feature_names:
        pa_or_pi = feature_name[-2:]
        prefix = 'Predicted Activity ' if pa_or_pi == 'PA' else 'Predicted Inactivity '
        for name, index in pass_feature_indexes.items():
            def to_map(pass_feature_string):
                if not type(pass_feature_string) == str:
                    if np.isnan(pass_feature_string):
                        return float('nan')
                    else:
                        print(pass_feature_string)
                        raise Exception
                return float(pass_feature_string.split()[index])

            pass_df[prefix + name] = list(map(to_map, pass_df[feature_name].values))
    return pass_df


def plot_feature(joined_df, feature_name, path, label):
    print('Saving {} to {}'.format(feature_name, path))
    values_for_coloring = joined_df[feature_name].values
    # for those long tailed distributions (like logp)
    vmin = np.percentile(values_for_coloring, 1)
    vmax = np.percentile(values_for_coloring, 99)
    plot_with_colorbar(joined_df, values_for_coloring, path, label=label, vmin=vmin, vmax=vmax)


def plot_hist(joined_df, feature_name, path):
    plt.hist(joined_df[feature_name].values)
    print('Saving {} histogram to {}'.format(feature_name, path))
    plt.savefig(path)
    plt.close()


def save_umap(umap_embedding, save_dir):
    save_path = os.path.join(save_dir, 'umap.pkl')
    with open(save_path, 'wb') as f:
        pickle.dump(umap_embedding, f)
    print('Saving UMAP to {}'.format(save_path))
    return save_path


def save_sdf_pandas(joined_df, save_path):
    with open(save_path, 'wb') as f:
        pickle.dump(joined_df, f)
    print('Saving SDF PASS Data Frame to {}'.format(save_path))
    return save_path


def get_feature_plots(joined_df, save_dir):
    for feature_name in plot_feature_names:
        dropped_df = joined_df.dropna(subset=[feature_name])
        plot_feature(dropped_df, feature_name, os.path.join(save_dir, feature_name + '.png'), label=feature_name)
        plot_hist(dropped_df, feature_name, os.path.join(save_dir, feature_name + '_histogram.png'))


def get_pass_feature_plots(joined_df, save_dir):
    for prefix in ['Predicted Activity ', 'Predicted Inactivity ']:
        for feature_name in pass_feature_indexes.keys():
            label = prefix + feature_name
            column_name = prefix + feature_name  # done in the split step
            # -1.0 is a null value
            no_neg1_df = joined_df[joined_df[column_name] != -1.0]
            dropped_df = no_neg1_df.dropna(subset=[column_name])
            filename = '_'.join(prefix.split() + [feature_name]) + '.png'
            plot_feature(dropped_df, column_name, os.path.join(save_dir, filename), label=label)


def plot_div5(embedding, is_div5, path):
    print('Saving div5 plot to {}'.format(path))
    plt.scatter(embedding[:, 0], embedding[:, 1],
                c=[sns.color_palette()[x] for x in is_div5],
                s=[2.0*x + 0.02*(1 - x) for x in is_div5])
    plt.savefig(path)
    plt.close()


def get_div5_binary_labels(embedding_data_dict, div5_nsc_ids_df):
    div5_nsc_ids_set = set(div5_nsc_ids_df['nsc_id'].values)
    is_div5 = []
    for filepath in embedding_data_dict.keys():
        if get_nsc_id_from_path(filepath) in div5_nsc_ids_set:
            is_div5.append(1)
        else:
            is_div5.append(0)
    return np.array(is_div5)


def get_all_feature_plots():
    for prefix, embedding_path in zip(prefixes, embedding_paths):
        print(prefix)
        print(embedding_path)
        save_dir = '/data/{}/'.format(prefix)
        try:
            os.mkdir(save_dir)
        except FileExistsError:
            # our directory is already there
            pass
        umap_load_path = os.path.join(save_dir, 'umap.pkl')
        pass_df_load_path = '/data/pass_data_800k.pkl'
        embedding_data_dict, umap_embedding, nsc_ids, pass_df, joined_df = load_data(embedding_path, umap_load_path=umap_load_path, pass_df_load_path=pass_df_load_path)
        # umap_load_path = save_umap(umap_embedding, save_dir)
        # print(umap_load_path)
        # pass_df_load_path = save_sdf_pandas(pass_df, all_nci_sdf_pandas_pickle_path)
        # get_feature_plots(joined_df, save_dir)
        # get_pass_feature_plots(joined_df, save_dir)

        # div5 stuff
        nsc_ids_path = '/data/div5_nsc_ids.txt'
        div5_nsc_ids_df = pd.read_csv(nsc_ids_path, header=None, names=['nsc_id'])
        is_div5 = get_div5_binary_labels(embedding_data_dict, div5_nsc_ids_df)
        plot_div5(umap_embedding, is_div5, os.path.join(save_dir, 'is_div5.png'))



