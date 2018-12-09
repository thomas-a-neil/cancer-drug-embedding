import os

from bokeh.plotting import figure, show, output_file
from bokeh.models import HoverTool, ColumnDataSource, CategoricalColorMapper, LinearColorMapper
from bokeh.palettes import Spectral10, Dark2
import pandas as pd

molecule_images_dir = 'molecule_images'


def get_molecule_images_png_path(nsc_id):
    return os.path.join(molecule_images_dir, '{}.png'.format(nsc_id))


def write_bokeh_plot(embedding_data_dict, umap_embedding, nsc_ids, is_div5, num_atoms_df, k=None, output_filepath='hover_plot_is_div5.html'):
    if not k:
        k = len(embedding_data_dict)
    output_file(output_filepath)

    molecule_df = pd.DataFrame(umap_embedding[:k], columns=('x', 'y'))
    molecule_df['is_div5'] = [str(x) for x in is_div5[:k]]
    molecule_df['nsc_id'] = [str(x) for x in nsc_ids[:k]]
    molecule_df['image'] = list(map(get_molecule_images_png_path, nsc_ids[:k]))

    datasource = ColumnDataSource(molecule_df)
    color_mapping = CategoricalColorMapper(factors=['0', '1'], palette=Dark2[3])

    plot_figure = figure(
        title='UMAP projection of the graph2vec embedded NCI open compound dataset. Colored by Open vs Diversity Set',
        plot_width=1000,
        plot_height=1000,
        tools=('pan, wheel_zoom, reset')
    )

    plot_figure.add_tools(HoverTool(tooltips="""
    <div>
        <div>
            <img src='@image' style='float: left; margin: 5px 5px 5px 5px'/>
        </div>
        <div>
            <span style='font-size: 16px; color: #224499'>NSC ID:</span>
            <span style='font-size: 18px'>@nsc_id</span>
        </div>
    </div>
    """))

    plot_figure.circle(
        'x',
        'y',
        source=datasource,
        color=dict(field='is_div5', transform=color_mapping),
        line_alpha=0.6,
        fill_alpha=0.6,
        size=4
    )
    show(plot_figure)


def write_continuous_bokeh_plot(embedding_data_dict, umap_embedding, nsc_ids, is_div5, num_atoms_df, k=None, output_filepath='hover_plot_nsc_id_continuous_coloring.html'):
    if k is None:
        k = len(embedding_data_dict)

    output_file(output_filepath)

    molecule_df = pd.DataFrame(umap_embedding[:k], columns=('x', 'y'))
    molecule_df['is_div5'] = [str(x) for x in is_div5[:k]]
    molecule_df['nsc_id'] = [x for x in nsc_ids[:k]]
    molecule_df['nsc_id_str'] = [str(x) for x in nsc_ids[:k]]
    molecule_df['image'] = list(map(get_molecule_images_png_path, nsc_ids[:k]))

    datasource = ColumnDataSource(molecule_df)
    color_mapping = LinearColorMapper(low=min(nsc_ids[:k]), high=max(nsc_ids[:k]), palette=Spectral10)

    plot_figure = figure(
        title='UMAP projection of the graph2vec embedded NCI open compound dataset. Colored by NSC ID',
        plot_width=1000,
        plot_height=1000,
        tools=('pan, wheel_zoom, reset')
    )

    plot_figure.add_tools(HoverTool(tooltips="""
    <div>
        <div>
            <img src='@image' style='float: left; margin: 5px 5px 5px 5px'/>
        </div>
        <div>
            <span style='font-size: 16px; color: #224499'>NSC ID:</span>
            <span style='font-size: 18px'>@nsc_id_str</span>
        </div>
    </div>
    """))

    plot_figure.circle(
        'x',
        'y',
        source=datasource,
        color=dict(field='nsc_id', transform=color_mapping),
        line_alpha=0.6,
        fill_alpha=0.6,
        size=4
    )
    show(plot_figure)
