import os
import json
import pandas
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import seaborn
import scipy.cluster.hierarchy as sch

from padmet.utils.sbmlPlugin import convert_from_coded_id
from typing import List, Tuple, Set, Dict


def extract_metabolites(input_dir: str) \
        -> Tuple[Dict[str, Set[str]], Dict[str, Set[str]], Dict[str, Set[str]], List[str]]:
    """ Extract the metabolites linked to the hosts, bacteria and holobionts from m2m metacom results

    Parameters
    ----------
    input_dir : str
        Path to the directory where are stored all the m2m metacom results for each host

    Returns
    -------
    dict[str, set[str]]
        bact_metabolites : Associate for each host the set of all metabolites producible by all the bacteria
    dict[str, set[str]]
        host_metabolites : Associate for each host the set of all metabolites producible by the host
    dict[str, set[str]]
        holo_metabolites : Associate for each host the set of all metabolites producible by his holobiont
    list[str]
        all_metabolites : List of all the metabolites found in all networks (host + bacteria)
    """
    bact_metabolites = dict()
    host_metabolites = dict()
    holo_metabolites = dict()
    all_metabolites = set()
    for host in os.listdir(input_dir):
        bact_scope_file = os.path.join(input_dir, host, 'indiv_scopes', 'indiv_scopes.json')
        host_scope_file = os.path.join(input_dir, host, 'community_analysis', 'comm_scopes.json')
        holo_scope_file = os.path.join(input_dir, host, 'community_analysis', 'addedvalue.json')
        with open(bact_scope_file, 'r') as bact_f, \
             open(host_scope_file, 'r') as host_f, \
             open(holo_scope_file, 'r') as holo_f:
            bact_scope = rename_metabolites(set().union(*[set(x) for x in dict(json.load(bact_f)).values()]))
            host_scope = rename_metabolites(set(dict(json.load(host_f))['host_scope']))
            holo_scope = rename_metabolites(set(dict(json.load(holo_f))['addedvalue']))

            all_metabolites = all_metabolites.union(bact_scope, host_scope, holo_scope)
            bact_metabolites[host] = bact_scope
            host_metabolites[host] = host_scope
            holo_metabolites[host] = holo_scope

    return bact_metabolites, host_metabolites, holo_metabolites, list(all_metabolites)


def fill_matrix(bact_metabolites: Dict[str, Set[str]], host_metabolites: Dict[str, Set[str]],
                holo_metabolites: Dict[str, Set[str]], all_metabolites: List[str]) -> pandas.DataFrame:
    """ Fill a matrix indicating how each metabolite is produced for each host following the encoding :
    0 = Not produced
    1 = Only Holobiont
    2 = Only Host
    3 = Only Bact
    4 = Host AND Bact

    Parameters
    ----------
    bact_metabolites : dict[str, set[str]]
        Associate for each host the set of all metabolites producible by all the bacteria
    host_metabolites : dict[str, set[str]]
        Associate for each host the set of all metabolites producible by the host
    holo_metabolites : dict[str, set[str]]
        Associate for each host the set of all metabolites producible by his holobiont
    all_metabolites : list[str]
        List of all the metabolites found in all networks (host + bacteria)

    Returns
    -------
    pandas.DataFrame
        Matrix indicating how each metabolite (row) is produced for each host (column)
    """
    host_list = list(host_metabolites.keys())
    df = pandas.DataFrame(index=all_metabolites, columns=host_list, dtype=int)
    for host in host_list:
        column = list()
        for metabolite in all_metabolites:
            if metabolite in host_metabolites[host]:
                if metabolite in bact_metabolites[host]:
                    column.append(4)
                else:
                    column.append(2)
            elif metabolite in holo_metabolites[host]:
                column.append(1)
            elif metabolite in bact_metabolites[host]:
                column.append(3)
            else:
                column.append(0)
        df[host] = column
    return df


def rename_metabolites(metabolites: Set[str]) -> Set[str]:
    """ Rename the metabolites from SBML format

    Parameters
    ----------
    metabolites : Set[str]
        Set of metabolites to rename

    Returns
    -------
    Set[str]
        Set of renamed metabolites
    """
    renamed_metabolites = set()
    for met in metabolites:
        renamed_metabolites.add(convert_from_coded_id(met)[0].replace('_C-BOUNDARY', ''))
    return renamed_metabolites


def heatmap(df: pandas.DataFrame, output_heatmap: str, output_clusters: str,  method: str, max_clust: int):
    """Generate a heatmap from the DataFrame indicating how each metabolite is produced for each host
    0 = Not produced (grey : #adb8be)
    1 = Only Holobiont (green : #4c5923)
    2 = Only Host (blue : #162e4e)
    3 = Only Bact (red : #4b0e0e)
    4 = Host AND Bact (purple : #673c84)

    Parameters
    ----------
    df : pandas.DataFrame.
        DataFrame indicating how each metabolite is produced for each host
    method : str
        Clustering method
    output_heatmap : str
        Output cluster-map png file path
    output_clusters : str
        Output clusters tsv file path
    max_clust: int
    """
    dist = sch.distance.pdist(df)
    linkage = sch.linkage(dist, method=method)
    clusters = sch.fcluster(linkage, max_clust, 'maxclust')
    clusters_d = {df.index[i]: clusters[i] for i in range(len(clusters))}

    assoc_clust_color = dict()
    colors = list()
    n = max(clusters)
    for c in clusters:
        rgb_color = plt.cm.twilight_shifted(c/n)
        colors.append(rgb_color)
        assoc_clust_color[c] = rgb_color
    clusters_color = pandas.Series(colors, index=df.index)
    #        grey       green      blue       red        purple
    cmap = ['#adb8be', '#4c5923', '#162e4e', '#4b0e0e', '#673c84']
    ticks_labels = ['0 Not produced', '1 Holobiont only', '2 Host only', '3 Bacteria only', '4 Host AND Bacteria']
    df = df.T
    plot = seaborn.clustermap(df, cmap=cmap, cbar_pos=(0.01, 0.885, 0.01, 0.1), figsize=(30, 15), method=method,
                              tree_kws=dict(linewidths=1), col_colors=clusters_color, dendrogram_ratio=0.1,
                              vmin=-0.5, vmax=4.5, linewidths=0)
    plt.setp(plot.ax_cbar.set_yticks([0, 1, 2, 3, 4]))
    plt.setp(plot.ax_cbar.set_yticklabels(ticks_labels))
    plt.setp(plot.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)
    labels = sorted(assoc_clust_color.keys())
    handles = [Patch(facecolor=assoc_clust_color[c]) for c in labels]
    plt.legend(handles, labels, title='Cluster', bbox_to_anchor=(1, 1), bbox_transform=plt.gcf().transFigure,
               loc='upper right')

    plot.savefig(output_heatmap)

    with open(output_clusters, 'w') as f:
        f.write('Cluster\tCompound')
        met = df.columns
        for ind in plot.dendrogram_col.reordered_ind:
            f.write(f'\n{clusters_d[met[ind]]}\t{met[ind]}')


def heatmap_host_bacteria(input_dir: str, output: str, method: str = 'ward', max_clust: int = 10):
    """Generate a matrix in tsv format from m2m metacom outputs indicating how each metabolite is produced for each host

    Parameters
    ----------
    input_dir : str
        Path to the directory where are stored all the m2m metacom results for each host
    output : str
        Path and name to the output tsv and heatmap file to generate
    method : str
        Method for the heatmap
    max_clust: int
    """
    bact_metabolites, host_metabolites, holo_metabolites, all_metabolites = extract_metabolites(input_dir)
    df = fill_matrix(bact_metabolites, host_metabolites, holo_metabolites, all_metabolites)
    df.to_csv(f'{output}_matrix.tsv', sep='\t')
    heatmap(df, f'{output}_heatmap.png', f'{output}_clusters.tsv', method, max_clust)
    return bact_metabolites, host_metabolites, holo_metabolites, all_metabolites
