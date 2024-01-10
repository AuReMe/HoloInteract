import json
import logging
import os.path
import pandas as pd
from ete3 import Tree
from typing import Dict, Tuple


SP_CHAR = ['-', '.', ' ']


# DIVERSE
# ==================================================================================================
def are_inputs_ok(community_sbml_path: str, hosts_sbml_path: str) -> bool:
    """ Checks if inputs file names coincides between community sbml path and hosts sbml path.

    Parameters
    ----------
    community_sbml_path: str
        Path to the directory containing SBML networks for all community of microorganisms
        classified in their natural host directory
    hosts_sbml_path: str
        Path to the directory containing SBML networks for all hosts

    Returns
    -------
    bool
        True if host species names coincides in the 2 input directories.

    Raises
    ------
    ValueError
        If at least 1 host is not in community directory
    ValueError
        If at least 1 host's community is not in hosts directory
    """
    hosts_sbml = set([x.split('.')[0] for x in os.listdir(hosts_sbml_path)])
    community_sbml = set(os.listdir(community_sbml_path))
    if hosts_sbml == community_sbml:
        return True
    elif hosts_sbml.difference(community_sbml) != set():
        raise ValueError(f'{hosts_sbml.difference(community_sbml)} not in community input path.')
    elif community_sbml.difference(hosts_sbml) != set():
        raise ValueError(f'{community_sbml.difference(hosts_sbml)} not in hosts input path.')


def create_new_dir(dir_path: str, verbose: bool = True):
    """ Create a directory checking if already existing.

    Parameters
    ----------
    dir_path: str
        Path of the directory to create
    verbose: str
        To print the directory creation
    """
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        if verbose:
            logging.info(f'{dir_path} directory created')
    else:
        if verbose:
            logging.info(f'{dir_path} directory already exists')


def get_abbr_name(name: str, name_assoc: Dict[str, str], prefix: str = None) -> str:
    """ Returns the abbreviated name of a species (host or microorganism) name

    Parameters
    ----------
    name: str
        Original name of the species
    name_assoc: Dict[str, str]
        Dictionary associating for each host and each microorganism name, its abbreviated name
    prefix: str
        Prefix to use for the abbreviated name (= natural host of a microorganism)

    Returns
    -------
    str
        Abbreviated name of the species
    """
    new_name = ''
    if prefix is not None:
        new_name += prefix + '_'
    for c in SP_CHAR:
        name = name.replace(c, '_')
    decomposition = name.split('_')
    for e in decomposition:
        if e == decomposition[0]:
            new_name += e[0].upper()
        elif e == decomposition[1]:
            if len(e) > 4:
                new_name += e[0].upper() + e[1:4].lower()
            else:
                new_name += e[0].upper() + e[1:].lower()
        elif new_name in name_assoc.values():
            if len(e) > 4:
                new_name += e[0].upper() + e[1:4].lower()
            else:
                new_name += e[0].upper() + e[1:].lower()
    j = 0
    while new_name in name_assoc.values():
        j += 1
        new_name = new_name[:-1] + str(j)
    return new_name


def create_abbreviation_names_dict(community_sbml_path: str, hosts_sbml_path: str,
                                   output_path: str) -> Dict[str, str]:
    """ Associating in a dictionary to each host and each community microorganism, its abbreviated
    name used for the analysis. Store the dictionary in a json file.

    Parameters
    ----------
    community_sbml_path: str
        Path to the directory containing SBML networks for all community of microorganisms
        classified in their natural host directory
    hosts_sbml_path: str
        Path to the directory containing SBML networks for all hosts
    output_path: str
        Path to store the output json file containing the dictionary of names associations.

    Returns
    -------
    Dict[str, str]
        Dictionary associating for each host and each microorganism name, its abbreviated name
    """
    are_inputs_ok(community_sbml_path, hosts_sbml_path)
    name_assoc = dict()
    host_list = [x.split('.')[0] for x in os.listdir(hosts_sbml_path)]
    for host in host_list:
        new_name = get_abbr_name(host, name_assoc)
        name_assoc[host] = new_name

    for comm_host in os.listdir(community_sbml_path):
        comm_list = [x.split('.')[0] for x in os.listdir(os.path.join(community_sbml_path, comm_host))]
        for comm in comm_list:
            new_name = get_abbr_name(comm, name_assoc, name_assoc[comm_host])
            name_assoc[f'{comm_host}_{comm}'] = new_name

    output_file = os.path.join(output_path, 'name_assoc.json')
    with open(output_file, 'w') as f:
        json.dump(name_assoc, fp=f, indent=4)
    return name_assoc


def load_name_assoc_file(name_assoc_directory_path: str) -> Dict[str, str]:
    """ Load the name_assoc.json file containing the dictionary associating for each host and each
    microorganism name, its abbreviated name.

    Parameters
    ----------
    name_assoc_directory_path: str
        Directory where is stored the name_assoc.json file.

    Returns
    -------
    Dict[str, str]
        Dictionary associating for each host and each microorganism name, its abbreviated name
    """
    name_assoc_file_path = os.path.join(name_assoc_directory_path, 'name_assoc.json')
    with open(name_assoc_file_path, 'r') as f:
        name_assoc_dict = json.load(f)
    return name_assoc_dict


# METABOLIC ANALYSIS
# ==================================================================================================
def merge_outputs(file_cluster: str, file_info: str):
    """ Merge the 2 outputs files :
        - <output_heatmap>_clusters.tsv from heatmap_host_bacteria function
        - <output_info>.tsv from proportion_workflow function
    according to the scope compounds column.
    Writes the merged matrix to a new file <name>_classes_cpd_info.tsv

    Parameters
    ----------
    file_cluster: str
        Path to <output_heatmap>_clusters.tsv from heatmap_host_bacteria function
    file_info: str
        Path to <output_info>.tsv from proportion_workflow function
    """
    df_clust = pd.read_csv(file_cluster, delimiter='\t', index_col='Compound')
    df_info = pd.read_csv(file_info, delimiter='\t', index_col='Compound')
    merge_df = df_clust.join(df_info)
    os.remove(file_cluster)
    os.remove(file_info)
    merge_df.to_csv(file_info, sep='\t')


def create_heatmap_output(output_path: str, output_name: str, analysis_method: str):
    """ Creates outputs directories for the heatmap analysis.

    Parameters
    ----------
    output_path: str
        Main output directory path
    output_name: str
        Output name of the analysis
    analysis_method: str
        Method used to generate scopes used as input for heatmap creation

    Returns
    -------
    str
        Directory to store heatmap analysis output files
    """
    output_heatmap = os.path.join(output_path, 'heatmap')
    create_new_dir(output_heatmap)
    output_heatmap = os.path.join(output_heatmap, analysis_method)
    create_new_dir(output_heatmap)
    output_heatmap = os.path.join(output_heatmap, output_name)
    return output_heatmap


# COEVOLUTION
# ==================================================================================================

def get_host_microorganism_from_name(name: str) -> Tuple[str, str]:
    """ Get the host AND microorganism names from the aggregate files names.

    Parameters
    ----------
    name: str
        Aggregate file name

    Returns
    -------
    str
        Name of the host
    str
        Name of the microorganism : natural host + '_' + microorganism
    """
    decomposition = name.split('_')
    host = decomposition[0]
    comm = decomposition[1] + '_' + decomposition[2]
    return host, comm


def col_normalization(complementarity_df: pd.DataFrame) -> pd.DataFrame:
    """ Apply normalization to columns of complementarity matrix.

    Parameters
    ----------
    complementarity_df: pd.DataFrame
        Complementarity matrix indicating for each couple host (col) / microorganism (row) its
        complementarity (not normalized)

    Returns
    -------
    pd.DataFrame
        Complementarity matrix indicating for each couple host (col) / microorganism (row) its
        complementarity (normalized)
    """
    for col in complementarity_df.columns:
        mean_col = complementarity_df[col].mean()
        complementarity_df[col] = complementarity_df[col]/mean_col
    return complementarity_df


def phylo_tree_names_ok(name_assoc: Dict[str, str], phylo_tree: str) -> bool:
    """ Verify if the phylogenetic tree file contains all the host species of the study.

    Parameters
    ----------
    name_assoc: Dict[str, str]
        Dictionary associating for each host and each microorganism name, its abbreviated name
    phylo_tree: str
        Path to phylogenetic tree of hosts at Newick format

    Returns
    -------
    bool
        True if all hosts of the study are in the phylogenetic tree leaves.

    Raises
    ------
    ValueError
        If one of the host of the study is not in the phylogenetic tree leaves.
    """
    tree = Tree(phylo_tree)
    tree_species = set(tree.get_leaf_names())
    for species, abbr_name in name_assoc.items():
        if '_' not in abbr_name and species not in tree_species:
            raise ValueError(f'Host {species} not in {phylo_tree} tree leaves. Careful of character'
                             f' cast and special characters.')
    return True

