import json
import os.path
import pandas as pd
from ete3 import Tree


# DIVERSE
# ======================================================================================================================
def create_new_dir(dir_path, verbose=True):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        if verbose:
            print(f'{dir_path} directory created')
    else:
        if verbose:
            print(f'{dir_path} directory already exists')


SP_CHAR = ['-', '.', ' ']


def get_abbr_name(name, dict_assoc, prefix=None):
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
        elif new_name in dict_assoc.values():
            if len(e) > 4:
                new_name += e[0].upper() + e[1:4].lower()
            else:
                new_name += e[0].upper() + e[1:].lower()
    j = 0
    while new_name in dict_assoc.values():
        j += 1
        new_name = new_name[:-1] + str(j)
    return new_name


def create_abbreviation_names_dict(community_path, host_path, output_path):
    name_assoc = dict()
    host_list = [x.split('.')[0] for x in os.listdir(host_path)]
    for host in host_list:
        new_name = get_abbr_name(host, name_assoc)
        name_assoc[host] = new_name

    for comm_host in os.listdir(community_path):
        comm_list = [x.split('.')[0] for x in os.listdir(os.path.join(community_path, comm_host))]
        for comm in comm_list:
            new_name = get_abbr_name(comm, name_assoc, name_assoc[comm_host])
            name_assoc[f'{comm_host}_{comm}'] = new_name

    output_file = os.path.join(output_path, 'name_assoc.json')
    with open(output_file, 'w') as f:
        json.dump(name_assoc, fp=f, indent=4)
    return name_assoc


create_abbreviation_names_dict('../../example2/inputs/community',
                               '../../example2/inputs/hosts',
                               '../../example2/outputs')

# METABOLIC ANALYSIS
# ======================================================================================================================
def merge_outputs(file_cluster, file_info):
    df_clust = pd.read_csv(file_cluster, delimiter='\t', index_col='Compound')
    df_info = pd.read_csv(file_info, delimiter='\t', index_col='Compound')
    merge_df = df_clust.join(df_info)
    os.remove(file_cluster)
    os.remove(file_info)
    merge_df.to_csv(file_info, sep='\t')


def create_heatmap_output(output_path, output_name, analysis_method):
    output_heatmap = os.path.join(output_path, 'heatmap')
    create_new_dir(output_heatmap)
    output_heatmap = os.path.join(output_heatmap, analysis_method)
    create_new_dir(output_heatmap)
    output_heatmap = os.path.join(output_heatmap, output_name)
    return output_heatmap


# COEVOLUTION
# ======================================================================================================================

def get_host_comm_from_name(name):
    decomposition = name.split('_')
    host = decomposition[0]
    comm = decomposition[1] + '_' + decomposition[2]
    return host, comm


def col_normalization(df):
    for col in df.columns:
        mean_col = df[col].mean()
        df[col] = df[col]/mean_col
    return df
