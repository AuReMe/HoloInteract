import os.path
import pandas as pd


# DIVERSE
# ======================================================================================================================
def create_new_dir(dir_path):
    if not os.path.exists(dir_path):
        os.makedirs(dir_path)
        # print(f'{dir_path} directory created')
    else:
        print(f'{dir_path} directory already exists')


def get_abbr_name(name, dict_assoc):
    new_name = ''
    decomposition = name.split('_')
    size = len(decomposition)
    for i in range(size):
        if i != size - 1:
            new_name += decomposition[i][0].upper()
        else:
            new_name += decomposition[i][:4].lower()
            j = 0
            while new_name in dict_assoc.keys():
                j += 1
                new_name = new_name[:-1] + str(j)
    return new_name


def create_abbreviation_names_dict(community_path, host_path):
    host_assoc = dict()
    host_list = [x.split('.')[0] for x in os.listdir(host_path)]
    for host in host_list:
        new_name = get_abbr_name(host, host_assoc)
        host_assoc[host] = new_name

    comm_assoc = dict()
    for comm_host in os.listdir(community_path):
        comm_list = [x.split('.')[0] for x in os.listdir(os.path.join(community_path, comm_host))]
        for comm in comm_list:
            new_name = get_abbr_name(comm, comm_assoc)
            new_name = f'{host_assoc[comm_host]}_{new_name}'
            comm_assoc[f'{comm_host}_{comm}'] = new_name
    return host_assoc, comm_assoc


COM_PATH = '../../example/inputs/community'
HOS_PATH = '../../example/inputs/hosts'
create_abbreviation_names_dict(COM_PATH, HOS_PATH)


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
