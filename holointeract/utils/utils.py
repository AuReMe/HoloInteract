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
