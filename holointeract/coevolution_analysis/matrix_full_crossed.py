import os
import json
import matplotlib.pyplot as plt
from holointeract.utils.utils import *

COEVOLUTION_STR = 'coevolution'


def col_normalization(df):
    for col in df.columns:
        mean_col = df[col].mean()
        df[col] = df[col]/mean_col
    return df


def complementarity_boxplot(df, output_name):
    fig, ax = plt.subplots()
    ax.boxplot(df.values)
    ax.set_xticklabels(df.columns, rotation=90)
    plt.subplots_adjust(top=0.95, bottom=0.40)
    plt.title('Comparison of each host complementarity')
    plt.ylabel('Normalized complementarity')
    plt.xlabel('Host')
    plt.savefig(output_name+".png")


def generate_added_value_df(host_path: str, comm_path: str, scopes_path, output, name):
    """
    """
    host_list = [x.split('.')[0] for x in os.listdir(host_path)]
    comm_list = []
    for comm_host in os.listdir(comm_path):
        comm_list += [x.split('.')[0] for x in os.listdir(os.path.join(comm_path, comm_host))]

    df = pd.DataFrame(index=comm_list, columns=host_list)
    for host in host_list:
        for comm in comm_list:
            added_value_file = os.path.join(scopes_path, f'{host}__{comm}', 'community_analysis', 'addedvalue.json')
            with open(added_value_file, 'r') as f:
                data = json.load(f)
                df.loc[comm, host] = str(len(data['addedvalue']))
    df = df[df.columns].astype(float)
    df = col_normalization(df)

    output = os.path.join(output, COEVOLUTION_STR)
    create_new_dir(output)
    complementarity_boxplot(df, os.path.join(output, f'{name}_boxplot.png'))
    df.to_csv(os.path.join(output, f'{name}_matrix.tsv'), sep='\t')

