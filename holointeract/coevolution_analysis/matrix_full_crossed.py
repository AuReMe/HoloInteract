import os
import json
import matplotlib.pyplot as plt
from holointeract.utils.utils import *

COEVOLUTION_STR = 'coevolution'


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


def complementarity_boxplot(df, output_name):
    fig, ax = plt.subplots()
    ax.boxplot(df.values)
    ax.set_xticklabels(df.columns, rotation=90)
    plt.subplots_adjust(top=0.95, bottom=0.40)
    plt.title('Comparison of each host complementarity')
    plt.ylabel('Normalized complementarity')
    plt.xlabel('Host')
    plt.savefig(output_name+".png")


def generate_added_value_df(scopes_path, output, name):
    """
    """
    host_set = set()
    comm_set = set()
    value_dict = {}
    for scope in os.listdir(scopes_path):
        host, comm = get_host_comm_from_name(scope)
        host_set.add(host)
        comm_set.add(comm)

        added_value_file = os.path.join(scopes_path, scope, 'community_analysis', 'addedvalue.json')
        with open(added_value_file, 'r') as f:
            data = json.load(f)
            value_dict[scope] = str(len(data['addedvalue']))

    df = pd.DataFrame(index=list(comm_set), columns=list(host_set))
    for couple, value in value_dict.items():
        host, comm = get_host_comm_from_name(couple)
        df.loc[comm, host] = value
    df = df[df.columns].astype(float)
    df = col_normalization(df)

    output = os.path.join(output, COEVOLUTION_STR)
    create_new_dir(output)
    complementarity_boxplot(df, os.path.join(output, f'{name}_boxplot.png'))
    df.to_csv(os.path.join(output, f'{name}_matrix.tsv'), sep='\t')
    return df

