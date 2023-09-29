import os
import plotly
import numpy as np
import plotly.graph_objs as go
import plotly.io as pio
import csv
import matplotlib.pyplot as plt
from holointeract.utils.utils import *
from statsmodels.stats.multitest import multipletests
from scipy.stats import linregress, spearmanr
from scipy.stats import f_oneway

# CONSTANT STR
# ======================================================================================================================
COEVOLUTION_STR = 'coevolution'
BENJAMINI = 'benjamini'
BONFERRONI = 'bonferroni'
SLOPE = 'slope'
P_VALUE = 'p_value'
SPEARMAN_C = 'spearman_coefficient'
X = 'x'
Y = 'y'
Y_PRED = 'y_pred'
LABEL = 'label'


# FUNCTIONS
# ======================================================================================================================
def coevolution(scopes_path, output, name, name_assoc, phylo_tree=None, correction=None):
    """
    """
    output = os.path.join(output, COEVOLUTION_STR)
    create_new_dir(output)

    complementarity_df = get_complementarity_df(scopes_path)
    complementarity_boxplot(complementarity_df, os.path.join(output, f'{name}_complementarity_boxplot.png'))
    complementarity_df.to_csv(os.path.join(output, f'{name}_complementarity_matrix.tsv'), sep='\t')
    if phylo_tree is not None:
        phylo_dist_df = get_phylo_dist_df(phylo_tree, name_assoc)
        linear_regression_complementarity_phylo_dist(complementarity_df, phylo_dist_df,
                                                     os.path.join(output, f'{name}_coevolution_regression'), correction)
        phylo_dist_df.to_csv(os.path.join(output, f'{name}_phylogenetic_dist_matrix.tsv'), sep='\t')


# MATRIX
# ======================================================================================================================
def get_complementarity_df(scopes_path):
    df = pd.DataFrame()
    for scope in os.listdir(scopes_path):
        host, comm = get_host_comm_from_name(scope)
        added_value_file = os.path.join(scopes_path, scope, 'community_analysis', 'addedvalue.json')
        with open(added_value_file, 'r') as f:
            data = json.load(f)
            df.loc[comm, host] = str(len(data['addedvalue']))
    df = df[df.columns].astype(float)
    df = col_normalization(df)
    return df


def get_phylo_dist_df(phylo_tree, name_assoc):
    df = pd.DataFrame()
    tree = Tree(phylo_tree)
    for host1 in tree.iter_leaves():
        host1_name = str(host1)[3:]
        for host2 in tree.iter_leaves():
            host2_name = str(host2)[3:]
            if host1_name in name_assoc and host2_name in name_assoc:
                if host1_name != host2_name:
                    df.loc[name_assoc[host1_name], name_assoc[host2_name]] = host1.get_distance(host2)
                else:
                    df.loc[name_assoc[host1_name], name_assoc[host2_name]] = float(0)
    return df


# FIGURES
# ======================================================================================================================

def complementarity_boxplot(df, output):
    fig, ax = plt.subplots()
    ax.boxplot(df.values)
    ax.set_xticklabels(df.columns, rotation=90)
    plt.subplots_adjust(top=0.95, bottom=0.40)
    plt.title('Comparison of each host complementarity')
    plt.ylabel('Normalized complementarity')
    plt.xlabel('Host')
    plt.savefig(output)


def linear_regression_complementarity_phylo_dist(complementarity_df, phylo_dist_df, output, correction=None):
    # ANOVA test
    print(f_oneway(*[complementarity_df[host].values for host in complementarity_df.columns]))

    n_micro, n_host = complementarity_df.shape
    colors = plotly.colors.qualitative.Dark2
    significant_slopes = 0

    with open(f'{output}.tsv', 'w', newline='') as tsv_file:
        writer = csv.writer(tsv_file, delimiter='\t')
        writer.writerow(['Microorganism', 'Slope', 'Spearman Correlation', 'P-value', 'Corrected P-value'])

        results, p_values = linear_regression_results(complementarity_df, phylo_dist_df)

        fig = go.Figure()
        fig.update_layout(xaxis_title="Phylogenetic distance (between host and natural microorganism host)",
                          yaxis_title="Metabolic complementarity (normalized)",
                          title="Linear regression for each microorganism")

        if correction == BENJAMINI:
            multiple_tests_cor = multipletests(p_values, alpha=0.05, method="fdr_bh", is_sorted=False,
                                               returnsorted=False)
            for i in range(len(complementarity_df.index)):
                writer.writerow([microorganism, results[microorganism][SLOPE], results[microorganism][SPEARMAN_C],
                                 results[microorganism][P_VALUE], multiple_tests_cor[1][i]])
                if multiple_tests_cor[0][i]:
                    microorganism = complementarity_df.index[i]
                    significant_slopes += 1
                    color = colors[significant_slopes % len(colors)]
                    add_slope(fig, results, microorganism, color)

        elif correction == BONFERRONI:
            for microorganism in complementarity_df.index:
                writer.writerow([microorganism, results[microorganism][SLOPE], results[microorganism][SPEARMAN_C],
                                 results[microorganism][P_VALUE], results[microorganism][P_VALUE] * n_host])
                if results[microorganism][P_VALUE] < (0.05 / n_host) and results[microorganism][SLOPE] < 0:
                    significant_slopes += 1
                    color = colors[significant_slopes % len(colors)]
                    add_slope(fig, results, microorganism, color)

        elif correction is None:
            for microorganism in complementarity_df.index:
                writer.writerow([microorganism, results[microorganism][SLOPE], results[microorganism][SPEARMAN_C],
                                 results[microorganism][P_VALUE], 'No correction'])
                if results[microorganism][P_VALUE] < 0.05 and results[microorganism][SLOPE] < 0:
                    significant_slopes += 1
                    color = colors[significant_slopes % len(colors)]
                    add_slope(fig, results, microorganism, color)

    pio.write_html(fig, f'{output}.html')


def linear_regression_results(complementarity_df, phylo_dist_df):
    results = {}
    p_values = []
    for microorganism in complementarity_df.index:
        microorganism_host = microorganism.split('_')[0]
        # Extract complementarity line values from Data Frame
        comm_line = complementarity_df.loc[microorganism]
        # Get phylogenetic distance values between each host and the comm host
        x_phylo_dist = [phylo_dist_df.loc[host, microorganism_host] for host in comm_line.index]
        y_complementarity = comm_line.values
        label = ['Host = ' + i for i in comm_line.index]
        # Linear regression
        slope, intercept, r_value, p_value, std_err = linregress(x_phylo_dist, y_complementarity)
        y_prediction = slope * np.array(x_phylo_dist) + intercept
        # Calculate Spearman coefficient and p-value
        spearman_coefficient, spearman_p_value = spearmanr(x_phylo_dist, y_complementarity)
        p_values.append(spearman_p_value)
        # Store result for the
        results[microorganism] = {P_VALUE: spearman_p_value,
                                  X: x_phylo_dist,
                                  Y: y_complementarity,
                                  Y_PRED: y_prediction,
                                  SLOPE: slope,
                                  LABEL: label,
                                  SPEARMAN_C: spearman_coefficient}
    return results, p_values


def add_slope(fig, results, microorganism, color):
    # Add regression line to plot
    fig.add_trace(go.Scatter(x=results[microorganism][X], y=results[microorganism][Y_PRED], mode='lines',
                             name=microorganism, line=dict(color=color)))
    # Add points to plot
    fig.add_trace(go.Scatter(x=results[microorganism][X], y=results[microorganism][Y],
                             text=results[microorganism][LABEL],
                             mode='markers', name=microorganism, marker=dict(color=color)))


SC_PATH = '../../example2/outputs/scopes/full'
PHY_FILE = '../../example/inputs/SpeciesTree_rooted.txt'
with open('../../example2/outputs/name_assoc.json', 'r') as f:
    NASS = json.load(f)
coevolution(SC_PATH, '', 'test', NASS, PHY_FILE)
