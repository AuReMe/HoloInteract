import os

import pandas as pd
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
from typing import Dict


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


# MAIN FUNCTION
# ======================================================================================================================
def coevolution(scopes_path: str, output: str, name: str, name_assoc: Dict[str, str], phylo_tree: str = None,
                correction: str = None):
    """ Performs coevolution analysis : build matrix and generate boxplot and linear regression figures. Looking for
    coevolution signs by searching correlation between complementarity for a microorganism and a host AND phylogenetic
    distance between the natural microorganism host and a host.

    Parameters
    ----------
    scopes_path: str
        Path to the directory with calculated scopes (with full method)
    output: str
        Path to the directory to store the output results
    name: str
        Name of the analysis for output files name
    name_assoc: Dict[str, str]
        Dictionary associating for each host and each microorganism name, its abbreviated name
    phylo_tree: str, optional (default=None)
        Path to phylogenetic tree of hosts at Newick format (if None, linear regression analysis will not be performed)
    correction: str, optional (default=None)
        Type of correction to apply to p-values for the multiple linear regression performed (bonferroni, benjamini or
        None for no correction)
    """
    output = os.path.join(output, COEVOLUTION_STR)
    create_new_dir(output)
    # Performs complementarity analysis (matrix + boxplot)
    complementarity_df = get_complementarity_df(scopes_path)
    complementarity_boxplot(complementarity_df, os.path.join(output, f'{name}_complementarity_boxplot.png'))
    complementarity_df.to_csv(os.path.join(output, f'{name}_complementarity_matrix.tsv'), sep='\t')
    # Performs linear regression if phylogenetic tree parameter given (+phylogenetic distance matrix creation)
    if phylo_tree is not None:
        phylo_dist_df = get_phylo_dist_df(phylo_tree, name_assoc)
        linear_regression_complementarity_phylo_dist(complementarity_df, phylo_dist_df,
                                                     os.path.join(output, f'{name}_coevolution_regression'), correction)
        phylo_dist_df.to_csv(os.path.join(output, f'{name}_phylogenetic_dist_matrix.tsv'), sep='\t')


# MATRIX FUNCTIONS
# ======================================================================================================================
def get_complementarity_df(scopes_path: str) -> pd.DataFrame:
    """ Generate the complementarity matrix. For each couple host / microorganism calculate its complementarity
    normalizing the number of compounds produces by added value (cooperation) in the scope.

    Parameters
    ----------
    scopes_path: str
        Path to the directory with calculated scopes (with full method)

    Returns
    -------
    pd.DataFrame
        Complementarity matrix
    """
    df = pd.DataFrame()
    for scope in os.listdir(scopes_path):
        host, comm = get_host_comm_from_name(scope)
        added_value_file = os.path.join(scopes_path, scope, 'community_analysis', 'addedvalue.json')
        with open(added_value_file, 'r') as f:
            data = json.load(f)
            # Store the number of compounds produces by added value (cooperation)
            df.loc[comm, host] = str(len(data['addedvalue']))
    df = df[df.columns].astype(float)
    # Normalize values in the matrix
    df = col_normalization(df)
    return df


def get_phylo_dist_df(phylo_tree: str, name_assoc: Dict[str, str]) -> pd.DataFrame:
    """ Creates a matrix calculating phylogenetic distance between each pair of host.

    Parameters
    ----------
    phylo_tree: str
        Path to phylogenetic tree of hosts at Newick format
    name_assoc: str
        Dictionary associating for each host and each microorganism name, its abbreviated name

    Returns
    -------
    pd.DataFrame
        matrix associating phylogenetic distance between each pair of host
    """
    df = pd.DataFrame()
    # Load tree with ete3
    tree = Tree(phylo_tree)
    for host1 in tree.iter_leaves():
        for host2 in tree.iter_leaves():
            # Checks if both host are in name assoc dict (<==> in the run analyse)
            if host1.name in name_assoc and host2.name in name_assoc:
                # Store distance between hosts in matrix
                if host1.name != host2.name:
                    df.loc[name_assoc[host1.name], name_assoc[host2.name]] = host1.get_distance(host2)
                else:
                    df.loc[name_assoc[host1.name], name_assoc[host2.name]] = float(0)
    return df


# FIGURES FUNCTIONS
# ======================================================================================================================

# BOXPLOT
def complementarity_boxplot(df, output):
    fig, ax = plt.subplots()
    ax.boxplot(df.values)
    ax.set_xticklabels(df.columns, rotation=90)
    plt.subplots_adjust(top=0.95, bottom=0.40)
    plt.title('Comparison of each host complementarity')
    plt.ylabel('Normalized complementarity')
    plt.xlabel('Host')
    plt.savefig(output)


# LINEAR REGRESSION
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
