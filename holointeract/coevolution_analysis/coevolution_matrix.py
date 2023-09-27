import os
from holointeract.utils.utils import *


import pandas
import numpy as np
import plotly.graph_objs as go
import plotly.io as pio
from statsmodels.stats.multitest import multipletests
import csv
from scipy.stats import linregress, spearmanr
from scipy.stats import f_oneway
import matplotlib.pyplot as plt

COEVOLUTION_STR = 'coevolution'


def coevolution(scopes_path, output, name, name_assoc, phylo_tree=None):
    """
    """
    output = os.path.join(output, COEVOLUTION_STR)
    create_new_dir(output)

    complementarity_df = get_complementarity_df(scopes_path)
    complementarity_boxplot(complementarity_df, os.path.join(output, f'{name}_complementarity_boxplot.png'))
    complementarity_df.to_csv(os.path.join(output, f'{name}_complementarity_matrix.tsv'), sep='\t')
    if phylo_tree is not None:
        phylo_dist_df = get_phylo_dist_df(phylo_tree, name_assoc)
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

def func(x, a, b, c):
    return a * np.exp(-b * np.array(x)) + c


def complementarity_boxplot(df, output):
    fig, ax = plt.subplots()
    ax.boxplot(df.values)
    ax.set_xticklabels(df.columns, rotation=90)
    plt.subplots_adjust(top=0.95, bottom=0.40)
    plt.title('Comparison of each host complementarity')
    plt.ylabel('Normalized complementarity')
    plt.xlabel('Host')
    plt.savefig(output)


def correlation_complementarity_phylo_dist(complementarity_df, phylo_dist_df):
    pass

# OLD
# ======================================================================================================================
def plot_regression_all(input_file, output, base_file, show_points=False, correction="", spacing_remove=0):
    donnees = pandas.read_csv(input_file, sep=",", header=0)

    added_value_df = pandas.read_csv(base_file, sep=",", header=0)

    pvalues = []

    # Vérification statistique des données
    # Effectuer l'ANOVA pour chaque colonne
    results = {}
    for host in added_value_df.columns:
        print(host)

        result = f_oneway(
            *[added_value_df[col].values for col in added_value_df.columns[1:]])
        results[host] = {'Statistique F': result.statistic,
                         'Valeur p': result.pvalue}

    # Afficher les résultats
    for host, result in results.items():
        print(result)
        print("Colonne :", host)
        print("Statistique F :", result['Statistique F'])
        print("Valeur p :", result['Valeur p'])
        print()

    #####

    spacing = added_value_df.shape[1] - 1
    print("spacing", spacing)
    n_groups = donnees.shape[0] // spacing

    fig = go.Figure()
    fig.update_layout(
        xaxis_title="Distance phylogénétique",
        yaxis_title="Complémentarité métabolique normalisée",
        title="Régression linéaire pour chaque bactérie"
    )

    # Définir une liste de couleurs
    colors = ['red', 'blue', 'green', 'purple', 'orange', "brown"]
    counter = 0
    with open(output + 'Regression_correlation.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            ['Couple', 'Pente', 'Corrélation de Spearman', 'P-value', 'P-value corrigée'])

        # Correction avant calculs

        if correction == "benjamini":
            for i in range(n_groups):
                start = i * spacing
                end = (i + 1) * spacing

                group_data = donnees.iloc[start:end]
                x_data = group_data["Distance"].to_numpy()
                y_data = group_data["Complementarite"].to_numpy()
                # Fit de la droite de régression
                slope, intercept, r_value, p_value, std_err = linregress(
                    x_data, y_data)
                y_pred = slope * x_data + intercept

                # Coefficient de corrélation de Spearman et valeur p
                pvalues.append(spearmanr(x_data, y_data)[1])
            resultat = multipletests(
                pvalues, alpha=0.05, method="fdr_bh", is_sorted=False, returnsorted=False)
            rejected = resultat[0]

            corrected_pvalues = resultat[1]

        for i in range(n_groups):
            print(i)
            start = i * spacing
            end = (i + 1) * spacing

            group_data = donnees.iloc[start:end]

            print("AAAAAA", group_data)
            x_data = group_data["Distance"].to_numpy()
            print("xdata", x_data)
            y_data = group_data["Complementarite"].to_numpy()
            print("ydata", y_data)
            # Fit de la droite de régression
            slope, intercept, r_value, p_value, std_err = linregress(
                x_data, y_data)
            y_pred = slope * x_data + intercept

            # Coefficient de corrélation de Spearman et valeur p
            spearman_coef, spearman_p = spearmanr(x_data, y_data)

            # Verification par un test de Bonferoni
            if correction == "bonferroni":
                bonferroni = len(added_value_df)
                # Écriture des informations dans le fichier CSV
                couple_name = " ".join(
                    donnees["Couple"].iloc[i * spacing].split("_")[0:-3])
                writer.writerow(
                    [couple_name, slope, spearman_coef, spearman_p, spearman_p * bonferroni])
                if spearman_p < (0.05 / bonferroni) and slope < 0:
                    counter += 1
                    # Plot de la droite de régression et des points correspondants
                    fig.add_trace(go.Scatter(x=x_data, y=y_pred, name=couple_name, line=dict(
                        color=colors[counter % len(colors)])))

                    if show_points:
                        fig.add_trace(go.Scatter(x=x_data, y=y_data, mode='markers', name=couple_name +
                                                                                          " points",
                                                 marker=dict(color=colors[counter % len(colors)])))
                """elif spearman_p < (0.05/bonferroni) and slope > 0:
                    print(couple_name)"""

            elif correction == "benjamini":
                couple_name = " ".join(
                    donnees["Couple"].iloc[i * spacing].split("_")[0:-3])
                writer.writerow(
                    [couple_name, slope, spearman_coef, spearman_p, corrected_pvalues[i]])
                if rejected[i] == True:
                    counter += 1
                    # Plot de la droite de régression et des points correspondants
                    fig.add_trace(go.Scatter(x=x_data, y=y_pred, name=couple_name, line=dict(
                        color=colors[counter % len(colors)])))

                    if show_points:
                        fig.add_trace(go.Scatter(x=x_data, y=y_data, mode='markers', name=couple_name +
                                                                                          " points",
                                                 marker=dict(color=colors[counter % len(colors)])))
                """elif corrected_pvalues[i] < (0.05) and slope > 0:
                    print(couple_name)"""

            elif correction == "":
                couple_name = " ".join(
                    donnees["Couple"].iloc[i * spacing].split("_")[0:-3])
                writer.writerow(
                    [couple_name, slope, spearman_coef, spearman_p])
                counter += 1
                print(spearman_p, type(spearman_p))
                if (spearman_p < (0.05) or np.isnan(spearman_p)) and slope < 0:
                    # Plot de la droite de régression et des points correspondants
                    fig.add_trace(go.Scatter(x=x_data, y=y_pred, name=couple_name, line=dict(
                        color=colors[counter % len(colors)])))

                    if show_points:
                        fig.add_trace(go.Scatter(x=x_data, y=y_data, mode='markers', name=couple_name +
                                                                                          " points",
                                                 marker=dict(color=colors[counter % len(colors)])))

    pio.write_html(fig, output + 'file.html', auto_open=True)
    print(spacing)


def job(input_file: str, ouput_file_for_graph: str, correction: str):
    """Create the html graph to analyse coevolution

    Args:
        input_file (str): matrix built by 'mat_dist_full_crossed.job' function
        ouput_file_for_graph (str): name of csv file for graph
        correction (str): multiple test correction
    """

    plot_regression_all(ouput_file_for_graph + ".csv",
                        base_file=input_file, show_points=True, correction=correction)
