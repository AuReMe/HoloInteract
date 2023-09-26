from scipy.optimize import curve_fit
import pandas
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import plotly.graph_objs as go
import plotly.io as pio
from statsmodels.stats.multitest import multipletests
import csv
from scipy.stats import linregress, spearmanr
from scipy.stats import f_oneway
import sys
import os
from ete3 import Tree


def func(x, a, b, c):
    return a * np.exp(-b * np.array(x)) + c


def get_dist_from_tree(tree: str):
    tree_file = os.path.join(tree)
    dict_tree = dict()
    tree = Tree(tree_file)
    for sp1 in tree.iter_leaves():
        dict_dist = dict()
        for sp2 in tree.iter_leaves():
            if sp1 != sp2:
                dict_dist[str(sp2)[3:]] = sp1.get_distance(sp2)
        dict_tree[str(sp1)[3:]] = dict_dist
    return dict_tree


def show_graph(input_file):

    fig, ax = plt.subplots(figsize=(10, 10))
    donnees = pandas.read_csv(input_file, sep=",", header=0)

    score_c = donnees["Complementarite"].to_list()
    distance = donnees["Distance"].to_list()

    slope, intercept = np.polyfit(distance, score_c, 1)
    print("pente: ", slope)
    print("intercept: ", intercept)

    popt, pcov = curve_fit(func, xdata=distance, ydata=score_c)

    # calcul de la distance des points à la courbe de régression
    """diff = np.abs(score_c - func(distance, *popt))"""

    # colormap pour colorer les points selon leur distance à la courbe de régression
    """cmap = plt.colormaps.get_cmap('viridis')
    norm = plt.Normalize(vmin=diff.min(), vmax=diff.max())"""

    #ax.scatter(distance, score_c,s=200, c=diff, cmap=cmap, norm=norm)
    ax.scatter(distance, score_c, s=2, c="black")
    ax.plot(distance, slope*np.array(distance) +
            intercept, color="red", linewidth=2)
    ax.plot(distance, func(distance, *popt),
            'r-', label='Courbe de régression')

    ax.set_xlabel("Distance phylogénétique")
    ax.set_ylabel("Complementarite métabolique normalisée")

    ax.set_title("Observation de la complémentarité métabolique de chaque bactérie avec chaque algue \nen fonction de la distance phylogénétique des algues")

    plt.savefig("noerror_courbe_nodup.png")


def get_graph_csv(input_file: str, output: str, phylogenetic_tree: str):
    data = pandas.read_csv(input_file, sep=",", header=0, index_col=0)
    print(data.columns)

    with open(output+".csv", "w") as fo:
        fo.write("Couple,Complementarite,Distance\n")

    dico_dist = get_dist_from_tree(phylogenetic_tree)

    for index, row in data.iterrows():
        print(index, row)
        for column, value in row.iteritems():
            if index.split("-")[0] == column:
                distance = 0
            else:
                distance = dico_dist[index.split(
                    "-")[0]][column]
            with open(output+".csv", "a") as fo:
                fo.write(f'{index}_{column},{float(value)},{distance}\n')


def plot_regression_all(input_file, output, base_file, show_points=False, correction="", spacing_remove=0):

    donnees = pandas.read_csv(input_file, sep=",", header=0)

    groups_data = pandas.read_csv(base_file, sep=",", header=0)

    pvalues = []

    # Vérification statistique des données
    # Effectuer l'ANOVA pour chaque colonne
    results = {}
    for col in groups_data.columns:
        print(col)

        result = f_oneway(
            *[groups_data[col].values for col in groups_data.columns[1:]])
        results[col] = {'Statistique F': result.statistic,
                        'Valeur p': result.pvalue}

    # Afficher les résultats
    for col, result in results.items():
        print(result)
        print("Colonne :", col)
        print("Statistique F :", result['Statistique F'])
        print("Valeur p :", result['Valeur p'])
        print()

    #####

    spacing = groups_data.shape[1]-1
    print("spacing", spacing)
    n_groups = donnees.shape[0]//spacing

    fig = go.Figure()
    fig.update_layout(
        xaxis_title="Distance phylogénétique",
        yaxis_title="Complémentarité métabolique normalisée",
        title="Régression linéaire pour chaque bactérie"
    )

    # Définir une liste de couleurs
    colors = ['red', 'blue', 'green', 'purple', 'orange', "brown"]
    counter = 0
    with open(output+'Regression_correlation.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            ['Couple', 'Pente', 'Corrélation de Spearman', 'P-value', 'P-value corrigée'])

    # Correction avant calculs

        if correction == "benjamini":
            for i in range(n_groups):
                start = i*spacing
                end = (i+1)*spacing

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
            start = i*spacing
            end = (i+1)*spacing

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
                bonferroni = len(groups_data)
                # Écriture des informations dans le fichier CSV
                couple_name = " ".join(
                    donnees["Couple"].iloc[i*spacing].split("_")[0:-3])
                writer.writerow(
                    [couple_name, slope, spearman_coef, spearman_p, spearman_p*bonferroni])
                if spearman_p < (0.05/bonferroni) and slope < 0:
                    counter += 1
                # Plot de la droite de régression et des points correspondants
                    fig.add_trace(go.Scatter(x=x_data, y=y_pred, name=couple_name, line=dict(
                        color=colors[counter % len(colors)])))

                    if show_points:
                        fig.add_trace(go.Scatter(x=x_data, y=y_data, mode='markers', name=couple_name +
                                                 " points", marker=dict(color=colors[counter % len(colors)])))
                """elif spearman_p < (0.05/bonferroni) and slope > 0:
                    print(couple_name)"""

            elif correction == "benjamini":
                couple_name = " ".join(
                    donnees["Couple"].iloc[i*spacing].split("_")[0:-3])
                writer.writerow(
                    [couple_name, slope, spearman_coef, spearman_p, corrected_pvalues[i]])
                if rejected[i] == True:
                    counter += 1
                # Plot de la droite de régression et des points correspondants
                    fig.add_trace(go.Scatter(x=x_data, y=y_pred, name=couple_name, line=dict(
                        color=colors[counter % len(colors)])))

                    if show_points:
                        fig.add_trace(go.Scatter(x=x_data, y=y_data, mode='markers', name=couple_name +
                                                 " points", marker=dict(color=colors[counter % len(colors)])))
                """elif corrected_pvalues[i] < (0.05) and slope > 0:
                    print(couple_name)"""

            elif correction == "":
                couple_name = " ".join(
                    donnees["Couple"].iloc[i*spacing].split("_")[0:-3])
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
                                                 " points", marker=dict(color=colors[counter % len(colors)])))

    pio.write_html(fig, output+'file.html', auto_open=True)
    print(spacing)


def stats_tests(input_file):
    donnees = pandas.read_csv(input_file, sep=",", header=0)
    score_c = donnees["Complementarite"].to_list()

    distance = donnees["Distance"].to_list()

    # Test de corrélation de Spearman
    rho, p_value = stats.spearmanr(distance, score_c)
    print("Coefficient de corrélation de Spearman : ", rho)
    print("Valeur-p : ", p_value)


def job(phylogenetic_tree: str, input_file: str, ouput_file_for_graph: str, correction: str, graph_name="graph_complementarite_distance"):
    """Create the html graph to analyse coevolution

    Args:
        phylogenetic_tree (str): phylogenetic tree of hosts (Newick format)
        input_file (str): matrix built by 'mat_dist_full_crossed.job' function
        ouput_file_for_graph (str): name of csv file for graph
        correction (str): multiple test correction 
        graph_name (str, optional): name for html file. Defaults to "graph_complementarite_distance".
    """
    get_graph_csv(input_file=input_file, output=ouput_file_for_graph,
                  phylogenetic_tree=phylogenetic_tree)
    plot_regression_all(ouput_file_for_graph+".csv", output=graph_name,
                        base_file=input_file, show_points=True, correction=correction)


if __name__ == "__main__":

    job(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
