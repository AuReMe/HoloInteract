from scipy.optimize import curve_fit
import pandas
import matplotlib.pyplot as plt
import numpy as np
from python_scripts.dist import get_dist
from scipy.stats import linregress
import scipy.stats as stats
import plotly.graph_objs as go
import plotly.io as pio
from statsmodels.stats.multitest import multipletests
import csv
from scipy.stats import linregress, spearmanr
from scipy.stats import f_oneway
import sys


def func(x, a, b, c):
    return a * np.exp(-b * np.array(x)) + c


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

    ####################

    with open(output+".csv", "w") as fo:
        fo.write("Couple,Complementarite,Distance\n")

    dico_dist = get_dist(phylogenetic_tree)

    for index, row in data.iterrows():
        print(index, row)
        for column, value in row.iteritems():
            if "-".join(index.split("-")[0].split("_")) == "-".join(column.split("_")):
                distance = 0
            else:
                distance = dico_dist["-".join(index.split("-")
                                              [0].split("_"))]["-".join(column.split("_"))]
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

    spacing = donnees.shape[1] - 1
    print("spacing", spacing)
    n_groups = donnees.shape[0]

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

        for i in range(n_groups-1):
            print(i)
            start = i*spacing
            end = (i+1)*spacing + 1

            group_data = donnees.iloc[start:end]
            print("AAAAAA")
            x_data = group_data["Distance"].to_numpy()
            print("xdata", x_data)
            y_data = group_data["Complementarite"].to_numpy()
            print("ydata", y_data)
            # Fit de la droite de régression
            if x_data != [] and y_data != []:
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
                    if spearman_p < (0.05) and slope < 0:
                        # Plot de la droite de régression et des points correspondants
                        fig.add_trace(go.Scatter(x=x_data, y=y_pred, name=couple_name, line=dict(
                            color=colors[counter % len(colors)])))

                        if show_points:
                            fig.add_trace(go.Scatter(x=x_data, y=y_data, mode='markers', name=couple_name +
                                                     " points", marker=dict(color=colors[counter % len(colors)])))

    pio.write_html(fig, output+'.html', auto_open=True)
    print(spacing)


def stats_tests(input_file):
    donnees = pandas.read_csv(input_file, sep=",", header=0)
    score_c = donnees["Complementarite"].to_list()

    distance = donnees["Distance"].to_list()

    # Test de corrélation de Spearman
    rho, p_value = stats.spearmanr(distance, score_c)
    print("Coefficient de corrélation de Spearman : ", rho)
    print("Valeur-p : ", p_value)


def job(phylogenetic_tree, input_file, ouput_file_for_graph, correction, graph_name="graph_complementarite_distance"):
    get_graph_csv(input_file=input_file, output=ouput_file_for_graph,
                  phylogenetic_tree=phylogenetic_tree)
    plot_regression_all(ouput_file_for_graph+".csv", output=graph_name,
                        base_file=input_file, show_points=True, correction=correction)


if __name__ == "__main__":

    job(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])

    # #space_to_remove=get_graph_csv("mat_dist_nonorm.csv", output="outlier_comp_dist_nonorm")
    # #print("nonorm_space_remove", space_to_remove)
    # #space_to_remove=get_graph_csv("mat_dist_norm_fullmicro.csv", output="outlier_comp_dist_norm_fullmicro")
    # #print("normfullmicro_space_remove", space_to_remove)
    # print("normmean_space_remove", space_to_remove)

    # #get_graph_csv("Short_dist_full_crossed.csv", output="Short_Comp_Dist")
    # #get_graph_csv("Long_dist_full_crossed.csv", output="Long_Comp_Dist")

    # space_to_remove=1+space_to_remove

    # #plot_regression_all("outlier_comp_dist_nonorm.csv",output="outlier_NoNorm_graph",base_file="mat_dist_nonorm.csv", show_points=True, correction="benjamini", spacing_remove=space_to_remove)
    # plot_regression_all("outlier_comp_dist_norm_mean.csv",output="outlier_Norm_Mean_graph",base_file="mat_dist_norm_mean.csv", show_points=True, correction="", spacing_remove=space_to_remove)
    # #plot_regression_all("outlier_comp_dist_norm_fullmicro.csv",output="outlier_Norm_Fullmicro_graph",base_file="mat_dist_norm_fullmicro.csv", show_points=True,correction="benjamini", spacing_remove=space_to_remove)

    # #plot_regression_all("Long_Comp_Dist.csv",output="Long_graph",base_file="Long_dist_full_crossed.csv", show_points=True)
    # #plot_regression_all("Short_Comp_Dist.csv",output = "Short_graph",base_file="Short_dist_full_crossed.csv", show_points=True)

    # #stats_tests("test_comp_dist_newnorm.csv")
    # #stats_tests("Long_Comp_Dist.csv")
    # #stats_tests("outlier_comp_dist_norm_meansans_outlier.csv")

    # #show_graph("outlier_comp_dist_norm_meansans_outlier.csv")
