import seaborn as sns
import sys
import pandas


def heatmap(input_file, method, output_file, color="tab10"):
    """Genere une clustermap du fichier csv donné en entree
    Args:
        input_file (str): fichier csv
        method (str): methode de clusterisation
        color (str): palette de couleur Seaborn à utiliser. Defauts : tab10
        output_file (str): nom du fichier PNG de sortie.
    """

    matrice = pandas.read_csv(input_file, sep=";", header=0, index_col=0)
    matrice = matrice.T
    palette = ["#1f77b4", "#2ca02c", "#9467bd", "#e377c2", "#17becf"]

    larg = matrice.shape[1]*0.2
    haut = matrice.shape[0]*6
    if matrice.shape[0] > 20:
        haut = 120

    clustermap = sns.clustermap(
        matrice,
        cmap=palette,
        method=method,
        cbar_pos=(0.01, 0.9, 0.05, 0.1),
        figsize=(larg, haut),  # taille de la figure en pouces
        dendrogram_ratio=(0.1, 0.1),
        yticklabels=True,
        xticklabels=True,

        cbar_kws={'label': 'Valeurs'},
        tree_kws=dict(linewidths=10)
    )
    clustermap.fig.suptitle(
        "Heatmap de la méthode de production des molécules")

    clustermap.ax_heatmap.set_yticklabels(
        clustermap.ax_heatmap.get_yticklabels(), rotation=0, fontsize=20)
    clustermap.ax_row_dendrogram.set_visible(True)
    clustermap.ax_col_dendrogram.set_visible(True)

    clustermap.savefig(output_file + ".png")


def analyse_tab_pourc(input_file, output_name):
    # Charger les données à partir du fichier CSV
    df = pandas.read_csv(input_file, sep=";", index_col=0)
    df = df.T
    # Créer une clustermap avec Seaborn
    g = sns.clustermap(df, cmap="viridis", method="ward", figsize=(280, 100), dendrogram_ratio=(0.2, 0.3),
                       yticklabels=True,
                       xticklabels=True,

                       cbar_kws={'label': 'Valeurs'},
                       tree_kws=dict(linewidths=10)
                       )

    # Personnaliser l'apparence du graphique
    g.fig.suptitle("Heatmap de la méthode de production des molécules")
    g.ax_heatmap.set_xlabel("Molécules", fontsize=200)
    g.ax_heatmap.set_ylabel("Méthode de production", fontsize=200)
    g.ax_heatmap.tick_params(axis="x", labelsize=10)
    g.ax_heatmap.tick_params(axis="y", labelsize=200)

    # Convertir le graphique en HTML zoomable
    g.savefig(output_name+".png")
    # plt.show()


if __name__ == "__main__":
    if sys.argv[1] == "clustermap":
        heatmap(input_file=sys.argv[2], method=sys.argv[3],
                output_file=sys.argv[4], color=sys.argv[5])
    elif sys.argv[1] == "pourcentage":
        analyse_tab_pourc(input_file=sys.argv[2], output_name=sys.argv[3])

# clustermap
# python python_scripts/heatmap.py clustermap table.csv ward Bacterie_cpd_clust tab10

# python python_scripts/heatmap.py clustermap long_read_solo_alg_mat.csv ward Long_Solo_cpd tab10
# python python_scripts/heatmap.py clustermap short_read_solo_alg_mat.csv ward Short_Solo_cpd tab10

# pourcentage
# python python_scripts/heatmap.py pourcentage All_solo_stats.csv All_solo_stats
