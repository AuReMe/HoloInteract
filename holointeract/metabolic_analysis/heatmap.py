import seaborn as sns
import pandas


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