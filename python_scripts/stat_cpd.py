import matplotlib.pyplot as plt
from python_scripts.onto_metacyc import get_classes, get_classes_opti
import pandas as pd
import plotly.graph_objs as go
import seaborn as sns
import sys
import matplotlib.pyplot as plt


def count_cat(table: str):
    """Count occurences of compounds for each type of production

    Args:
        table (str): csv file

    Returns:
        dict: dictionnary
    """
    dico = dict()
    matrix = []
    with open(table, 'r') as fi:
        for line in fi:
            matrix.append(line.strip().split(";"))

    for line in range(len(matrix)):
        if matrix[line][0] != "":
            dico[matrix[line][0]] = {"1": 0, "2": 0, "3": 0, "4": 0, "5": 0}
            for col in range(len(matrix[line])):
                if col > 0 and matrix[line][col] in dico[matrix[line][0]]:
                    dico[matrix[line][0]][matrix[line][col]] += 1

    return dico


def hist_global(dico: dict, output_fig_cpd):
    """histogramme des molecules

    Args:
        dico (dict): dico de la fonction count_cat
    """
    compounds = []
    valeur = []
    for cpd in dico:
        if dico[cpd]["3"] >= (sum(dico[cpd].values())/3):
            compounds.append(cpd)
            valeur.append(dico[cpd]["3"])

    # trier les catégories par ordre décroissant
    sorted_compounds, sorted_valeur = zip(
        *sorted(zip(compounds, valeur), key=lambda x: x[1], reverse=False))

    fig = go.Figure(
        [go.Bar(x=sorted_valeur, y=sorted_compounds, orientation='h')])
    fig.update_layout(
        title="Histogramme des molécules",
        xaxis_title="Fréquence",
        yaxis_title="Catégories"
    )
    fig.write_image(output_fig_cpd+".png")


def write_tab(dico: dict, output: str):
    with open(output+".csv", "w") as fo:
        fo.write(f"Compound;Lien Metacyc;Nb holobiontes;Ontology\n")
    data = ""
    for cpd in dico:
        print(cpd)
        if dico[cpd]["3"] > 0:
            data += f'{cpd};https://metacyc.org/compound?orgid=META&id={cpd};{dico[cpd]["3"]};\
                {get_classes(cpd, "/scratch/clucas/HoloInteract/toy_example/metacyc_26.0.padmet")}\n'
    with open(output+".csv", "a") as fo:
        fo.write(data)


def write_tab_opti(dico: dict, output: str, padmet: str):
    with open(output+".csv", "w") as fo:
        fo.write(f"Compound;Lien Metacyc;Nb holobiontes;Ontology\n")
    data = ""
    dico_onto = dict()
    with open(padmet, "r") as fi:
        for line in fi:
            if "is_a_class" in line:
                liste = line.strip().split("	")
                dico_onto[liste[0]] = liste[2]

    for cpd in dico:
        print(":)", cpd)
        if dico[cpd]["3"] > 0:
            data += f'{cpd};https://metacyc.org/compound?orgid=META&id={cpd};{dico[cpd]["3"]};\
                {get_classes_opti(cpd, dico_onto)}\n'
    with open(output+".csv", "a") as fo:
        fo.write(data)


def classing_cpd(input_file, output_file):
    with open(input_file) as f:
        lines = f.readlines()

    with open(output_file, 'w') as out:
        out.write("ID;0;1;2;3;4;5\n")  # En-tête du fichier de sortie
        for line in lines:
            values = line.strip().split(';')
            categories = values[1:]
            category_counts = {}
            for category in categories:
                if category in category_counts:
                    category_counts[category] += 1
                else:
                    category_counts[category] = 1
            total_categories = len(categories)
            proportions = {}
            for i in range(6):
                if str(i) in category_counts:
                    proportions[str(i)] = category_counts[str(i)
                                                          ] / total_categories * 100
                else:
                    proportions[str(i)] = 0
            out.write(
                f"{values[0]};{proportions['0']:.2f};{proportions['1']:.2f};{proportions['2']:.2f};{proportions['3']:.2f};{proportions['4']:.2f};{proportions['5']:.2f}\n")


def mole_count(input_file):
    # Charger les données à partir d'un fichier CSV
    data = pd.read_csv(input_file, sep=";", header=0, index_col=0)

    sns.violinplot(data=data, inner="box")
    plt.ylim(bottom=0)
    plt.title("Boxplot")
    plt.show()
    sns.boxplot(data=data)
    plt.show()


def get_category_production(categories: list, pourcentage: list, input_file: str):
    if len(categories) != len(pourcentage):
        raise TypeError

    data = pd.read_csv(input_file, sep=",", header=0, index_col=0)
    if len(categories) > 1:
        print(categories[0], pourcentage[0])
        print(categories[1], pourcentage[1])
        for index, row in data.iterrows():
            if type(pourcentage[0]) == tuple and type(pourcentage[1]) == tuple:
                if float(pourcentage[0][0]) <= float(row[categories[0]]) <= float(pourcentage[0][1]) and float(pourcentage[1][0]) <= float(row[categories[1]]) <= float(pourcentage[1][1]):
                    print(index)

            elif float(row[categories[0]]) > float(pourcentage[0]) and float(row[categories[1]]) > float(pourcentage[1]):
                print(index)
    else:
        for index, row in data.iterrows():
            if type(pourcentage[0]) == tuple:
                if float(pourcentage[0][0]) <= float(row[categories[0]]) <= float(pourcentage[0][1]):
                    print(index, row[categories[0]])

            elif float(row[categories[0]]) == float(pourcentage[0]):
                print(index, row[categories[0]], ">", pourcentage[0])


def job(matrice, output_fig_cpd, output_file_cpd, padmet):
    dico = count_cat(matrice)
    hist_global(dico, output_fig_cpd)
    write_tab_opti(dico, output_file_cpd, padmet)


if __name__ == "__main__":
    job(sys.argv[1], sys.argv[2], sys.argv[3])
