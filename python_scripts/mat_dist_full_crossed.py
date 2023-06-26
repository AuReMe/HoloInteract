from glob import glob
import subprocess as sub
import sys
import json
import pandas as pd
from scipy.stats import ttest_ind
import matplotlib.pyplot as plt
import pprint


def create_matrix(list_algue: str, list_bact: str):
    """_summary_

    Args:
        path_alg (str): _description_
        path_bact (str): _description_

    Returns:
        _type_: _description_
    """

    matrix = [["" for _ in range(len(list_algue)+1)]
              for _ in range((len(list_bact)+1))]

    list_algue.insert(0, " ")

    matrix[0] = list_algue
    print(matrix)
    print(list_bact)
    for i in range(len(matrix)):
        if i > 0:
            print(list_bact[i-1].split("/")[-1])
            matrix[i][0] = list_bact[i-1].split("/")[-1]
    return matrix


def fill_mat(matrix: list, path: str):
    """_summary_

    Args:
        matrix (list): _description_
        path (str): _description_
        path_fullmicrobiota (str): _description_

    Returns:
        _type_: _description_
    """

    for col in range(len(matrix[0])):
        for line in range(len(matrix)):
            if col > 0 and line > 0 and len(glob(path+matrix[0][col]+"/"+matrix[line][0].split("-")[-1]+"/community_analysis/*")) > 1:
                with open(path+matrix[0][col]+"/"+matrix[line][0].split("-")[-1]+"/community_analysis/addedvalue.json") as fi:
                    data = json.load(fi)
                matrix[line][col] = str(len(data["addedvalue"]))

    return matrix


def clean_df(matrix):
    index = [x[0] for x in matrix]
    data = pd.DataFrame(matrix[1:], columns=matrix[0], index=index[1:])
    same_cols = []

    # Boucle pour comparer chaque colonne avec toutes les autres colonnes
    for i, col in enumerate(data.columns):
        # Vérification si la colonne a déjà été ajoutée à la liste des colonnes identiques
        if col not in same_cols:
            # Ajout de la colonne à la liste des colonnes identiques
            # Comparaison de la colonne avec toutes les autres colonnes
            for j in range(i+1, len(data.columns)):
                if data[col].equals(data.iloc[:, j]):
                    # Ajout des colonnes identiques à la liste des colonnes identiques
                    if data.columns[j] not in same_cols:
                        same_cols.append(data.columns[j])

    # Affichage de la liste des colonnes qui ont les mêmes valeurs
    data = data.drop(same_cols, axis=1)
    indices_to_remove = []

    for i in same_cols:
        indices_to_remove = data[data.index.str.startswith(i)].index
        data = data.drop(indices_to_remove)
    #
    data = data.set_index(" ")
    data.index
    data = data[~(data == "").all(axis=1)]
    #data = data[~(data == "0").all(axis=1)]

    #data = data.reset_index(drop=True)
    data = data[data.columns].astype(float)

    return data


def normalisation(data):
    for col in data.columns:
        if col != " ":
            mean_col = data[col].mean()

            data[col] = data[col]/mean_col
    return data


def boxploter(data, output_name):
    # Sélection des colonnes qui ont toutes les mêmes valeurs
    # Création d'une liste pour stocker tous les noms de colonnes qui ont les mêmes valeurs
    # Création du boxplot
    fig, ax = plt.subplots(figsize=(10, 5))
    bp = ax.boxplot(data.values, positions=range(len(data.columns)))

    # Configuration de l'axe X
    ax.set_xticks(range(len(data.columns)))

    ax.set_xticklabels(data.columns, rotation=90)
    plt.subplots_adjust(top=0.95, bottom=0.25)
    #plt.title("Comparaison des données de chaque algue")
    plt.ylabel("Complémentarité normalisée")
    # Affichage du graphique
    plt.savefig(output_name+".png")


def job(list_algue, list_bact, scopes_bacteries_path, output_name="matrice_complementarite_normalized"):
    mat = create_matrix(list_algue, list_bact)

    data = fill_mat(mat, scopes_bacteries_path)

    data = clean_df(data)

    data = normalisation(data)

    data.to_csv(output_name+".csv")

    boxploter(data, output_name)


if __name__ == "__main__":

    job(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

# python python_scripts/mat_dist_full_crossed.py /groups/phaeo_bact/SBML_EggNog/ /scratch/clucas/all_scopes_matrix/ /scratch/clucas/scopes_full_microbiota/ mat_dist_full_crossed_nonorm

# python python_scripts/mat_dist_full_crossed.py /groups/phaeo_bact/Long_reads/ /scratch/clucas/Long_all_scopes/ /scratch/clucas/Long_full_scopes/ Long_dist_full_crossed
# python python_scripts/mat_dist_full_crossed.py /groups/phaeo_bact/Short_reads/ /scratch/clucas/Short_all_scopes/ /scratch/clucas/Short_full_scopes/ Short_dist_full_crossed
