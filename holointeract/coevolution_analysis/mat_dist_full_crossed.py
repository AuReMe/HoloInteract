import os
from glob import glob
import json
import pandas as pd
import matplotlib.pyplot as plt


def create_matrix(host_path: str, comm_path: str):
    """
    """
    host_list = [x.split('.')[0] for x in os.listdir(host_path)]
    comm_list = []
    for comm_host in os.listdir(comm_path):
        comm_list += [x.split('.')[0] for x in os.listdir(os.path.join(comm_path, comm_host))]

    print(host_list)
    print(comm_list)

    # matrix = [["" for _ in range(len(list_algue)+1)]
    #           for _ in range((len(list_bact)+1))]
    #
    # list_algue.insert(0, " ")
    #
    # matrix[0] = list_algue
    # print(matrix)
    # print(list_bact)
    # for i in range(len(matrix)):
    #     if i > 0:
    #         print(list_bact[i-1].split("/")[-1])
    #         matrix[i][0] = list_bact[i-1].split("/")[-1]


COMM_PATH = '../../example/inputs/community'
HOST_PATH = '../../example/inputs/hosts'
create_matrix(HOST_PATH, COMM_PATH)


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


def job(list_algue: list, list_bact: list, scopes_bacteries_path: str, output_name: str = "matrice_complementarite_normalized"):
    """Create the csv file of the scopes of each host with each bacteria

    Args:
        list_algue (list): list of hosts
        list_bact (list): list of bacteria
        scopes_bacteries_path (str): directory of each host with each bacteria
        output_name (str, optional): output name for csv file . Defaults to "matrice_complementarite_normalized".
    """
    mat = create_matrix(list_algue, list_bact)

    data = fill_mat(mat, scopes_bacteries_path)

    data = clean_df(data)

    data = normalisation(data)

    data.to_csv(output_name+".csv")

    boxploter(data, output_name)
