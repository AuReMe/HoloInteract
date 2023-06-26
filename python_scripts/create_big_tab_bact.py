from glob import glob
import subprocess as sub
import sys
import json




def get_full_metabo_list(path_alg, path_bact, path_holo):
    """Parcours les résultats m2m metacom et va ressortir la liste de tous les métabolites ainsi que les métabolites
    produits par les algues, les bactéries et les holobiontes

    Args:
        path_alg (str): chemin vers les résultats de scope (comm_scope)
        path_bact (str): chemin vers les added value
        path_holo (str): chemin vers les added value
        method (str, optional): Auttorise interactions entre bactéries ou pas ('solo' si non). Defaults to "coop".

    Returns:
        list(list_metabo) (_type_): liste de tous les metabolites retrouves dans les scopes
        dico_alg_metabo (_type_): dictionnaire des métabolites associés à chaque algue
        dico_bact_metabo (_type_): dictionnaire des métabolites associés à chaque bactérie
        dico_holo_metabo (_type_): dictionnaire des métabolites associés à chaque holobionte (added value)
    """
    dico_alg_metabo = dict()
    list_metabo=set()
    dico_bact_metabo = dict()
    dico_holo_metabo = dict()
    print("start",len(list_metabo))

    files = glob(path_alg+"*")
    algues = [x for x in files if x not in ["m2m_recon.log","pgdb", "recon_stats.tsv"]]
    #algue
    for i in sorted(algues):
        scope_file = i+"/community_analysis/comm_scopes.json"
        with open(scope_file, "r") as fi:
            data = json.load(fi)
        for value in data.values():
            list_metabo.update(set(value))
        dico_alg_metabo[i.split("/")[-1]] = data["host_scope"]


    #bact

    for i in sorted(glob(path_bact+"*")):  
        for j in sorted(glob(i+"/*")):
            if len(glob(j+"/indiv_scopes/*")) >0:
                with open(j+"/indiv_scopes/indiv_scopes.json") as fi:
                    data = json.load(fi)
                for value in data.values():
                        list_metabo.update(set(value))
                dico_bact_metabo[i.split("/")[-1]+"-"+j.split("/")[-1]] = data[j.split("/")[-1]]
            else:
                print("probleme with :", i)

    #holobionte
    for i in sorted(glob(path_holo+"*")):
            for j in sorted(glob(i+"/*")):
                if len(sorted(glob(j+"/community_analysis/*"))) >1:
                    with open(j+"/community_analysis/addedvalue.json") as fi:
                        data = json.load(fi)
                    for value in data.values():
                        list_metabo.update(set(value))
                    dico_holo_metabo[i.split("/")[-1]+"-"+j.split("/")[-1]] = data['addedvalue']
                else:
                    print("probleme with :", i)


    return list(list_metabo), dico_alg_metabo, dico_bact_metabo, dico_holo_metabo

def associate_alg_bact(path):

    list_holo = []
    algaes = sorted(glob(path+"*"))
    dico_presence = dict()
    liste_algues = set()
    for i in algaes:
        if i.split("/")[-1] not in dico_presence:
            bacteries = [j.split("/")[-1].split(".")[0] for j in sorted(glob(i+"/sbml/*")) if j.split("/")[-1] != "taxon_id.tsv"]
            dico_presence[i.split("/")[-1]] = bacteries

    nb_columns = 0
    for i in dico_presence:
        list_holo.extend([i+"-"+x for x in dico_presence[i]])
        nb_columns +=1

    return dico_presence,list_holo

def matrix(list_metabo, list_holo,dico_algue,dico_bact,dico_holo):

    # Creation de matrice
    matrice = [["" for _ in range(len(list_holo)+1)]for _ in range((len(list_metabo)-1))] 
    list_holo.insert(0," ")
    matrice.insert(0, list_holo)

    for i in range(len(matrice)):
        if i>0:
            matrice[i][0] = list_metabo[i-1]


    for i in range(len(matrice)):
        for j in range(len(matrice[0])):
            if i> 0 and j>0 and matrice[0][j] in dico_bact and matrice[0][j] in dico_holo:
                if matrice[i][0] in dico_bact[matrice[0][j]]:
                    matrice[i][j] += "2"
                if matrice[i][0] in dico_holo[matrice[0][j]]:
                    matrice[i][j] += "3"

                if matrice[i][0] in dico_algue[matrice[0][j].split("-")[0]]:
                    matrice[i][j] += "1"
    for i in range(len(matrice)):
        for j in range(len(matrice[0])):
            if matrice[i][j] == "":
                matrice[i][j] = "0"
            if matrice[i][j] == "21":
                matrice[i][j] = "4"
    matrix_added = []
    for i in matrice:
        if "3" in i[1:] or i == matrice[0]:
            matrix_added.append(i)

    return matrice, matrix_added

            

def modif_cpd_name(matrix):
    for i in range(len(matrix)):
        name = matrix[i][0]
        if "__" in matrix[i][0]:
            name = matrix[i][0].split("__")
            for j in range(len(name)):
                if name[j] == "45":
                    name[j] = "-"
                elif name[j] == "46":
                    name[j] = "."
                elif name[j] == "43":
                    name[j] = "+"
            name = "".join(name)
        matrix[i][0] = name[2:-2]
    return matrix

def write_mat(mat, name_mat):
    with open(name_mat+".csv", "w")as fo:
                for i in mat:
                    fo.write(";".join(i)+"\n")
    

def job(path_alg, path_bact, path_holo, path_metabolic_networks, output_name):
    list_holo = ""
    metabos, dico_algue, dico_bact, dico_holo= get_full_metabo_list(path_alg, path_bact, path_holo)

    dico, list_holo= associate_alg_bact(path_metabolic_networks)
    mat,mat_added = matrix( metabos, list_holo, dico_algue, dico_bact, dico_holo)
    new_mat = modif_cpd_name(mat)
            
    write_mat(new_mat, output_name)
    return dico_algue, dico_bact

if __name__ == "__main__":
    #argv1 : chemin vers scope des algues
    #argv2 : chemin vers scope des bactérie (jusqu'a /Added_value)
    #argv3 : chemin vers SBML_EggNog
    job(sys.argv[1], sys.argv[2], sys.argv[2], sys.argv[3],sys.argv[4])

    

# python python_scripts/create_big_tab.py /scratch/clucas/algaes_scopes/ /scratch/clucas/Added_value/ /groups/phaeo_bact/SBML_EggNog/
