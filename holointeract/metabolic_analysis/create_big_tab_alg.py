from glob import glob
import subprocess as sub
import json
import sys
import pprint




def get_full_metabo_list(path_alg, path_bact, path_holo, method = "coop"):
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

    if method == "solo":
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
    
    elif method == "coop":
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
        print("alg",len(list_metabo))

        #bact

        for i in sorted(glob(path_bact+"*")):
            if len(sorted(glob(i+"/indiv_scopes/*"))) >1:
                with open(i+"/indiv_scopes/indiv_scopes.json") as fi:
                    data = json.load(fi)
                for value in data.values():
                    list_metabo.update(set(value))
                for cle in data:
                    dico_bact_metabo[i.split("/")[-1]+"-"+cle] = data[cle]
            else:
                print("probleme with :", i)
        print("bact",len(list_metabo))


        #holobionte
        for i in sorted(glob(path_holo+"*")):
            if len(sorted(glob(i+"/community_analysis/*"))) >1:
                with open(i+"/community_analysis/addedvalue.json") as fi:
                    data = json.load(fi)
                for value in data.values():
                    list_metabo.update(set(value))
                dico_holo_metabo[i.split("/")[-1]+"-"+i.split("/")[-1]] = data['addedvalue']
            else:
                print("probleme with :", i)
        print("holo",len(list_metabo))


    return list(list_metabo), dico_alg_metabo, dico_bact_metabo, dico_holo_metabo


def associate_alg_bact(path):
    """retrouver les bacteries présentes dans le microbiote de chaque algue

    Args:
        path (str): chemin vers les réseaux métaboliques des bactéries triés par algue (Eggnog ou prokka)

    Returns:
        dict: microbiotes
        int: nombre de bactéries au total
        list: nom des holobiontes
    """
    list_holo = []
    algaes = sorted(glob(path+"*"))
    dico_presence = dict()
    for i in algaes:
        if i.split("/")[-1] not in dico_presence:
            bacteries = [j.split("/")[-1].split(".")[0] for j in sorted(glob(i+"/sbml/*")) if j.split("/")[-1] != "taxon_id.tsv"]
            dico_presence[i.split("/")[-1]] = bacteries

    nb_columns = 0
    for i in dico_presence:
        #if i != "Clin_U_DA3" and  i != "Esil_F_AA2" and i !="Fjap_U_CB1" and i != "Ecro_F_AB2":
        nb_columns += len(dico_presence[i])
        list_holo.extend([i+"-"+x for x in dico_presence[i]])

    return dico_presence, nb_columns, list_holo



def matrix_alg(path_alg, dico_alg, dico_bact, dico_holo, metabos, dico_presence, method="coop"):
    """Creer une matrice répertoriant la production de metabolites par les holobiontes

    Args:
        path_alg (_type_): chemin vers les réseaux métaboliques
        dico_alg (_type_): dictionnaire de la fonction get_full_metabo_list
        dico_bact (_type_): dictionnaire de la fonction get_full_metabo_list
        dico_holo (_type_): dictionnaire de la fonction get_full_metabo_list
        metabos (_type_): liste des metabolites
        dico_presence (_type_): dictionnaire des microbiote
        method (str, optional): Autorise les cooperation ou non entre bacteries ('solo' pour non). Defaults to "coop".

    Returns:
        list: matrice
    """
    liste_alg = []
    if method == "coop":
        for alg in sorted(glob(path_alg+"*")):
            if alg+"/sbml" in sorted(glob(alg+"/*")):
                liste_alg.append(alg.split("/")[-1])

        liste_alg.insert(0," ")
        
        matrix = [["" for _ in range(len(liste_alg)+1)]for _ in range((len(metabos)-1))] 
        matrix.insert(0, liste_alg)
        
        print("alg", dico_alg.keys())
        print("holo",dico_holo.keys())
        print("bact",dico_bact.keys())

        for i in range(len(matrix)):
            if i>0:
                matrix[i][0] = metabos[i-1]

        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                if i >0 and j>0 and matrix[i][0] in dico_alg[matrix[0][j]] :
                    if "1" not in matrix[i][j]:
                        matrix[i][j] += "1"
                if i >0 and j>0:
                    for bact in dico_presence[matrix[0][j]]:
                        if i >0 and j >0 and matrix[i][0] in dico_holo[matrix[0][j]+"-"+matrix[0][j]] or matrix[i][0] in dico_bact[matrix[0][j]+"-"+bact]:
                            #bact
                            if matrix[0][j]+"-"+bact in dico_bact and matrix[i][0] in dico_bact[matrix[0][j]+"-"+bact] and "2" not in matrix[i][j]:
                                matrix[i][j] += "2"
                            #holo
                            if matrix[0][j]+"-"+matrix[0][j] in dico_holo and matrix[i][0] in dico_holo[matrix[0][j]+"-"+matrix[0][j]] and "3" not in matrix[i][j]:
                                matrix[i][j] += "3"
                            
                            

        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                if i>0 and j > 0 and matrix[i][j] == "":
                    matrix[i][j] = "0"
                if i>0 and j > 0 and matrix[i][j] == "12":
                    matrix[i][j] = "5"
                if i > 0 and j > 0 and ("32" in matrix[i][j] or "23" in matrix[i][j]):
                    matrix[i][j] = "4"
    
    ########## Solo

    elif method == "solo":
        print("alg", dico_alg.keys())
        print("holo",dico_holo.keys())
        print("bact",dico_bact.keys())

        for alg in sorted(glob(path_alg+"*")):
            if  alg+"/sbml" in glob(alg+"/*"):
                liste_alg.append(alg.split("/")[-1])


        liste_alg.insert(0," ")
        
        matrix = [["" for _ in range(len(liste_alg)+1)]for _ in range((len(metabos)-1))] 
        matrix.insert(0, liste_alg)

        for i in range(len(matrix)):
            if i>0:
                matrix[i][0] = metabos[i-1]

        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                if i >0 and j>0 and matrix[i][0] in dico_alg[matrix[0][j]] :
                    matrix[i][j] += "1"

                if i >0 and j>0:
                    for bact in dico_presence[matrix[0][j]]:
                        if i >0 and j >0 and matrix[i][0] in dico_holo[matrix[0][j]+"-"+bact] or matrix[i][0] in dico_bact[matrix[0][j]+"-"+bact]:
                            #bact
                            if matrix[0][j]+"-"+bact in dico_bact and matrix[i][0] in dico_bact[matrix[0][j]+"-"+bact] and "2" not in matrix[i][j]:
                                matrix[i][j] += "2"
                            #holo
                            if matrix[0][j]+"-"+bact in dico_holo and matrix[i][0] in dico_holo[matrix[0][j]+"-"+bact] and "3" not in matrix[i][j]:
                                matrix[i][j] += "3"

        for i in range(len(matrix)):
            for j in range(len(matrix[0])):
                if i>0 and j > 0 and matrix[i][j] == "":
                    matrix[i][j] = "0"
                if i>0 and j > 0 and matrix[i][j] == "12":
                    matrix[i][j] = "5"
                if i > 0 and j > 0 and ("32" in matrix[i][j] or "23" in matrix[i][j]):
                    matrix[i][j] = "4"

    return matrix

def modif_cpd_name(matrix):
    """Ecris les bons noms des composes

    Args:
        matrix (list): matrice de la fonction matrix_alg

    Returns:
        list: matrice
    """
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
    """ecris la matrice dans un fichier csv

    Args:
        mat (list): maatrice
        name_mat (str)): nom du fichier à creer
    """
    with open(name_mat+".csv", "w")as fo:
        for i in mat:
            fo.write(";".join(i)[:-1]+"\n")


def job(alg_scopes, bact_scopes, sbml_path, output_name, method="coop"):
    metabos, dico_algue, dico_bact, dico_holo= get_full_metabo_list(alg_scopes, bact_scopes, bact_scopes,method=method)

    dico_presence, nb_columns, list_holo = associate_alg_bact(sbml_path)
    metabos = sorted(metabos, reverse=True)

    matrice = matrix_alg(sbml_path, dico_algue,dico_bact, dico_holo, metabos, dico_presence, method=method)
    matrice = modif_cpd_name(matrice)
    write_mat(matrice, output_name)



if __name__ == "__main__":
    #argv1 : chemin vers scope des algues
    #argv2 : chemin vers scope des bactéries (jusqu'a /Added_value)
    #argv3 : chemin vers SBML_EggNog
    #argv4 : nom du fichier de sortie
    #argv5 :  coop ou solo


    job(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],)

    """metabos, dico_algue, dico_bact, dico_holo= get_full_metabo_list(sys.argv[1], sys.argv[2], sys.argv[2],method=sys.argv[5])

    dico_presence, nb_columns, list_holo = associate_alg_bact(sys.argv[3])
    metabos = sorted(metabos, reverse=True)
    meme = []
    for i in dico_bact:
        if i not in dico_holo:
            meme.append(i)
    for i in meme:
        del dico_bact[i]

    matrice = matrix_alg(sys.argv[3], dico_algue,dico_bact, dico_holo, metabos, dico_presence, method=sys.argv[5])
    matrice = modif_cpd_name(matrice)
    write_mat(matrice, sys.argv[4])"""

# python holointeract/create_big_tab_alg.py /scratch/clucas/algaes_scopes/ /scratch/clucas/Added_value/ /groups/phaeo_bact/SBML_EggNog/ alg_mat solo
# python holointeract/create_big_tab_alg.py /scratch/clucas/algaes_scopes/ /scratch/clucas/Coop_Added_value/ /groups/phaeo_bact/SBML_EggNog/ coop_alg_mat coop
