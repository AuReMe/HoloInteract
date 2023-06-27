from glob import glob
import sys


def get_classes(cpd: str, file: str):
    """get classes f the compound 

    Args:
        cpd (str): compound
        file (str): metacyc padmet file

    Returns:
        str: classes of the compound
    """
    with open(file, "r") as fi:
        for line in fi:
            if cpd == "Compounds" and line.split("\t")[0] == cpd and "is_a_class" in line:
                return ""
            elif line.split("\t")[0] == cpd and "is_a_class" in line:
                return get_classes(line.strip().split("\t")[2], file)+"/"+line.strip().split("\t")[2]
    return ""


def get_classes_opti(cpd: str, dico: str):
    """get classes f the compound 

    Args:
        cpd (str): compound
        file (str): metacyc padmet file

    Returns:
        str: classes of the compound
    """
    if cpd.split("/")[-1] in dico and cpd.split("/")[-1] != "Compounds":
        return get_classes_opti(cpd+"/"+dico[cpd.split("/")[-1]], dico)
    return cpd


if __name__ == "__main__":
    a = get_classes("CPD-21306", "/udd/colucas/Documents/metacyc_26.0.padmet")
    print('a')
