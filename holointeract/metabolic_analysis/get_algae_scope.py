from glob import glob
import subprocess as sub
import sys
from holointeract.metabolic_analysis.create_big_tab_alg import modif_cpd_name
import pandas


def run_miscoto(path: str, seeds: str, output: str):
    """Get hosts scope

    Args:
        path (str): path to host metabolic networks
        seeds (str): Seeds files
        output (str): Scope directory
    """

    algaes = sorted(glob(path+"*"))
    sub.run([f'mkdir {output}'], shell=True)
    for algae in algaes:
        out = output + algae.split("/")[-1].split(".")[0]
        sub.run([f'mkdir ./tempo'], shell=True)
        sub.run([f'cp {algae} ./tempo/'], shell=True)
        sub.run(
            [f'cp {algae} ./tempo/{algae.split("/")[-1].split(".")[0]}_1.sbml'], shell=True)
        sub.run(
            [f'm2m metacom -n ./tempo/ -s {seeds} -m {algae} -o {out} -c 30'], shell=True)
        sub.run([f'rm -r ./tempo'], shell=True)
