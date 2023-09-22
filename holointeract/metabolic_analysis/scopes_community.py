import os
import shutil
import subprocess as sub

from glob import glob
from metage2metabo.m2m.m2m_workflow import metacom_analysis
from holointeract.utils.utils import *

SOLO_METHOD = 'solo'
COOP_METHOD = 'coop'
SCOPES_STR = 'scopes'


def community_scopes(community_sbml_path, hosts_sbml_path, output_dir, seeds, method, cpu=1):
    temp_path = os.path.join(output_dir, '.temp')
    create_new_dir(temp_path)
    output_dir = os.path.join(output_dir, SCOPES_STR)
    create_new_dir(output_dir)
    output_dir = os.path.join(output_dir, method)
    create_new_dir(output_dir)

    hosts_sbml = os.listdir(hosts_sbml_path)

    for h_sbml in hosts_sbml:
        host_name = h_sbml.split('.')[0]
        h_comm_path = os.path.join(community_sbml_path, host_name)
        h_sbml_path = os.path.join(hosts_sbml_path, h_sbml)

        if method == SOLO_METHOD:
            comm_sbml = os.listdir(h_comm_path)
            for c_sbml in comm_sbml:
                c_name = c_sbml.split('.')[0]
                c_sbml_path = os.path.join(h_comm_path, c_sbml)
                h_c_output = os.path.join(output_dir, f'{host_name}__{c_name}')
                tmp_com_path = duplicate_networks_temp(c_sbml_path, temp_path, c_name)
                metacom_analysis(tmp_com_path, h_c_output, seeds, h_sbml_path, None, cpu)
                shutil.rmtree(tmp_com_path)

        elif method == COOP_METHOD:
            h_output = os.path.join(output_dir, host_name)
            create_new_dir(h_output)
            if len(os.listdir(h_comm_path)) == 1:
                c_sbml = os.listdir(h_comm_path)[0]
                c_name = c_sbml.split('.')[0]
                c_sbml_path = os.path.join(h_comm_path, c_sbml)
                tmp_com_path = duplicate_networks_temp(c_sbml_path, temp_path, c_name)
                metacom_analysis(tmp_com_path, h_output, seeds, h_sbml_path, None, cpu)
                shutil.rmtree(tmp_com_path)
            else:
                metacom_analysis(h_comm_path, h_output, seeds, h_sbml_path, None, cpu)

    shutil.rmtree(temp_path)


def duplicate_networks_temp(com_sbml_path, temp_path, c_name):
    input_dir = os.path.join(temp_path, 'input')
    create_new_dir(input_dir)
    copy_1 = os.path.join(input_dir, f'{c_name}_1.sbml')
    copy_2 = os.path.join(input_dir, f'{c_name}_2.sbml')
    shutil.copy(com_sbml_path, copy_1)
    shutil.copy(com_sbml_path, copy_2)
    return input_dir


# COEVOLUTION SCOPES
# ======================================================================================================================

def get_all_full_scope(list_algue: str, path_all_bact: str, path_sbml_algue, output_dir: str, seeds: str):
    """Compute scopes of each host with each bacteria

    Args:
        list_algue (str): list of hosts
        path_all_bact (str): path to directory of all bacteria metabolic networks
        path_sbml_algue (_type_): path to host metabolic networks
        output_dir (str): path to all scopes
        seeds (str): seeds file
    """
    bacts = glob(path_all_bact+"*")

    for algue in list_algue:
        sub.run([f'mkdir {output_dir}{algue.split("/")[-1]}'], shell=True)
        for bact in bacts:
            output = output_dir + \
                algue.split("/")[-1] + "/" + bact.split("-")[1]

            net_alg = path_sbml_algue + algue.split('/')[-1]+".sbml"

            sub.run(
                [f'mkdir {output_dir}{algue.split("/")[-1]}/{bact.split("-")[1]}'], shell=True)
            sub.run(
                [f'm2m metacom -n {bact} -s {seeds} -m {net_alg} -o {output} -c 30'], shell=True)


def coev_scopes(path_metabolic_networks, path_all_bact, list_algues, path_sbml_algues, output_dir_scopes, seeds):
    """create folder with all bacteria metabolic networks and computes scopes of each host with each bacteria

    Args:
        path_metabolic_networks (_type_): _description_
        path_all_bact (_type_): _description_
        list_algues (_type_): _description_
        path_sbml_algues (_type_): _description_
        output_dir_scopes (_type_): _description_
        seeds (_type_): _description_
    """
    get_all_full_scope(list_algues, path_all_bact,
                       path_sbml_algues, output_dir_scopes, seeds=seeds)
