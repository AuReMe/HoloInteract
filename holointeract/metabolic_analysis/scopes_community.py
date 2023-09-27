import os
import shutil

from metage2metabo.m2m.m2m_workflow import metacom_analysis
from holointeract.utils.utils import *

SOLO_METHOD = 'solo'
COOP_METHOD = 'coop'
FULL_METHOD = 'full'
SCOPES_STR = 'scopes'


def community_scopes(community_sbml_path, hosts_sbml_path, output_dir, seeds, method, name_assoc, cpu=1):
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
            solo_method_scopes(h_comm_path, output_dir, host_name, temp_path, seeds, h_sbml_path, name_assoc, cpu)

        elif method == COOP_METHOD:
            coop_method_scopes(output_dir, host_name, h_comm_path, temp_path, seeds, h_sbml_path, name_assoc, cpu)

        elif method == FULL_METHOD:
            all_comm_path = os.listdir(community_sbml_path)
            for comm_host in all_comm_path:
                comm_path = os.path.join(community_sbml_path, comm_host)
                full_method_scopes(comm_path, output_dir, host_name, temp_path, seeds, h_sbml_path, name_assoc,
                                   comm_host, cpu)

    shutil.rmtree(temp_path)


def coop_method_scopes(output_dir, host_name, h_comm_path, temp_path, seeds, h_sbml_path, name_assoc, cpu):
    h_output = os.path.join(output_dir, name_assoc[host_name])
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


def solo_method_scopes(h_comm_path, output_dir, host_name, temp_path, seeds, h_sbml_path, name_assoc, cpu):
    comm_sbml = os.listdir(h_comm_path)
    for c_sbml in comm_sbml:
        c_name = c_sbml.split('.')[0]
        c_sbml_path = os.path.join(h_comm_path, c_sbml)
        h_c_output = os.path.join(output_dir, name_assoc[f'{host_name}_{c_name}'])
        tmp_com_path = duplicate_networks_temp(c_sbml_path, temp_path, c_name)
        metacom_analysis(tmp_com_path, h_c_output, seeds, h_sbml_path, None, cpu)
        shutil.rmtree(tmp_com_path)


def full_method_scopes(h_comm_path, output_dir, host_name, temp_path, seeds, h_sbml_path, name_assoc, comm_host,  cpu):
    comm_sbml = os.listdir(h_comm_path)
    for c_sbml in comm_sbml:
        c_name = c_sbml.split('.')[0]
        c_sbml_path = os.path.join(h_comm_path, c_sbml)
        h_c_output = os.path.join(output_dir, name_assoc[host_name] + '_' + name_assoc[f'{comm_host}_{c_name}'])
        tmp_com_path = duplicate_networks_temp(c_sbml_path, temp_path, c_name)
        metacom_analysis(tmp_com_path, h_c_output, seeds, h_sbml_path, None, cpu)
        shutil.rmtree(tmp_com_path)


def duplicate_networks_temp(com_sbml_path, temp_path, c_name):
    input_dir = os.path.join(temp_path, 'input')
    create_new_dir(input_dir)
    copy_1 = os.path.join(input_dir, f'{c_name}_1.sbml')
    copy_2 = os.path.join(input_dir, f'{c_name}_2.sbml')
    shutil.copy(com_sbml_path, copy_1)
    shutil.copy(com_sbml_path, copy_2)
    return input_dir

