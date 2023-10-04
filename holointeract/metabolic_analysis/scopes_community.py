import os
import shutil
from metage2metabo.m2m.m2m_workflow import metacom_analysis
from holointeract.utils.utils import *
from typing import Dict


# CONSTANT STR
# ======================================================================================================================
SOLO_METHOD = 'solo'
COOP_METHOD = 'coop'
FULL_METHOD = 'full'
SCOPES_STR = 'scopes'


# MAIN FUNCTION
# ======================================================================================================================
def community_scopes(community_sbml_path: str, hosts_sbml_path: str, output_dir: str, seeds: str, method: str,
                     name_assoc: Dict[str, str], cpu: int = 1):
    """ Performs scopes for hosts and communities of microorganisms. Pairs selected depending on method:
        - coop: each host with all of its community microorganisms (#host scopes : linear)
        - solo: each host with each of its community microorganism (#microorganism scopes : linear)
        - full: each host with each microorganism (its own + other host's) (#host x #microorganism scopes : quadratic)

    Parameters
    ----------
    community_sbml_path: str
        Path to the directory containing SBML networks for all community of microorganisms classified in their natural
        host directory
    hosts_sbml_path: str
        Path to the directory containing SBML networks for all hosts
    output_dir: str
        Output directory to store the scopes files results.
    seeds: str
        Path to SBML seeds (growth medium compounds) file for performing scopes
    method: str
        Method of the scope analyse:
            - coop: each host with all of its community microorganisms (#host scopes)
            - solo: each host with each of its community microorganism (#microorganism scopes)
            - full: each host with each microorganism (its own + other host's) (#host x #microorganism scopes)
    name_assoc: Dict[str, str]
        Dictionary associating for each host and each microorganism name, its abbreviated name
    cpu: int, optional (default=1)
        Number of cpu for metacom analyse (useful for coop, not for solo and full)
    """
    are_inputs_ok(community_sbml_path, hosts_sbml_path)
    # Create directories
    temp_path = os.path.join(output_dir, '.temp')
    create_new_dir(temp_path, verbose=False)
    output_dir = os.path.join(output_dir, SCOPES_STR)
    create_new_dir(output_dir)
    output_dir = os.path.join(output_dir, method)
    create_new_dir(output_dir)
    # Get host sbml files
    hosts_sbml = os.listdir(hosts_sbml_path)
    for h_sbml in hosts_sbml:
        host_name = h_sbml.split('.')[0]
        h_comm_path = os.path.join(community_sbml_path, host_name)
        h_sbml_path = os.path.join(hosts_sbml_path, h_sbml)

        # Select method
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
    # Remove temporary files
    shutil.rmtree(temp_path)


# FUNCTIONS
# ======================================================================================================================

def coop_method_scopes(output_dir: str, host_name: str, h_comm_path: str, temp_path: str, seeds: str, h_sbml_path: str,
                       name_assoc: Dict[str, str], cpu: int):
    """ Performs scopes for coop method : each host with all of its community microorganisms (#host scopes : linear)

    Parameters
    ----------
    output_dir: str
        Output directory to store the scopes files results.
    host_name: str
        Name of the host
    h_comm_path: str
        Path to the host's community SBML files
    temp_path: str
        Path to the temporary files directory
    seeds: str
        Path to SBML seeds (growth medium compounds) file for performing scopes
    h_sbml_path: str
        Path to the host's SBML files
    name_assoc: Dict[str, str]
        Dictionary associating for each host and each microorganism name, its abbreviated name
    cpu: int, optional
        Number of cpu for metacom analyse
    """
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


def solo_method_scopes(h_comm_path: str, output_dir: str, host_name: str, temp_path: str, seeds: str, h_sbml_path: str,
                       name_assoc: Dict[str, str], cpu: int):
    """ Performs scopes for solo method : each host with each of its community microorganism (#microorganism scopes :
    linear)

    Parameters
    ----------
    h_comm_path: str
        Path to the host's community SBML files
    output_dir: str
        Output directory to store the scopes files results.
    host_name: str
        Name of the host
    temp_path: str
        Path to the temporary files directory
    seeds: str
        Path to SBML seeds (growth medium compounds) file for performing scopes
    h_sbml_path: str
        Path to the host's SBML files
    name_assoc: Dict[str, str]
        Dictionary associating for each host and each microorganism name, its abbreviated name
    cpu: int, optional
        Number of cpu for metacom analyse
    """
    comm_sbml = os.listdir(h_comm_path)
    for c_sbml in comm_sbml:
        c_name = c_sbml.split('.')[0]
        c_sbml_path = os.path.join(h_comm_path, c_sbml)
        h_c_output = os.path.join(output_dir, name_assoc[f'{host_name}_{c_name}'])
        tmp_com_path = duplicate_networks_temp(c_sbml_path, temp_path, c_name)
        metacom_analysis(tmp_com_path, h_c_output, seeds, h_sbml_path, None, cpu)
        shutil.rmtree(tmp_com_path)


def full_method_scopes(h_comm_path: str, output_dir: str, host_name: str, temp_path: str, seeds: str, h_sbml_path: str,
                       name_assoc: Dict[str, str], comm_host: str,  cpu: int):
    """ Performs scopes for full method : each host with each microorganism (its own + other host's)
    (#host x #microorganism scopes : quadratic)

    Parameters
    ----------
    h_comm_path: str
        Path to the host's community SBML files
    output_dir: str
        Output directory to store the scopes files results.
    host_name: str
        Name of the host
    temp_path: str
        Path to the temporary files directory
    seeds: str
        Path to SBML seeds (growth medium compounds) file for performing scopes
    h_sbml_path: str
        Path to the host's SBML files
    name_assoc: Dict[str, str]
        Dictionary associating for each host and each microorganism name, its abbreviated name
    comm_host: str
        Name of the natural microorganism host
    cpu: int, optional
        Number of cpu for metacom analyse
    """
    comm_sbml = os.listdir(h_comm_path)
    for c_sbml in comm_sbml:
        c_name = c_sbml.split('.')[0]
        c_sbml_path = os.path.join(h_comm_path, c_sbml)
        h_c_output = os.path.join(output_dir, name_assoc[host_name] + '_' + name_assoc[f'{comm_host}_{c_name}'])
        tmp_com_path = duplicate_networks_temp(c_sbml_path, temp_path, c_name)
        metacom_analysis(tmp_com_path, h_c_output, seeds, h_sbml_path, None, cpu)
        shutil.rmtree(tmp_com_path)


def duplicate_networks_temp(microorganism_sbml_path: str, temp_path: str, microorganism_name: str) -> str:
    """ Duplicate a microorganism sbml network file in the temporary directory. Used in case of metacom analyse between
    1 host and exactly 1 microorganism. Metacom allowing only community size > 1.

    Parameters
    ----------
    microorganism_sbml_path: str
        Path to the SBML network file for a microorganism
    temp_path: str
        Path to the temporary files directory
    microorganism_name: str
        Name of the microorganism to duplicate sbml network file

    Returns
    -------
    str
        Path to the input directory for metacom containing the duplicated sbml networks files
    """
    input_dir = os.path.join(temp_path, 'input')
    create_new_dir(input_dir, verbose=False)
    copy_1 = os.path.join(input_dir, f'{microorganism_name}_1.sbml')
    copy_2 = os.path.join(input_dir, f'{microorganism_name}_2.sbml')
    shutil.copy(microorganism_sbml_path, copy_1)
    shutil.copy(microorganism_sbml_path, copy_2)
    return input_dir
