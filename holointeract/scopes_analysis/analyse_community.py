from glob import glob
import subprocess as sub
import sys


def test_analyse(path: str, method: str, path_output: str, path_algues_network: str, seeds: str):
    """Compute scopes analysis of microbiota with or without bacteria cooperation depending on the 'method' argument

    Args:
        path (str): Path to host names containing a sbml folder containing bacteria metabolic networks
        method (str): 'solo' or 'coop' (if coop, bacteria of a single microbiota can cooperate)
        path_output (str): bacteria scope directory
        path_algues_network (str): directory of host metabolic networks
        seeds (str): seeds file
    """
    if method == "solo":

        algues = sorted(glob(path+"*"))
        for algue in algues:
            sub.run([f'mkdir {path_output}{algue.split("/")[-1]}'], shell=True)

            sub.run([f'echo {algue.split("/")[-1]}'], shell=True)

            bacts = glob(algue+'/sbml/*')
            for bact in bacts:
                sub.run(
                    [f'mkdir {path_output}{algue.split("/")[-1]}/{bact.split("/")[-1].split(".")[0]}'], shell=True)

                sub.run([f'mkdir ./tempo'], shell=True)
                sub.run(
                    [f'cp {algue}/sbml/{bact.split("/")[-1].split(".")[0]}* ./tempo'], shell=True)
                sub.run(
                    [f'cp {algue}/sbml/{bact.split("/")[-1].split(".")[0]}* ./tempo/{bact.split("/")[-1].split(".")[0]}_1.sbml'], shell=True)

                output = path_output + \
                    algue.split("/")[-1]+'/'+bact.split("/")[-1].split(".")[0]
                net_alg = path_algues_network + algue.split('/')[-1] + ".sbml"

                sub.run(
                    [f'm2m metacom -n ./tempo/ -s {seeds} -m {net_alg} -o {output} -c 30'], shell=True)

                sub.run([f'rm -r ./tempo'], shell=True)

    elif method == "coop":

        algues = sorted(glob(path+"*"))

        for algue in algues:
            sub.run([f'mkdir {path_output}{algue.split("/")[-1]}'], shell=True)

            sub.run([f'echo {algue.split("/")[-1]}'], shell=True)

            sub.run([f'mkdir {path_output}{algue.split("/")[-1]}'], shell=True)

            sub.run([f'mkdir ./tempo'], shell=True)

            if len(glob(f'{algue}/sbml/*')) == 1:
                file = glob(f'{algue}/sbml/*')[0]
                sub.run([f'cp {file} ./tempo'], shell=True)
                sub.run(
                    [f'cp {file} ./tempo/{file.split("/")[-1].split(".")[0]}_1.sbml'], shell=True)
            else:
                sub.run(
                    [f'cp {algue}/sbml/* ./tempo'], shell=True)

            output = path_output + algue.split("/")[-1]
            net_alg = path_algues_network + algue.split('/')[-1] + ".sbml"

            sub.run(
                [f'm2m_c.sif m2m metacom -n ./tempo/ -s {seeds} -m {net_alg} -o {output} -c 30'], shell=True)

            sub.run([f'rm -r ./tempo'], shell=True)


def prepare_all_solo_bact(input_path: str, output_dir: str):
    """Prepare the folder with all the bacteria to compute scopes of each host with each bacteria

    Args:
        input_path (str): Path to host names containing a sbml folder containing bacteria metabolic networks 
        output_dir (str): directory of all bacteria metabolic networks
    """
    algues = sorted(glob(input_path+"*"))
    all_bact = []
    for algue in algues:
        for bact in sorted(glob(algue+"/sbml/*")):
            all_bact.append(algue.split(
                "/")[-1]+"-"+bact.split("/")[-1].split(".")[0])

    sub.run([f'mkdir {output_dir}'], shell=True)
    for couple in all_bact:
        sub.run([f'mkdir {output_dir}{couple}'], shell=True)
        sub.run(
            [f'cp {input_path}{couple.split("-")[0]}/sbml/{couple.split("-")[1]}.sbml {output_dir}{couple}/{couple.split("-")[1]}.sbml'], shell=True)
        sub.run(
            [f'cp {input_path}{couple.split("-")[0]}/sbml/{couple.split("-")[1]}.sbml {output_dir}{couple}/{couple.split("-")[1]}_1.sbml'], shell=True)


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
    prepare_all_solo_bact(path_metabolic_networks, path_all_bact)
    get_all_full_scope(list_algues, path_all_bact,
                       path_sbml_algues, output_dir_scopes, seeds=seeds)


if __name__ == "__main__":

    test_analyse(path=sys.argv[1], method=sys.argv[2],
                 path_output=sys.argv[3], seeds=sys.argv[4])
