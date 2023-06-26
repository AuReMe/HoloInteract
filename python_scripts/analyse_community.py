from glob import glob
import subprocess as sub
import sys


def test_analyse(path: str, method: str, path_output: str, path_algues_network: str, seeds: str):
    """_summary_

    Args:
        path (str): _description_
        method (str): _description_
        path_output (str): _description_
        path_algues_network (str): _description_
        seeds (str): _description_
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


def prepare_all_solo_bact(input_path, output_dir):
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
    """_summary_

    Args:
        path_alg (str): _description_
        path_all_bact (str): _description_
    """

    #sub.run([f'mkdir /scratch/clucas/full_microbiota'], shell=True)
    #bacts = sorted(glob(path_all_bact+"*"))
    # for bact in bacts:
    #    sub.run(
    #        [f'cp {bact}/{bact.split(".")[-1]}.sbml /scratch/clucas/full_microbiota/{bact.split("/")[-1]}.sbml'], shell=True)

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
            #sub.run([f'rm -r {output_dir}'], shell=True)


def coev_scopes(path_metabolic_networks, path_all_bact, list_algues, path_sbml_algues, output_dir_scopes, seeds):
    prepare_all_solo_bact(path_metabolic_networks, path_all_bact)
    get_all_full_scope(list_algues, path_all_bact,
                       path_sbml_algues, output_dir_scopes, seeds=seeds)


if __name__ == "__main__":

    test_analyse(path=sys.argv[1], method=sys.argv[2], path_output=sys.argv[3], seeds=sys.argv[4])
