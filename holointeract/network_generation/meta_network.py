import glob
import subprocess as sub
import sys
from holointeract.rename_gbk_id import rename_id


def build_network_eggnog(input_dir: str, output_dir: str, singularity_path: str):
    """Generate metabolic network

    Args:
        input_dir (str): path to gbk files
        output_dir (str): directory for metabolic networks
        singularity_path (str): path to singularity with pathway tools
    """
    new_name = " "
    liste = ['-', '|', '(', ')', '\'', '=', '#', '*', ':',
             '!', '+', '[', ']', ',', " "]

    list_algue = sorted(glob.glob(input_dir+"*"))
    sub.run([f'mkdir {output_dir}'], shell=True)

    for algue in list_algue:
        genome = algue
        sub.run([f'mkdir {output_dir}{algue.split("/")[-1]}'], shell=True)

        cpu = "30"
        bacts = sorted(glob.glob(algue+"/*"))

        for b in bacts:
            for ch in b.split("/")[-1]:
                if ch in liste:
                    new_name = rename_id(b.split("/")[-1].split(".")[0], b)

            if new_name != " ":
                sub.run([f'mkdir {b.split(".")[0]}'], shell=True)
                sub.run([f'mv {new_name} {b.split(".")[0]}'], shell=True)
                genome = new_name
                new_name = " "

            else:
                print(b, "is ok")

        sub.run(
            [f'singularity exec -B {singularity_path}:{singularity_path} {singularity_path}/m2m_c.sif m2m recon --clean --noorphan -g {genome} -o {output_dir}{algue.split("/")[-1]} -c {cpu} '], shell=True)


def build_network_prokka():
    new_name = " "

    liste = ['-', '|', '/', '(', ')', '\'', '=', '#',
             '*', ':', '!', '+', '[', ']', ',', " "]

    list_algue = sorted(glob.glob(sys.argv[1]+"/*"))
    sub.run([f'mkdir /scratch/clucas/sbml_prokka'], shell=True)

    for algue in list_algue:
        genome = algue
        output = "/scratch/clucas/sbml_prokka/"+algue.split("/")[-1]
        sub.run(
            [f'mkdir /scratch/clucas/sbml_prokka/{algue.split("/")[-1]}'], shell=True)

        cpu = "30"
        bacts = sorted(glob.glob(algue+"/*"))
        if len(bacts) > 1:
            for b in bacts:
                for ch in b.split("/")[-1]:
                    if ch in liste:
                        new_name = rename_id(b.split("/")[-1].split(".")[0], b)

                if new_name != " ":
                    sub.run([f'mkdir {b.split(".")[0]}'], shell=True)
                    sub.run([f'mv {new_name} {b.split(".")[0]}'], shell=True)
                    genome = new_name
                    new_name = " "

                else:
                    sub.run([f'mkdir {b.split(".")[0]}'], shell=True)
                    sub.run([f'mv {b} {b.split(".")[0]}'], shell=True)
            sub.run(
                [f'singularity exec -B /scratch/clucas:/scratch/clucas /scratch/clucas/m2m_c.sif m2m recon --noorphan --clean -g {genome} -o {output} -c {cpu} '], shell=True)


if __name__ == "__main__":
    build_network_eggnog(input_dir=sys.argv[1], output_dir=sys.argv[2])
    # build_network_prokka()
