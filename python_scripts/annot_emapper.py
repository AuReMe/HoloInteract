import subprocess
import glob
import sys


def annot_eggnog(input_dir: str, output_dir: str):
    """Performs eggnogmapper annotation

    Args:
        input_dir (str): directory to host names directories containing microbiota genomes
        output_dir (str): directory to genome annotation
    """
    files = glob.glob(input_dir+"*")
    for algue in files:
        liste_bact = glob.glob(algue+"/*")
        for bact in liste_bact:
            out = output_dir+algue.split("/")[-1]

            subprocess.run(
                [f'emapper.py --data_dir /db/eggnog/5.0.2 --cpu 30 -i {bact} --itype genome --genepred prodigal --output_dir {out} -o {bact.split("/")[-1]}'], shell=True)


if __name__ == "__main__":
    annot_eggnog(input_dir=sys.argv[1], output_dir=sys.argv[2])
