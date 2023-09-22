import sys
import glob
import subprocess as sub
from holointeract.network_generation.rename_gbk_id import rename_id


def annot_eggnog(input_dir:str, output_dir:str):
    files = glob.glob(input_dir+"*")
    #subprocess.run([f'mkdir {output_dir}'], shell=True)
    for algue in files:
        liste_bact = glob.glob(algue+"/*")
        #subprocess.run([f'mkdir {output_dir}{algue.split("/")[-1]}'], shell=True)
        for bact in liste_bact:
            out = output_dir+algue.split("/")[-1]

            sub.run([f'emapper.py --data_dir /db/eggnog/5.0.2 --cpu 30 -i {bact} --itype genome --genepred prodigal --output_dir {out} -o {bact.split("/")[-1]}'], shell=True)

def egg2gbk(input_dir:str, output_dir:str):
    """Prepare les resultats eggnogmapper pour les passer à pathway tools

    Args:
        input_dir(str): chemin vers les genomes
    """
    directory = sorted(glob.glob(input_dir+"*"))
    algues = [x.split("/")[-1] for x in directory]
    sub.run([f'mkdir {output_dir}'], shell = True)
    for index in range(len(algues)):
        bacts = glob.glob(input_dir+ algues[index]+"/*")

        bact_list = [x.split("/")[-1].split(".")[0] for x in bacts]

        bact_list = list(set(bact_list))

        sub.run([f'mkdir {output_dir}{algues[index]}'], shell = True)

       

        for bact in range(len(bact_list)):
            out = output_dir + algues[index]+"/"+bact_list[bact]
            sub.run([f'mkdir {output_dir}{algues[index]}/{bact_list[bact]}'], shell = True)

            genome = input_dir + algues[index] + "/" + bact_list[bact]+ "/" + bact_list[bact]+ ".fna"
            annotations = input_dir + algues[index] + "/" + bact_list[bact]+ "/"+ bact_list[bact] + ".fna.emapper.annotations"
            prot = input_dir + algues[index] + "/" + bact_list[bact] + "/"+bact_list[bact]+ ".faa"
            gff = input_dir + algues[index] + "/" + bact_list[bact] + "/"+bact_list[bact]+ ".fna.emapper.genepred.gff"


            sub.run([f'emapper2gbk genomes -fn {genome} -fp {prot} --gff {gff} -gt eggnog -o {out}/{bact_list[bact]}.gbk -c 30 -a {annotations}'], shell=True)


def build_network_eggnog(input_dir:str, output_dir:str, singularity_path:str):
    new_name = " "
    liste = ['-', '|', '(', ')', '\'', '=', '#', '*', ':', '!', '+', '[', ']', ',', " "]

    list_algue = sorted(glob.glob(input_dir+"*"))
    sub.run([f'mkdir {output_dir}'], shell=True)

    for algue in  list_algue:
        genome = algue
        sub.run([f'mkdir {output_dir}{algue.split("/")[-1]}'], shell=True)

        cpu="30"
        bacts = sorted(glob.glob(algue+"/*"))
        
        for b in bacts:
            for ch in b.split("/")[-1]:
                if ch in liste:
                    new_name = rename_id(b.split("/")[-1].split(".")[0], b)

            if new_name != " ":
                    sub.run([f'mkdir {b.split(".")[0]}'], shell=True)
                    sub.run([f'mv {new_name} {b.split(".")[0]}'], shell=True)
                    genome = new_name
                    new_name= " "
                    
            else:
                print(b, "is ok")

        sub.run([f'singularity exec -B {singularity_path}:{singularity_path} {singularity_path}/m2m_c.sif m2m recon --clean --noorphan -g {genome} -o {output_dir}{algue.split("/")[-1]} -c {cpu} '], shell = True)
        
def networks_from_genomes(genomes_path,gbk_files,metabolic_networks_path,singularity_path):
    print("Annotation done")
    annot_eggnog(input_dir=genomes_path, output_dir=genomes_path)

    print("Start emapper2GBK")
    egg2gbk(input_dir=genomes_path, output_dir=gbk_files)

    print("Start metabolic network building")
    build_network_eggnog(
        input_dir=gbk_files, output_dir=metabolic_networks_path,singularity_path=singularity_path)

if __name__ == "__main__":
    networks_from_genomes(genomes_path= sys.argv[1],gbk_files=sys.argv[2],metabolic_networks_path=sys.argv[3],singularity_path=sys.argv[4])

