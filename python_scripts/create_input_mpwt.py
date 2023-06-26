from glob import glob
import subprocess as sub
import sys


def egg2gbk(input_dir:str, output_dir:str):
    """Prepare les resultats eggnogmapper pour les passer à pathway tools

    Args:
        input_dir(str): chemin vers les genomes
    """
    directory = sorted(glob(input_dir+"*"))
    algues = [x.split("/")[-1] for x in directory]
    sub.run([f'mkdir {output_dir}'], shell = True)
    for index in range(len(algues)):
        bacts = glob(input_dir+ algues[index]+"/*")

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

            

    
        
if __name__ == "__main__":
    egg2gbk(sys.argv[1])   

# input_dir = /groups/phaeo_bact/mstam/genomes/ 