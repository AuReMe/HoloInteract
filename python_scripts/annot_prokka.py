import glob
import subprocess as sub

restant = ["Acidovorax_sp_OV235_new",  "Devosia_sp_I507",  "Nocardioides_alkalitolerans_DSM_16699",  "Rhodobacter_sp_SW2_new",\
"Acidovorax_sp_SD340,Flagellimonas_aquimarina",  "Halomonas_sp_KO116_new",  "Parvibaculum_lavamentivorans_DS-1_new",  "Sinorhizobium_sp_RAC02",\
"Aeromicrobium_sp_Leaf289_new", "Fuerstia_marisgermanicae",  "Maribacter_aquivivus",  "Phenylobacterium_hankyongense_new",  "Sphingopyxis_sp_UBA6723",\
"Caulobacteraceae_bacterium,Gammaproteobacteria_bacterium,Maribacter_caenipelagi",  "Pseudolabrys_sp_Root1462_new",  "Sulfitobacter_mediterraneus",\
"Caulobacteraceae_bacterium_UBA3198","Gemmata_obscuriglobus_UQM_2246",  "Microthrixaceae_bacterium_new", "Psychrobacter_pacificensis"]

liste_algues = sorted(glob.glob("/groups/phaeo_bact/mstam/genomes/*"))
out = "/scratch/clucas/prokka_gff/"


for algue in liste_algues:
    sub.run([f'mkdir {out}{algue.split("/")[-1]}'], shell=True)
    liste_bact = glob.glob(algue+"/*")
    liste_simple = list(set([x.split('.')[0] for x in liste_bact]))
    
    for j in liste_simple:
        if j.split("/")[-1] in restant:
            liste = glob.glob("/scratch/clucas/finish_prokka/"+j.split("/")[-1]+"/*")
            for file in liste:
                if j.split("/")[-1]+"new_nodup.fna" in file:
                    genome = file
                elif j.split("/")[-1]+"new.fna" in file:
                    genome = file
                elif j.split("/")[-1]+".fna" in file:
                    genome = file

                if j.split("/")[-1]+"new.faa" in file:
                    prot = file
                elif j.split("/")[-1]+".faa" in file:
                    prot = file

        else:
            if j+"new_nodup.fna" in liste_bact:
                genome = j+"new_nodup.fna"
            elif j+"new.fna" in liste_bact:
                genome = j+"new.fna"
            elif j+".fna" in liste_bact:
                genome = j+".fna"

            if j+"new.faa" in liste_bact:
                prot = j+"new.faa"
            elif j+".faa" in liste_bact:
                prot = j+".faa" 

        print(genome, prot)
        sub.run([f'prokka --proteins {prot} --outdir {out}{algue.split("/")[-1]}/{j.split("/")[-1]} --cpus 30 --prefix {genome.split("/")[-1].split(".")[0]} {genome}'], shell = True)



    """for j in liste_bact:
        if j.split("/")[-1] in restant:
            new_file = ""
            with open(j+'.fna', 'r') as fi:
                for line in fi:
                    if line[0] == ">" and "NODE" in line:
                        newline = line.split("NODE")[0]+line.split("NODE")[1].split(".")[-1]
                        new_file += newline
                    else:
                        new_file += line
            with open(out+j.split("/")[-1]+'_new.fna',"w") as fo:
                fo.write(new_file)
            new_file = ""
            with open(j+'.faa', 'r') as fi:
                for line in fi:
                    if line[0] == ">" and "NODE" in line:
                        
                       
                        newline = line.split("NODE")[0]+line.split("NODE")[1].split(".")[-1]
                        new_file += newline
                    else:
                        new_file += line
            with open(out+j.split("/")[-1]+'_new.faa',"w") as fo:
                fo.write(new_file)"""

