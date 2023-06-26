#!/bin/bash

# file submission.SBATCH
#SBATCH -c 30
#SBATCH --job-name="holointeract"
#SBATCH --output=workflow.log
#SBATCH --mail-user=corentin.lucas@inria.fr
#SBATCH --mail-type=BEGIN,FAIL,END

. /local/env/envconda.sh
conda activate .env_metabo/
. /local/env/enveggnog-mapper-2.1.9.sh

#python python_scripts/network_from_genomes.py /scratch/clucas/HoloInteract/toy_example/Bact_genomes/ \
                                            /scratch/clucas/HoloInteract/toy_example/gbk_files/ \
                                            /scratch/clucas/HoloInteract/toy_example/Bact_networks/ \
                                            /scratch/clucas/



#holointeract metabolic_analysis \
#                -n /scratch/clucas/HoloInteract/toy_example/Bact_networks/ \
#                -b /scratch/clucas/HoloInteract/toy_example/Bact_scopes/ \
#                -m solo \
#                --host_metabolic_networks_path /scratch/clucas/HoloInteract/toy_example/alg_networks/ \
#                -s /scratch/clucas/HoloInteract/toy_example/seeds_seawater_artefact.sbml \
#                -a /scratch/clucas/HoloInteract/toy_example/alg_scopes/ \
#                -i histogramme_added_value \
#                -c added_value_compound \
#                --clustermap_name clustermap \
#                -t ward \
#                -x matrice

holointeract coevolution_analysis \
                -n /scratch/clucas/HoloInteract/toy_example/Bact_networks/ \
                -b /scratch/clucas/HoloInteract/toy_example/Bact_scopes/ \
                -m solo \
                --host_metabolic_networks_path /scratch/clucas/HoloInteract/toy_example/alg_networks/ \
                -s /scratch/clucas/HoloInteract/toy_example/seeds_seawater_artefact.sbml \
                -a /scratch/clucas/HoloInteract/toy_example/alg_scopes/ \
                -i histogramme_added_value \
                -c added_value_compound \
                -t ward \
                -x matrice \
                --all_bact /scratch/clucas/HoloInteract/toy_example/all_bact/ \
                -f coevolution_graph \
                --all_bact_scopes /scratch/clucas/HoloInteract/toy_example/all_bact_scopes/ \
                --phylogenetic_tree /scratch/clucas/HoloInteract/toy_example/SpeciesTree_rooted.txt \
                --output_values_coev_graph output_values_coevolution


