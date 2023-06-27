import python_scripts.annot_emapper
import python_scripts.create_input_mpwt
import python_scripts.meta_network
import python_scripts.analyse_community
import python_scripts.get_algae_scope
import python_scripts.create_big_tab_alg
import python_scripts.stat_cpd
import python_scripts.create_big_tab_bact
import python_scripts.heatmap
import python_scripts.mat_dist_full_crossed
import python_scripts.Graph_dist
from sys import argv
from argparse import ArgumentParser
from rich.traceback import install
from rich import print
install(show_locals=True)


parser: ArgumentParser = ArgumentParser(
    description='HoloInteract', add_help=True)
subparsers = parser.add_subparsers(
    help='Available subcommands', dest="subcommands")

parser._positionals.title = 'Subcommands'
parser._optionals.title = 'Global Arguments'

# Annotation
parser_annot: ArgumentParser = subparsers.add_parser(
    'annot_emapper',
    help="Performs emapper annotations\n"
    "The bacteria genomes must be in a type repository structure like : ../Genomes/Host/Bacteria/Bacteria_genome.fna . You have to give the path to ../Genomes/"
)
parser_annot.add_argument("--repository", "-r", type=str, required=True,
                          help="Path to the repository ../Genomes/ containing hosts folders.")
parser_annot.add_argument("--out", "-o", type=str,
                          required=True, help="Ouput directory path")

# GBK creation

parser_gbk: ArgumentParser = subparsers.add_parser(
    'create_gbk',
    help="Performs emapper2gbk on emapper annotations\n"
    "Creates Gbk files."
)
parser_gbk.add_argument("--repository", "-r", type=str, required=True,
                        help="Path to the repository containing eggnogmapper annotations.")
parser_gbk.add_argument("--out", "-o", type=str,
                        required=True, help="Ouput directory path")

# Reseau metabolic

parser_meta_network: ArgumentParser = subparsers.add_parser(
    'meta_network',
    help="Creates metabolic Networks of the generated annotations\n"
    "Creates SBML files."
)
parser_meta_network.add_argument("--repository", "-r", type=str,
                                 required=True, help="Path to the repository containing gbk files.")
parser_meta_network.add_argument(
    "--out", "-o", type=str, required=True, help="Ouput directory path")
parser_meta_network.add_argument(
    "--singularity_path", "-s", type=str, required=True, help="path to singularity image")


# Alg_scope

parser_alg_scope: ArgumentParser = subparsers.add_parser(
    'get_alg_scope',
    help="Computes the scopes of the metabolic network of all the algaes, using Metage2Metabo.\n"
)
parser_alg_scope.add_argument("--repository", "-r", type=str, required=True,
                              help="Path to the repository containing algaes metabolic networks.")
parser_alg_scope.add_argument("--seeds", "-s", type=str, required=True,
                              help="Seeds file to initiate metabolic network computing.")
parser_alg_scope.add_argument(
    "--out", "-o", type=str, required=True, help="Ouput directory path")

# Analyse community
parser_community_analysis: ArgumentParser = subparsers.add_parser(
    'community_analysis',
    help="Computes the scopes of the metabolic network of all the algaes, using Metage2Metabo.\n"
)
parser_community_analysis.add_argument("--repository", "-r", type=str, required=True,
                                       help="Path to the metabolic network repository at the level of algaes.")
parser_community_analysis.add_argument(
    "--seeds", "-s", type=str, required=True, help="Seeds file to initiate metabolic network computing.")
parser_community_analysis.add_argument("--method", type=str, required=True, choices=[
                                       "solo", "coop"], help="Authorize or not bacteria cooperation inside a microbiota")
parser_community_analysis.add_argument("--alg_networks_dir", "-a", type=str,
                                       required=True, help="Path to the repository containing algaes metabolic networks.")
parser_community_analysis.add_argument(
    "--out", "-o", type=str, required=True, help="Ouput directory path")

# Stats on added value

parser_stat_cpd: ArgumentParser = subparsers.add_parser(
    'stat_cpd',
    help="Compute some statistics on metabolites produced by addedvalue in holobionts\n"
)
parser_stat_cpd.add_argument("--matrice", "-i", type=str, required=True,
                             help="Path to the matrice produced by the function create_big_tab_alg or create_big_tab_bact.")
parser_stat_cpd.add_argument(
    "--histogram_name", "-n", type=str, required=True, help="Name of the graph")
parser_stat_cpd.add_argument(
    "--padmet_file", "-p", type=str, required=True, help="Padmet database file")
parser_stat_cpd.add_argument(
    "--output_file_name", "-o", type=str, required=True, help="Name of the csv file")

# Clustermap

parser_clustermap: ArgumentParser = subparsers.add_parser(
    'clustermap',
    help="Generates a cluster map of the scopes matrice\n"
)
parser_clustermap.add_argument("--matrice", "-i", type=str, required=True,
                               help="Path to the matrice produced by the function create_big_tab_alg or create_big_tab_bact.")
parser_clustermap.add_argument(
    "--clustermap_name", "-n", type=str, required=True, help="Name of the graph")
parser_clustermap.add_argument("--method", type=str, required=True, choices=[
                               'single', 'complete', 'average', 'weighted', 'centroid', 'ward'], help="Clustering method")

# Coevolution scopes

parser_coevolution_scopes: ArgumentParser = subparsers.add_parser(
    'coevolution_scopes',
    help="Computes the scopes of the metabolic network of all the bacteria with all algaes, using Metage2Metabo.\n"
)
parser_coevolution_scopes.add_argument("--repository", "-r", type=str, required=True,
                                       help="Path to the metabolic network repository at the level of algaes.")
parser_coevolution_scopes.add_argument(
    "--seeds", "-s", type=str, required=True, help="Seeds file to initiate metabolic network computing.")
parser_coevolution_scopes.add_argument(
    "--list_alg", "-l", type=str, required=True, help="A list of all algaes name.")
parser_coevolution_scopes.add_argument("--all_bact", "-p", type=str, required=True,
                                       help="Path to create the pool of all bacterias of the dataset (in the singularity path)")
parser_coevolution_scopes.add_argument("--alg_networks_dir", "-a", type=str,
                                       required=True, help="Path to the repository containing algaes metabolic networks.")
parser_coevolution_scopes.add_argument("--out", "-o", type=str, required=True,
                                       help="Ouput directory path to stock the scopes of all bacteria with all algaes")

# Coevolution matrix

parser_mat_full_crossed: ArgumentParser = subparsers.add_parser(
    'coevolution_matrix',
    help="Generates the matrix using the scopes of all bacterias with all algaes\n"
)
parser_mat_full_crossed.add_argument(
    "--list_alg", "-l", type=str, required=True, help="A list of all algaes name.")
parser_mat_full_crossed.add_argument(
    "--list_bact", "-b", type=str, required=True, help="A list of all bact name.")
parser_mat_full_crossed.add_argument("--all_scopes", "-r", type=str, required=True,
                                     help="path to the directory containing the scopes of all bacteria with all algaes")
parser_mat_full_crossed.add_argument(
    "--out", "-o", type=str, required=True, help="Ouput name of the matrix")

# Coevolution graph

parser_coevolution_graph: ArgumentParser = subparsers.add_parser(
    'coevolution_graph',
    help="Genrates a graph of the metabolic complementarity as a function of phylogenetic distance to host.\n"
)
parser_coevolution_graph.add_argument("--phylogenetic_tree", "-t", type=str,
                                      required=True, help="A phylogenetic tree (in Newick format) of the hosts.")
parser_coevolution_graph.add_argument(
    "--coevolution_matrix", "-m", type=str, required=True, help="Coevolution matrix")
parser_coevolution_graph.add_argument(
    "--csv_file_name", "-r", type=str, required=True, help="Name of csv file")
parser_coevolution_graph.add_argument(
    "--graph_name", "-o", type=str, required=True, help="Name of the graph")
parser_coevolution_graph.add_argument("--correction", "-c", type=str, required=True, choices=[
                                      "bonferroni", "benjamini", ""], help="Multiple tests correction")


# Metabolic network reconstruction
# parser_networks_from_genomes: ArgumentParser = subparsers.add_parser(
#     'networks_from_genomes',
#     help="Reconstruction of the metabolic networks from the bacteria genomes\n"
# )
# parser_full_metabolic_analysis.add_argument("-g","--genomes_path", type=str, required=True,help="The bacteria genomes must be in a type repository structure like : ../Genomes/Host/Bacteria/Bacteria_genome.fna . You have to give the path to ../Genomes/")
# parser_full_metabolic_analysis.add_argument("-n", "--metabolic_networks_path",type=str, required=True,help="path to metabolic networks")
# parser_full_metabolic_analysis.add_argument("-p","--singularity_path", type=str, required=True,help="path to singularity")
# parser_full_metabolic_analysis.add_argument("-k", "--gbk_files", type=str, required=True,help="path to gbk files")


# METABOLIC ANALYSIS ONLY

# parser_full_metabolic_analysis: ArgumentParser = subparsers.add_parser(
#     'full_metabolic_analysis',
#     help="Perform every steps to analysis the metabolic interaction in holobionts.\n"
# )
# parser_full_metabolic_analysis.add_argument("-g","--genomes_path", type=str, required=True,help="The bacteria genomes must be in a type repository structure like : ../Genomes/Host/Bacteria/Bacteria_genome.fna . You have to give the path to ../Genomes/")
# parser_full_metabolic_analysis.add_argument("-n", "--metabolic_networks_path",type=str, required=True,help="path to metabolic networks")
# parser_full_metabolic_analysis.add_argument("-b", "--bact_scopes_path", type=str, required=True,help="path to scopes")
# parser_full_metabolic_analysis.add_argument("-m","--analysis_method",choices=["solo", "coop"], type=str, required=True,help="Authorize cooperation between bacteria or not")
# parser_full_metabolic_analysis.add_argument("-p","--singularity_path", type=str, required=True,help="path to singularity")
# parser_full_metabolic_analysis.add_argument("--host_metabolic_networks_path", type=str, required=True,help="host_metabolic network")
# parser_full_metabolic_analysis.add_argument("-s", "--seeds", type=str, required=True,help="seeds file")
# parser_full_metabolic_analysis.add_argument("-a","--host_scopes_path", type=str, required=True,help="path to hosr scopes")
# parser_full_metabolic_analysis.add_argument("-k", "--gbk_files", type=str, required=True,help="path to gbk files")
# parser_full_metabolic_analysis.add_argument("-i","--histogram_cpd_name", type=str, required=True,help="histogram added value compounds name")
# parser_full_metabolic_analysis.add_argument("-c","--csv_cpd_name", type=str, required=True,help="name of the csv file containing added value compounds")
# parser_full_metabolic_analysis.add_argument("--clustermap_name" ,type=str, required=True,help="name of the clustermap")
# parser_full_metabolic_analysis.add_argument("-t","--clustering_method", type=str, required=True, choices=['single', 'complete', 'average','weighted', 'centroid','ward'],help="clustering method")
# parser_full_metabolic_analysis.add_argument("-x", "--matrix_name",type=str, required=True,help="name of the csv_file permitting to create the clustermap")

parser_metabolic_analysis: ArgumentParser = subparsers.add_parser(
    'metabolic_analysis',
    help="Perform every steps to analysis the metabolic interaction in holobionts if you already have the metabolic networks.\n"
)
parser_metabolic_analysis.add_argument(
    "-n", "--metabolic_networks_path", type=str, required=True, help="path to metabolic networks")
parser_metabolic_analysis.add_argument(
    "-b", "--bact_scopes_path", type=str, required=True, help="path to scopes")
parser_metabolic_analysis.add_argument("-m", "--analysis_method", choices=[
                                       "solo", "coop"], type=str, required=True, help="Authorize cooperation between bacteria or not")
parser_metabolic_analysis.add_argument(
    "--host_metabolic_networks_path", type=str, required=True, help="host_metabolic network")
parser_metabolic_analysis.add_argument(
    "-s", "--seeds", type=str, required=True, help="seeds file")
parser_metabolic_analysis.add_argument(
    "-a", "--host_scopes_path", type=str, required=True, help="path to hosr scopes")
parser_metabolic_analysis.add_argument(
    "-i", "--histogram_cpd_name", type=str, required=True, help="histogram added value compounds name")
parser_metabolic_analysis.add_argument(
    "-c", "--csv_cpd_name", type=str, required=True, help="name of the csv file containing added value compounds")
parser_metabolic_analysis.add_argument(
    "--padmet_file", type=str, required=True, help="padmet database file")
parser_metabolic_analysis.add_argument(
    "--clustermap_name", type=str, required=True, help="name of the clustermap")
parser_metabolic_analysis.add_argument("-t", "--clustering_method", type=str, required=True, choices=[
                                       'single', 'complete', 'average', 'weighted', 'centroid', 'ward'], help="clustering method")
parser_metabolic_analysis.add_argument(
    "-x", "--matrix_name", type=str, required=True, help="name of the csv_file permitting to create the clustermap")

# FULL WORKFLOW

# parser_full_coevolution_analysis: ArgumentParser = subparsers.add_parser(
#     'full_coevolution_analysis',
#     help="Perform every steps to analysis the metabolic interaction in holobionts.\n"
# )
# parser_full_coevolution_analysis.add_argument("-g","--genomes_path", type=str, required=True,help="The bacteria genomes must be in a type repository structure like : ../Genomes/Host/Bacteria/Bacteria_genome.fna . You have to give the path to ../Genomes/")
# parser_full_coevolution_analysis.add_argument("-n", "--metabolic_networks_path",type=str, required=True,help="path to metabolic networks")
# parser_full_coevolution_analysis.add_argument("-b", "--bact_scopes_path", type=str, required=True,help="path to scopes")
# parser_full_coevolution_analysis.add_argument("-m","--analysis_method",choices=["solo", "coop"], type=str, required=True,help="Authorize cooperation between bacteria or not")
# parser_full_coevolution_analysis.add_argument("--singularity_path", type=str, required=True,help="path to singularity")
# parser_full_coevolution_analysis.add_argument("--host_metabolic_networks_path", type=str, required=True,help="host_metabolic network")
# parser_full_coevolution_analysis.add_argument("-s", "--seeds", type=str, required=True,help="seeds file")
# parser_full_coevolution_analysis.add_argument("-a","--host_scopes_path", type=str, required=True,help="path to hosr scopes")
# parser_full_coevolution_analysis.add_argument("-k", "--gbk_files", type=str, required=True,help="path to gbk files")
# parser_full_coevolution_analysis.add_argument("-i","--histogram_cpd_name", type=str, required=True,help="histogram added value compounds name")
# parser_full_coevolution_analysis.add_argument("-c","--csv_cpd_name", type=str, required=True,help="name of the csv file containing added value compounds")
# parser_full_coevolution_analysis.add_argument("--clustermap_name" ,type=str, required=True,help="name of the clustermap")
# parser_full_coevolution_analysis.add_argument("-t","--clustering_method", type=str, required=True, choices=['single', 'complete', 'average','weighted', 'centroid','ward'],help="clustering method")
# parser_full_coevolution_analysis.add_argument("-x", "--matrix_name",type=str, required=True,help="name of the csv_file permitting to create the clustermap")
# parser_full_coevolution_analysis.add_argument("--all_bact","-p", type=str, required=True,help="Path to create the pool of all bacterias of the dataset (in the singularity path)")
# parser_full_coevolution_analysis.add_argument("--coevolution_matrix",type=str, required=True,help="name of the csv_file permitting to create the clustermap")
# parser_full_coevolution_analysis.add_argument("-f", "--coevolution_graph",type=str, required=True,help="name of the csv_file permitting to create the clustermap")
# parser_full_coevolution_analysis.add_argument("--all_bact_scopes", "-r", type=str, required=True, help="path to the directory containing the scopes of all bacteria with all algaes")
# parser_full_coevolution_analysis.add_argument("--coevolution_matrix_name", "-e", type=str, required=True, help="name of the coevolution matrix")
# parser_full_coevolution_analysis.add_argument("--phylogenetic_tree", "-y", type=str, required=True, help="phylogenetic_tree (Newick format) of the hosts")
# parser_full_coevolution_analysis.add_argument("--output_values_coev_graph", "-z", type=str, required=True, help="phylogenetic_tree (Newick format) of the hosts")

parser_coevolution_analysis: ArgumentParser = subparsers.add_parser(
    'coevolution_analysis',
    help="Perform every steps to analysis the metabolic interaction in holobionts.\n"
)
parser_coevolution_analysis.add_argument(
    "-n", "--metabolic_networks_path", type=str, required=True, help="path to metabolic networks")
parser_coevolution_analysis.add_argument(
    "-b", "--bact_scopes_path", type=str, required=True, help="path to scopes")
parser_coevolution_analysis.add_argument("-m", "--analysis_method", choices=[
                                         "solo", "coop"], type=str, required=True, help="Authorize cooperation between bacteria or not")
parser_coevolution_analysis.add_argument(
    "--host_metabolic_networks_path", type=str, required=True, help="host_metabolic network")
parser_coevolution_analysis.add_argument(
    "-s", "--seeds", type=str, required=True, help="seeds file")
parser_coevolution_analysis.add_argument(
    "-a", "--host_scopes_path", type=str, required=True, help="path to hosr scopes")
parser_coevolution_analysis.add_argument(
    "-i", "--histogram_cpd_name", type=str, required=True, help="histogram added value compounds name")
parser_coevolution_analysis.add_argument(
    "-c", "--csv_cpd_name", type=str, required=True, help="name of the csv file containing added value compounds")
parser_coevolution_analysis.add_argument("-t", "--clustering_method", type=str, required=True, choices=[
                                         'single', 'complete', 'average', 'weighted', 'centroid', 'ward'], help="clustering method")
parser_coevolution_analysis.add_argument(
    "-x", "--matrix_name", type=str, required=True, help="name of the csv_file permitting to create the clustermap")
parser_coevolution_analysis.add_argument("--all_bact", "-p", type=str, required=True,
                                         help="Path to create the pool of all bacterias of the dataset (in the singularity path)")
parser_coevolution_analysis.add_argument("-f", "--coevolution_graph", type=str,
                                         required=True, help="name of the csv_file permitting to create the clustermap")
parser_coevolution_analysis.add_argument("--all_bact_scopes", "-r", type=str, required=True,
                                         help="path to the directory containing the scopes of all bacteria with all algaes")
parser_coevolution_analysis.add_argument(
    "--phylogenetic_tree", "-y", type=str, required=True, help="phylogenetic_tree (Newick format) of the hosts")
parser_coevolution_analysis.add_argument(
    "--output_values_coev_graph", "-z", type=str, required=True, help="phylogenetic_tree (Newick format) of the hosts")
parser_coevolution_analysis.add_argument(
    "--correction", type=str, help="Mutiple test correction", choices=["", "bonferroni", "benjamini"], default="")
parser_metabolic_analysis.add_argument(
    "--padmet_file", type=str, required=True, help="padmet database file")


#############################################################################################

args = parser.parse_args()


def networks_from_genomes(genomes_path, gbk_files, metabolic_networks_path, singularity_path):
    print("Annotation done")
    python_scripts.annot_emapper.annot_eggnog(
        input_dir=genomes_path, output_dir=genomes_path)

    print("Start emapper2GBK")
    python_scripts.create_input_mpwt.egg2gbk(
        input_dir=genomes_path, output_dir=gbk_files)

    print("Start metabolic network building")
    python_scripts.meta_network.build_network_eggnog(
        input_dir=gbk_files, output_dir=metabolic_networks_path, singularity_path=singularity_path)

# def full_metabolic_analysis(genomes_path, metabolic_networks_path, scopes_path, analysis_method, singularity_path, algaes_network_path, seeds, alg_scopes, gbk_files, output_fig_cpd, output_file_cpd, clustering_method,matrice_name="matrice"):
#     # print("Annotation done")
#     # python_scripts.annot_emapper.annot_eggnog(input_dir=genomes_path, output_dir=genomes_path)

#     # print("Start emapper2GBK")
#     # python_scripts.create_input_mpwt.egg2gbk(input_dir=genomes_path, output_dir=gbk_files)

#     # print("Start metabolic network building")
#     # python_scripts.meta_network.build_network_eggnog(
#     #     input_dir=gbk_files, output_dir=metabolic_networks_path,singularity_path=singularity_path)

#     print("Start analysing networks")
#     # Algues
#     python_scripts.get_algae_scope.run_miscoto(
#         path=algaes_network_path, seeds=seeds, output=alg_scopes)
#     # Bacts
#     python_scripts.analyse_community.test_analyse(path=metabolic_networks_path, method=analysis_method,
#                                    path_output=scopes_path, path_algues_network=algaes_network_path, seeds=seeds)

#     print("Start matrix creation")
#     if analysis_method == "solo":
#         dico_algue, dico_bact = python_scripts.create_big_tab_bact.job(
#             path_alg=alg_scopes, path_bact=scopes_path, path_holo=scopes_path, path_metabolic_networks=metabolic_networks_path, output_name=matrice_name)
#     else:
#         dico_algue, dico_bact = python_scripts.create_big_tab_alg.job(
#             alg_scopes=alg_scopes, bact_scopes=scopes_path, sbml_path=metabolic_networks_path, output_name=matrice_name, method=analysis_method)

#     python_scripts.stat_cpd.job(matrice_name+".csv", output_fig_cpd=output_fig_cpd, output_file_cpd=output_file_cpd)

#     print("Start heatmap")
#     python_scripts.heatmap.heatmap(matrice_name+".csv", clustering_method,
#                     analysis_method=analysis_method, output_file=matrice_name, color="tab10")


def metabolic_analysis(metabolic_networks_path, scopes_path, analysis_method, algaes_network_path, seeds, alg_scopes, output_fig_cpd, output_file_cpd, clustering_method, padmet, matrice_name="matrice"):

    print("Start analysing networks")
    # Algues
    python_scripts.get_algae_scope.run_miscoto(
        path=algaes_network_path, seeds=seeds, output=alg_scopes)
    # Bacts
    python_scripts.analyse_community.test_analyse(path=metabolic_networks_path, method=analysis_method,
                                                  path_output=scopes_path, path_algues_network=algaes_network_path, seeds=seeds)

    print("Start matrix creation")
    if analysis_method == "solo":
        dico_algue, dico_bact = python_scripts.create_big_tab_bact.job(
            path_alg=alg_scopes, path_bact=scopes_path, path_holo=scopes_path, path_metabolic_networks=metabolic_networks_path, output_name=matrice_name)
    else:
        dico_algue, dico_bact = python_scripts.create_big_tab_alg.job(
            alg_scopes=alg_scopes, bact_scopes=scopes_path, sbml_path=metabolic_networks_path, output_name=matrice_name, method=analysis_method)

    python_scripts.stat_cpd.job(
        matrice_name+".csv", output_fig_cpd=output_fig_cpd, output_file_cpd=output_file_cpd, padmet=padmet)

    print("Start heatmap")
    python_scripts.heatmap.heatmap(matrice_name+".csv", clustering_method,
                                   output_file=matrice_name, color="tab10")


def coevolution_analysis(metabolic_networks_path, scopes_path, analysis_method, algaes_network_path, seeds, alg_scopes, all_scopes, all_bact, phylogenetic_tree, output_fig_cpd, output_file_cpd, clustering_method, csv_file_for_graph, coevolution_graph_name, correction, padmet, matrice_name="matrice"):

    print("Start analysing networks")
    # Algues
    python_scripts.get_algae_scope.run_miscoto(
        path=algaes_network_path, seeds=seeds, output=alg_scopes)
    # Bacts
    python_scripts.analyse_community.test_analyse(path=metabolic_networks_path, method=analysis_method,
                                                  path_output=scopes_path, path_algues_network=algaes_network_path, seeds=seeds)

    print("Start matrix creation")
    if analysis_method == "solo":
        dico_algue, dico_bact = python_scripts.create_big_tab_bact.job(
            path_alg=alg_scopes, path_bact=scopes_path, path_holo=scopes_path, path_metabolic_networks=metabolic_networks_path, output_name=matrice_name)
    else:
        dico_algue, dico_bact = python_scripts.create_big_tab_alg.job(
            alg_scopes=alg_scopes, bact_scopes=scopes_path, sbml_path=metabolic_networks_path, output_name=matrice_name, method=analysis_method)

    python_scripts.stat_cpd.job(
        matrice_name+".csv", output_fig_cpd=output_fig_cpd, output_file_cpd=output_file_cpd, padmet=padmet)

    list_bact = [x for x in dico_bact.keys()]
    list_algue = [x for x in dico_algue.keys()]
    print(list_bact)
    print("Start heatmap")
    python_scripts.heatmap.heatmap(matrice_name+".csv", clustering_method,
                                   output_file=matrice_name, color="tab10")

    python_scripts.analyse_community.coev_scopes(metabolic_networks_path, path_all_bact=all_bact,
                                                 list_algues=list_algue, path_sbml_algues=algaes_network_path, output_dir_scopes=all_scopes, seeds=seeds)
    # On continue vers la coévolution
    python_scripts.mat_dist_full_crossed.job(list_algue=list_algue, list_bact=list_bact,
                                             scopes_bacteries_path=all_scopes, output_name=matrice_name+"_coevolution")

    # Construction du graph de coévolution

    python_scripts.Graph_dist.job(phylogenetic_tree, input_file=matrice_name+"_coevolution.csv",
                                  ouput_file_for_graph=csv_file_for_graph, graph_name=coevolution_graph_name, correction=correction)


# def full_coevolution_analysis(genomes_path, metabolic_networks_path, scopes_path, analysis_method, singularity_path,all_bact,all_scopes, algaes_network_path, seeds, alg_scopes, gbk_files, phylogenetic_tree, output_fig_cpd, output_file_cpd,clustering_method, out_file_graph,matrice_name="matrice"):
#     print("Annotation done")
#     python_scripts.annot_emapper.annot_eggnog(input_dir=genomes_path, output_dir=genomes_path)

#     print("Start emapper2GBK")
#     python_scripts.create_input_mpwt.egg2gbk(input_dir=genomes_path, output_dir=gbk_files)

#     print("Start metabolic network building")
#     python_scripts.meta_network.build_network_eggnog(
#         input_dir=gbk_files, output_dir=metabolic_networks_path, singularity_path=singularity_path)

#     print("Start analysing networks")
#     # Algues
#     python_scripts.get_algae_scope.run_miscoto(
#         path=algaes_network_path, seeds=seeds, output=alg_scopes)
#     # Bacts
#     python_scripts.analyse_community.test_analyse(path=metabolic_networks_path, method=analysis_method,
#                                    path_output=scopes_path, path_algues_network=algaes_network_path, seeds=seeds)

#     print("Start matrix creation")
#     if analysis_method == "solo":
#         dico_algue, dico_bact = python_scripts.create_big_tab_bact.job(
#             path_alg=alg_scopes, path_bact=scopes_path, path_holo=scopes_path, path_metabolic_networks=metabolic_networks_path, output_name=matrice_name)
#     else:
#         dico_algue, dico_bact = python_scripts.create_big_tab_alg.job(
#             alg_scopes=alg_scopes, bact_scopes=scopes_path, sbml_path=metabolic_networks_path, output_name=matrice_name, method=analysis_method)

#     python_scripts.stat_cpd.job(matrice_name+".csv", output_fig_cpd=output_fig_cpd, output_file_cpd=output_file_cpd)

#     list_bact = [x for x in dico_bact.keys()]
#     list_algue = [x for x in dico_algue.keys()]

#     print("Start heatmap")
#     python_scripts.heatmap.heatmap(matrice_name+".csv", clustering_method,
#                     analysis_method=analysis_method, output_file=matrice_name, color="tab10")

#     python_scripts.analyse_community.coev_scopes(metabolic_networks_path,path_all_bact=all_bact,
#                                   list_algues=list_algue,path_sbml_algues=algaes_network_path,output_dir_scopes=all_scopes, seeds=seeds)
#     # On continue vers la coévolution
#     python_scripts.mat_dist_full_crossed.job(list_algue=list_algue, list_bact=list_bact,
#                               scopes_bacteries_path=all_scopes, output_name=matrice_name+"_coevolution")

#     # Construction du graph de coévolution

#     python_scripts.Graph_dist.job(phylogenetic_tree, input_file=matrice_name+"_coevolution.csv",
#                    ouput_file_for_graph=out_file_graph, graph_name="Coevolution_distance")


def main():
    "Call for subprograms"
    if len(argv) == 1:
        print(
            "[dark_orange]You need to provide a command and its arguments for the program to work.\n"
            "Try to use -h or --help to get list of available commands.")
        exit()

    elif args.subcommands == 'annot_emapper':
        python_scripts.annot_emapper.annot_eggnog(
            input_dir=args.repository, output_dir=args.out)
    elif args.subcommands == "create_gbk":
        python_scripts.create_input_mpwt.egg2gbk(
            input_dir=args.repository, output_dir=args.out)
    elif args.subcommands == "meta_network":
        python_scripts.meta_network.build_network_eggnog(
            input_dir=args.repository, output_dir=args.out, singularity_path=args.singularity_path)
    elif args.subcommands == "get_alg_scope":
        python_scripts.get_algae_scope.run_miscoto(
            path=args.repository, seeds=args.seeds, output=args.out)
    elif args.subcommands == "community_analysis":
        python_scripts.analyse_community.test_analyse(
            path=args.repository, seeds=args.seeds, method=args.method, path_algues_network=args.alg_networks_dir, path_output=args.out)
    elif args.subcommands == "stat_cpd":
        python_scripts.stat_cpd.job(
            matrice=args.matrice, output_fig_cpd=args.histogram_name, output_file_cpd=args.output_file_name, padmet=args.padmet_file)
    elif args.subcommands == "clustermap":
        python_scripts.heatmap.heatmap(
            input_file=args.matrice, method=args.method, output_file=args.clustermap_name)
    elif args.subcommands == "coevolution_scopes":
        python_scripts.analyse_community.coev_scopes(path_metabolic_networks=args.repository, path_all_bact=args.all_bact,
                                                     list_algues=args.list_alg, path_sbml_algues=args.alg_networks_dir, output_dir_scopes=args.out, seeds=args.seeds)
    elif args.subcommands == 'coevolution_matrix':
        python_scripts.mat_dist_full_crossed.job(
            list_algue=args.list_alg, list_bact=args.list_bact, scopes_bacteries_path=args.all_scopes, output_name=args.out)
    elif args.subcommands == 'coevolution_graph':
        python_scripts.Graph_dist.job(phylogenetic_tree=args.phylogenetic_tree, input_file=args.coevolution_matrix,
                                      ouput_file_for_graph=args.csv_file_name, correction=args.correction, graph_name=args.graph_name)

    elif args.subcommands == 'metabolic_analysis':
        metabolic_analysis(metabolic_networks_path=args.metabolic_networks_path, scopes_path=args.bact_scopes_path, analysis_method=args.analysis_method,
                           algaes_network_path=args.host_metabolic_networks_path, seeds=args.seeds,
                           alg_scopes=args.host_scopes_path, output_fig_cpd=args.histogram_cpd_name, output_file_cpd=args.csv_cpd_name,
                           clustering_method=args.clustering_method, matrice_name=args.matrix_name, padmet=args.padmet_file
                           )

    elif args.subcommands == 'coevolution_analysis':
        coevolution_analysis(metabolic_networks_path=args.metabolic_networks_path, scopes_path=args.bact_scopes_path, analysis_method=args.analysis_method,
                             algaes_network_path=args.host_metabolic_networks_path, seeds=args.seeds, alg_scopes=args.host_scopes_path, output_fig_cpd=args.histogram_cpd_name, output_file_cpd=args.csv_cpd_name,
                             clustering_method=args.clustering_method, matrice_name=args.matrix_name, phylogenetic_tree=args.phylogenetic_tree,
                             csv_file_for_graph=args.output_values_coev_graph, correction=args.correction, all_bact=args.all_bact,
                             coevolution_graph_name=args.coevolution_graph, all_scopes=args.all_bact_scopes, padmet=args.padmet_file
                             )
    else:
        print(
            "[dark_orange]Unknown command. Please use the help (-h) to see available commands.")
        exit(1)
