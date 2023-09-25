"""
HoloInteract
"""
# ### IMPORTS
# ======================================================================================================================
import os.path
from ontosunburst.class_metabolites import proportion_workflow

import holointeract.network_generation.annot_emapper
import holointeract.network_generation.create_input_mpwt
import holointeract.network_generation.meta_network

from holointeract.metabolic_analysis.scopes_community import *
from holointeract.metabolic_analysis.heatmap_host_bacteria import *

import holointeract.coevolution_analysis.mat_dist_full_crossed
import holointeract.coevolution_analysis.Graph_dist

from sys import argv
from argparse import ArgumentParser
from rich.traceback import install
from rich import print

# ### PARSING
# ======================================================================================================================

install(show_locals=True)
LINKAGE_METHODS = ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']


def get_command_line_args():
    parser = ArgumentParser(prog='HoloInteract', description='', epilog='')
    subparsers = parser.add_subparsers(help='Available subcommands', dest='subcommands')
    args_metabolic_analysis(subparsers)
    args = parser.parse_args()
    return args


def args_metabolic_analysis(subparsers):
    parser_metabolic_analysis = subparsers.add_parser('metabolic_analysis',
                                                      help='Performs every steps to analysis the metabolic interactions'
                                                           ' in holobionts. Metabolic networks required.')
    parser_metabolic_analysis.add_argument('-cn', '--community_networks', type=str, required=True,
                                           help='path to community networks in SBML')
    parser_metabolic_analysis.add_argument('-hn', '--host_networks', type=str, required=True,
                                           help='path to hosts networks in SBML')
    parser_metabolic_analysis.add_argument('-o', '--output', type=str, required=True,
                                           help='path to output directory')
    parser_metabolic_analysis.add_argument('-s', '--seeds', type=str, required=True,
                                           help='path to seeds SBML file')
    parser_metabolic_analysis.add_argument('-n', '--name', type=str, required=False, default='run',
                                           help='output files name')
    parser_metabolic_analysis.add_argument('-am', '--analysis_method', type=str, required=False,
                                           choices=[SOLO_METHOD, COOP_METHOD], default=COOP_METHOD,
                                           help='method of analysis')
    parser_metabolic_analysis.add_argument('-cm', '--clustering_method', type=str, required=False,
                                           choices=LINKAGE_METHODS, default='ward',
                                           help='method for linkage in clustering')
    parser_metabolic_analysis.add_argument('--max_clust', type=int, required=False, default=10,
                                           help='path to seeds SBML file')
    parser_metabolic_analysis.add_argument('-cpu', type=int, required=False, default=1,
                                           help='number of cpu to use')
    return subparsers


# ### FUNCTIONS
# ======================================================================================================================


def networks_from_genomes(genomes_path, gbk_files, metabolic_networks_path, singularity_path):
    print("Annotation done")
    holointeract.network_generation.annot_emapper.annot_eggnog(
        input_dir=genomes_path, output_dir=genomes_path)

    print("Start emapper2GBK")
    holointeract.network_generation.create_input_mpwt.egg2gbk(
        input_dir=genomes_path, output_dir=gbk_files)

    print("Start metabolic network building")
    holointeract.network_generation.meta_network.build_network_eggnog(
        input_dir=gbk_files, output_dir=metabolic_networks_path, singularity_path=singularity_path)


def metabolic_analysis(community_networks_path, host_networks_path, output_path, seeds, output_name, analysis_method,
                       clustering_method, max_clust, cpu):

    print("Start calculate community scopes")
    community_scopes(community_sbml_path=community_networks_path, hosts_sbml_path=host_networks_path,
                     output_dir=output_path, seeds=seeds, method=analysis_method, cpu=cpu)

    print("Start generating clustermap")
    input_heatmap = os.path.join(output_path, SCOPES_STR, analysis_method)
    output_heatmap = create_heatmap_output(output_path, output_name, analysis_method)
    bact_metabolites, host_metabolites, holo_metabolites, all_metabolites = heatmap_host_bacteria(
        input_dir=input_heatmap, output=output_heatmap, method=clustering_method, max_clust=max_clust)

    proportion_workflow(set(all_metabolites), output=output_heatmap + '_classes_cpd_info')


def create_heatmap_output(output_path, output_name, analysis_method):
    output_heatmap = os.path.join(output_path, 'heatmap')
    create_new_dir(output_heatmap)
    output_heatmap = os.path.join(output_heatmap, analysis_method)
    create_new_dir(output_heatmap)
    output_heatmap = os.path.join(output_heatmap, output_name)
    return output_heatmap


# def coevolution_analysis(metabolic_networks_path, scopes_path, analysis_method, algaes_network_path, seeds, alg_scopes,
#                          all_scopes, all_bact, phylogenetic_tree, output_fig_cpd, output_file_cpd, clustering_method,
#                          csv_file_for_graph, coevolution_graph_name, correction, padmet, matrice_name="matrice"):
#     metabolic_analysis()
#
#     print("Start matrix creation")
#     if analysis_method == "solo":
#         dico_algue, dico_bact = holointeract.metabolic_analysis.create_big_tab_bact.job(
#             path_alg=alg_scopes, path_bact=scopes_path, path_holo=scopes_path,
#             path_metabolic_networks=metabolic_networks_path, output_name=matrice_name)
#     else:
#         dico_algue, dico_bact = holointeract.metabolic_analysis.create_big_tab_alg.job(
#             alg_scopes=alg_scopes, bact_scopes=scopes_path, sbml_path=metabolic_networks_path, output_name=matrice_name,
#             method=analysis_method)
#
#     holointeract.metabolic_analysis.stat_cpd.job(
#         matrice_name + ".csv", output_fig_cpd=output_fig_cpd, output_file_cpd=output_file_cpd, padmet=padmet)
#
#     list_bact = [x for x in dico_bact.keys()]
#     list_algue = [x for x in dico_algue.keys()]
#     print(list_bact)
#     print("Start heatmap")
#     holointeract.metabolic_analysis.heatmap.heatmap(matrice_name + ".csv", clustering_method,
#                                                     output_file=matrice_name, color="tab10")
#
#     holointeract.metabolic_analysis.analyse_community.coev_scopes(metabolic_networks_path, path_all_bact=all_bact,
#                                                                   list_algues=list_algue,
#                                                                   path_sbml_algues=algaes_network_path,
#                                                                   output_dir_scopes=all_scopes, seeds=seeds)
#     # On continue vers la coévolution
#     holointeract.coevolution_analysis.mat_dist_full_crossed.job(list_algue=list_algue, list_bact=list_bact,
#                                                                 scopes_bacteries_path=all_scopes,
#                                                                 output_name=matrice_name + "_coevolution")
#
#     # Construction du graph de coévolution
#
#     holointeract.coevolution_analysis.Graph_dist.job(phylogenetic_tree, input_file=matrice_name + "_coevolution.csv",
#                                                      ouput_file_for_graph=csv_file_for_graph,
#                                                      graph_name=coevolution_graph_name, correction=correction)


def main():
    args = get_command_line_args()

    if len(argv) == 1:
        print('[dark_orange]You need to provide a command and its arguments for the program to work.\n'
              'Try to use -h or --help to get list of available commands.')
        exit()

    elif args.subcommands == 'metabolic_analysis':
        metabolic_analysis(community_networks_path=args.community_networks, host_networks_path=args.host_networks,
                           output_path=args.output, seeds=args.seeds, output_name=args.name,
                           analysis_method=args.analysis_method, clustering_method=args.clustering_method,
                           max_clust=args.max_clust, cpu=args.cpu)

    else:
        print('[dark_orange]Unknown command. Please use the help (-h) to see available commands.')
        exit(1)


# holointeract metabolic_analysis -cn example/inputs/community/ -hn example/inputs/hosts/ -o example/outputs/ -s example/inputs/seeds/seeds_seawater_artefact.sbml
