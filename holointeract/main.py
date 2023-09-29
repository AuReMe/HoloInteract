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
from holointeract.metabolic_analysis.heatmap_host_microorganism import *

from holointeract.coevolution_analysis.coevolution_analysis import *

from sys import argv
from argparse import ArgumentParser
from rich.traceback import install
from rich import print

# ======================================================================================================================
install(show_locals=True)
LINKAGE_METHODS = ['single', 'complete', 'average', 'weighted', 'centroid', 'median', 'ward']


# ### PARSING
# ======================================================================================================================

def get_command_line_args():
    parser = ArgumentParser(prog='HoloInteract', description='', epilog='')
    subparsers = parser.add_subparsers(help='Available subcommands', dest='subcommands')
    args_metabolic_analysis(subparsers)
    args_coevolution_analysis(subparsers)
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
                                           choices=[SOLO_METHOD, COOP_METHOD, FULL_METHOD], default=COOP_METHOD,
                                           help='method of analysis')
    parser_metabolic_analysis.add_argument('-cm', '--clustering_method', type=str, required=False,
                                           choices=LINKAGE_METHODS, default='ward',
                                           help='method for linkage in clustering')
    parser_metabolic_analysis.add_argument('--max_clust', type=int, required=False, default=10,
                                           help='path to seeds SBML file')
    parser_metabolic_analysis.add_argument('-cpu', type=int, required=False, default=1,
                                           help='number of cpu to use')


def args_coevolution_analysis(subparsers):
    parser_coevolution = subparsers.add_parser('coevolution',
                                               help='')
    parser_coevolution.add_argument('-cn', '--community_networks', type=str, required=True,
                                    help='path to community networks in SBML')
    parser_coevolution.add_argument('-hn', '--host_networks', type=str, required=True,
                                    help='path to hosts networks in SBML')
    parser_coevolution.add_argument('-o', '--output', type=str, required=True,
                                    help='path to output directory')
    parser_coevolution.add_argument('-s', '--seeds', type=str, required=True,
                                    help='path to seeds SBML file')
    parser_coevolution.add_argument('-n', '--name', type=str, required=False, default='run',
                                    help='output files name')
    parser_coevolution.add_argument('-p', '--phylo_tree', type=str, required=False, default=None,
                                    help='path to phylogenetic tree (Newick format)')
    parser_coevolution.add_argument('-cor', '--correction', type=str, required=False, default=None,
                                    choices=[None, BONFERRONI, BENJAMINI],
                                    help='correction to apply to p-values')
    parser_coevolution.add_argument('-cm', '--clustering_method', type=str, required=False,
                                    choices=LINKAGE_METHODS, default='ward',
                                    help='method for linkage in clustering')
    parser_coevolution.add_argument('--max_clust', type=int, required=False, default=10,
                                    help='path to seeds SBML file')
    parser_coevolution.add_argument('-cpu', type=int, required=False, default=1,
                                    help='number of cpu to use')


# ### FUNCTIONS
# ======================================================================================================================


def networks_from_genomes(genomes_path, gbk_files, metabolic_networks_path, singularity_path):
    pass
    # print("Annotation done")
    # holointeract.network_generation.annot_emapper.annot_eggnog(
    #     input_dir=genomes_path, output_dir=genomes_path)
    #
    # print("Start emapper2GBK")
    # holointeract.network_generation.create_input_mpwt.egg2gbk(
    #     input_dir=genomes_path, output_dir=gbk_files)
    #
    # print("Start metabolic network building")
    # holointeract.network_generation.meta_network.build_network_eggnog(
    #     input_dir=gbk_files, output_dir=metabolic_networks_path, singularity_path=singularity_path)


def metabolic_analysis(community_networks_path, host_networks_path, output_path, seeds, output_name, analysis_method,
                       clustering_method, max_clust, cpu):
    name_assoc = create_abbreviation_names_dict(community_networks_path, host_networks_path, output_path)

    print("Start calculate community scopes")
    community_scopes(community_sbml_path=community_networks_path, hosts_sbml_path=host_networks_path,
                     output_dir=output_path, seeds=seeds, method=analysis_method, name_assoc=name_assoc, cpu=cpu)

    print("Start generating clustermap")
    input_heatmap = os.path.join(output_path, SCOPES_STR, analysis_method)
    output_heatmap = create_heatmap_output(output_path, output_name, analysis_method)
    bact_metabolites, host_metabolites, holo_metabolites, all_metabolites = heatmap_host_bacteria(
        input_dir=input_heatmap, output=output_heatmap, method=clustering_method, max_clust=max_clust)

    output_info = output_heatmap + '_classes_cpd_info'
    proportion_workflow(set(all_metabolites), output=output_info)
    merge_outputs(f'{output_heatmap}_clusters.tsv', f'{output_info}.tsv')

    return name_assoc


def coevolution_analysis(community_networks_path, host_networks_path, output_path, seeds, output_name,
                         clustering_method, max_clust, phylo_tree, correction, cpu):
    scopes_path = os.path.join(output_path, SCOPES_STR, FULL_METHOD)
    if not os.path.exists(scopes_path):
        name_assoc = metabolic_analysis(community_networks_path=community_networks_path,
                                        host_networks_path=host_networks_path, output_path=output_path, seeds=seeds,
                                        output_name=output_name, analysis_method=FULL_METHOD,
                                        clustering_method=clustering_method, max_clust=max_clust, cpu=cpu)
    else:
        name_assoc = load_name_assoc_file(output_path)

    coevolution(scopes_path=scopes_path, output=output_path, name=output_name, name_assoc=name_assoc,
                phylo_tree=phylo_tree, correction=correction)


# MAIN
# ======================================================================================================================
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

    elif args.subcommands == 'coevolution':
        coevolution_analysis(community_networks_path=args.community_networks, host_networks_path=args.host_networks,
                             output_path=args.output, seeds=args.seeds, output_name=args.name,
                             clustering_method=args.clustering_method, max_clust=args.max_clust,
                             phylo_tree=args.phylo_tree, correction=args.correction, cpu=args.cpu)

    else:
        print('[dark_orange]Unknown command. Please use the help (-h) to see available commands.')
        exit(1)

# holointeract metabolic_analysis -cn example/inputs/community/ -hn example/inputs/hosts/ -o example/outputs/ -s example/inputs/seeds/seeds_seawater_artefact.sbml -am solo
# holointeract coevolution -cn example/inputs/community/ -hn example/inputs/hosts/ -o example/outputs/ -s example/inputs/seeds/seeds_seawater_artefact.sbml -p example/inputs/SpeciesTree_rooted.txt
