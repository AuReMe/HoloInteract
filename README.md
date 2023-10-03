[![](https://img.shields.io/badge/python-3.9-blue.svg)]()
[![](https://img.shields.io/badge/python-3.10-blue.svg)]()
[![](https://img.shields.io/badge/python-3.11-blue.svg)]()
[![](https://img.shields.io/badge/documentation-unfinished-orange.svg)]()
[![](https://img.shields.io/badge/wiki-nonexistent-red.svg)]()

# **HoloInteract** - Metabolic Interaction in Holobionts

Uses the [Metage2Metabo library]() to generate and manipulate metabolic networks.

## Installation

### Requirements

#### PiPy `requirements.txt` 

- python >=3.9
- setuptools>=65.5.1
- rich>=13.5.3
- matplotlib>=3.8.0
- pandas>=1.5.3
- plotly>=5.17.0
- seaborn>=0.12.2
- Metage2Metabo>=1.5.4
- scipy>=1.11.2
- padmet>=5.0.1
- ete3>=3.1.3
- numpy>=1.26.0
- statsmodels>=0.14.0

#### Local

- ontosunburst : https://github.com/PaulineGHG/Ontology_sunburst.git

### Quick install

In HoloInteract directory :

```bash
bash install_dependencies.sh
pip install -e .
```

## Commands
Subcommands available through `holointeract`:

### Metabolic Analysis : `holointeract metabolic_analysis`

#### Arguments :

- `--comm`, `--community_networks` : path to community networks in SBML
- `--host`, `--host_networks` : path to hosts networks in SBML
- `-o`, `--output` : path to output directory
- `-s`, `--seeds` : path to seeds SBML file
- `-n`, `--name` : output files name (default=run)
- `-m`, `--scopes_method` : method of scopes generation 
`['solo', 'coop', 'full']` (default=`'coop'`)
- `--cm`, `--clustering_method` : method for linkage in clustering (default=`'ward'`) 
- `--max_clust` : maximal number of cluster for the dendrogram division (default=`10`)
- `--cpu` : number of cpu to use (default=`1`)

#### Test example :

```commandline
holointeract metabolic_analysis --comm small_example/inputs/community/ 
--host small_example/inputs/hosts/ -o small_example/outputs/ 
-s small_example/inputs/seeds/seeds_seawater_artefact.sbml -n test -m coop 
```

### Coevolution Analysis : `holointeract coevolution`

#### Arguments :

- `--comm`, `--community_networks` : path to community networks in SBML
- `--host`, `--host_networks` : path to hosts networks in SBML
- `-o`, `--output` : path to output directory
- `-s`, `--seeds` : path to seeds SBML file
- `-n`, `--name` : output files name (default=run)
- `--cm`, `--clustering_method` : method for linkage in clustering (default=`'ward'`) 
- `--max_clust` : maximal number of cluster for the dendrogram division (default=`10`)
- `-p`, `--phylo_tree` : path to phylogenetic tree (Newick format) (default=`None`)
- `--cor`, `--correction` : correction to apply to p-values 
`['bonferroni', 'benjamini', None]` (default=`None`)
- `--cpu` : number of cpu to use (default=`1`)

#### Test example :

```commandline
holointeract coevolution --comm small_example/inputs/community/ 
--host small_example/inputs/hosts/ -o small_example/outputs/ 
-s small_example/inputs/seeds/seeds_seawater_artefact.sbml -n test 
-p small_example/inputs/SpeciesTree_rooted.txt 
```

### Help available

```bash
holointeract -h
```
or
```bash
holointeract 'subcommand' -h
```

## Required Inputs

### Community Networks (`--comm`, `--community_networks`)

Directory containing SBML networks files for each community organism.

Each organism must be placed in the directory named after its natural host.

Host names in community networks directory and host networks directory must coincide.

#### Example :

```
├── Community_networks
│ ├── Host_1
│ │ ├── Microorganism_1.sbml  
│ │ ├── Microorganism_2.sbml  
│ ├── Host_2
│ │ ├── Microorganism_3.sbml  
│ │ ├── Microorganism_4.sbml
```

### Host Networks (`--host`, `--host_networks`)

Directory containing SBML networks files for each host.

Host names in community networks directory and host networks directory must coincide.

#### Example :

```
├── Host_networks
│ ├── Host_1.sbml 
│ ├── Host_2.sbml
```

### Seeds (`-s`, `--seeds`)

Seeds are list of compounds corresponding to a common growth medium shared between all hosts.

Artefacts compounds should be added to seeds to avoid cycles + cofactors.

The seeds file must be given at SBML format.

### Phylogenetic Tree (`-p`, `--phylo_tree`)

The phylogenetic tree is used to calculate phylogenetic distance between each pair of hosts.

The file must be given at Newick format.

All the hosts must be present in the tree.

All the host names must coincide to the names in community networks directory 
and host networks directory.