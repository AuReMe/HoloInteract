[![](https://img.shields.io/badge/python-3.9-blue.svg)]()
[![](https://img.shields.io/badge/python-3.10-blue.svg)]()
[![](https://img.shields.io/badge/python-3.11-blue.svg)]()
[![](https://img.shields.io/badge/documentation-unfinished-orange.svg)]()
[![](https://img.shields.io/badge/wiki-nonexistent-red.svg)]()

# **HoloInteract** - Metabolic Interaction in Holobionts

Uses the [Metage2Metabo library]() to generate and manipulate metabolic networks.

## Installation

Requires **python >=3.9**.

```bash
git clone git@github.com:CoLucas22/HoloInteract
cd HoloInteract
bash env.sh
```

## Quick start : provided commands

First you have to load the conda environment

```bash
conda activate .env_metabo/
```

Subcommands available through `holointeract`:

- **get_alg_scope** : Compute the scope of hosts.
- **community_analysis** : Compute the scope of the bacteria community.
- **stat_cpd** : Compute statistics on the scopes matrix.
- **clustermap** : Generates a clustermap of the scopes matrix.
- **coevolution_scopes** : Computes the scopes of the bacterias in cooperation with each algae.
- **coevolution_matrix** : Generates a matrix of the coevolution scopes.
- **coevolution_graph** : Generates a graph of the metabolic. complementarity as a function of phylogenetic distance from the host organism.

These commands below launch a full workflow :

- **metabolic_analysis** All the scope analysis of holobionts.

  - **get_alg_scope** : Compute the scope of hosts.
  - **community_analysis** : Compute the scope of the bacteria community.
  - **stat_cpd** : Compute statistics on the scopes matrix.
  - **clustermap** : Generates a clustermap of the scopes matrix.

- **coevolution_analysis** All the scope analysis of holobionts and the coevolution research.
  - **get_alg_scope** : Compute the scope of hosts.
  - **community_analysis** : Compute the scope of the bacteria community.
  - **stat_cpd** : Compute statistics on the scopes matrix.
  - **clustermap** : Generates a clustermap of the scopes matrix.
  - **coevolution_scopes** : Computes the scopes of the bacterias in cooperation with each algae.
  - **coevolution_matrix** : Generates a matrix of the coevolution scopes.
  - **coevolution_graph** : Generates a graph of the metabolic. complementarity as a function of phylogenetic distance from the host organism.

### Required :

You have to put a '/' at the end of the paths you give as an argument of the function.
ex:

```bash
holointeract get_alg_scopes -r /repository/to/host/networks/ --seeds /path/to/seeds.sbml -o /path/to/out/directory/
```

The repository structure of the bacterias metabolic networks has to be :

├── Bacteria_networks/  
│ ├── Host1/  
│ │ ├── sbml/  
│ │ │ ├── Bacteria1.sbml  
│ │ │ ├── Bacteria2.sbml  
│ ├── Host2/  
│ │ ├── sbml/  
│ │ │ ├── Bacteria3.sbml  
│ │ │ ├── Bacteria4.sbml

The repository structure of the host metabolic networks has to be :

├── Host_networks/  
│ ├── Host1.sbml  
│ ├── Host2.sbml

## Reconstruct metabolic networks :

Need a singularity image with Pathway Tools installed and use the Genouest Bioinformatics Cluster.

```bash
. /local/env/enveggnog-mapper-2.1.9.sh

python python_scripts/network_from_genomes.py /path/to/bact/genomes/ /path/to/gbk/files/ /path/to/bact/networks/ /path/to/singularity
```
