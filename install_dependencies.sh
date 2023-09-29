#!/bin/bash

pip install -r requirements.txt

mkdir local_packages
cd local_packages/
git clone https://github.com/PaulineGHG/Ontology_sunburst.git

cd Ontology_sunburst/
pip install -e .

cd ../../