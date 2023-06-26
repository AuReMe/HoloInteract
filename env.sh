#!/bin/sh
. /local/env/envconda.sh


WD=$(pwd)
export CONDA_ALWAYS_YES="true"

# init env
# conda create -p $WD"/.env_eggnog" python=3.8

# activate env to install packages
# conda activate $WD"/.env_eggnog"


# python -m pip install --upgrade pip

# python -m pip install metage2metabo

# python -m pip install --force-reinstall -v "eggnog-mapper" 

# conda deactivate



#2e
conda create -p $WD"/.env_metabo" python=3.8
conda activate $WD"/.env_metabo"

python -m pip install --upgrade pip

python -m pip install . 



# installing required python packages





unset CONDA_ALWAYS_YES
conda deactivate