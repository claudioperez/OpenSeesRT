wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash "Mambaforge-$(uname)-$(uname -m).sh"
exec bash

mamba activate base
mamba install jupyterlab

mamba create -n opensees -c conda-forge -c opensees python=3.9 ipykernel opensees
mamba activate opensees
python -m ipykernel install --user --name opensees --display-name "Python (opensees)"
