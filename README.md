<img align="left" src="https://github.com/BRACE2/OpenSeesRT/blob/master/docs/figures/peer-black.svg" width="250px">

# G3 / OpenSees

<br>

--------------------------------------------------------------------

<br>

<div style="align:center">

<!--
[![PyPI Downloads](https://img.shields.io/pypi/dm/opensees?style=for-the-badge)](https://pypi.org/project/opensees)
[![Latest conda-forge version](https://img.shields.io/conda/vn/conda-forge/opensees?logo=conda-forge&style=for-the-badge)](https://anaconda.org/conda-forge/opensees)
-->

[![Latest PyPI version](https://img.shields.io/pypi/v/opensees?logo=pypi&style=for-the-badge)](https://pypi.python.org/pypi/opensees)
[![](https://img.shields.io/conda/v/opensees/opensees?color=%23660505&style=for-the-badge)](https://anaconda.org/opensees/opensees)

</div>

--------------------------------------------------------------------


<!--

## Installing / Compiling

### I. Simple user install

1. Install `mambaforge`
2. `mamba install -c opensees opensees`

### II. Editable (developer) install

1. Install `mambaforge`
2. `pip install -e .`
3. To recompile : `python setup.py build_ext`

### III. Distribution build

Pre-requisites:
1. Install `conda-build`

Steps for `opensees`:
1. `conda-build -c conda-forge etc/conda/ --python 3.7`
2. `anaconda upload <path>`

Steps for `opensees-intel`
1. `conda-build -c intel etc/conda-intel/ --python 3.7`
2. `anaconda upload <path>`

### IV. Python-only build

```
python setup.py [install|develop] --skip-cmake
```

## Build Environment

```shell
  conda create -c conda-forge -n skbuild python==3.8 cmake'>=3.18'  \
    scikit-build pybind11 setuptools
```

-->


## Installing
- For Windows, first activate WSL2 and install a Linux distribution.

- Install `mambaforge`, a small Anaconda distribution.
  - Go to [https://github.com/conda-forge/miniforge#mambaforge](https://github.com/conda-forge/miniforge#mambaforge)
  - Click the `Mambaforge-Linux-x86_64` link to download an install script.
  - Run the downloaded script.

  This can all be done at once with the following commands:

  ```bash
  wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
  bash Mambaforge-$(uname)-$(uname -m).sh
  ```

- Create an environment to run `opensees` and `jupyter-lab` out of
  
  ```bash
  conda create -n opensees -c conda-forge python=3.9 jupyterlab matplotlib numpy scipy pyyaml

  # the --skip-deps flag ensures the packages from conda-forge arent clobbered by PyPI ones.
  pip install opensees --skip-deps
  ```

Once this is done, you can run the commands

```bash
conda activate opensees
jupyter-lab
```

from the WSL terminal, and a url will be printed which you can open in a browser to work in
Jupyter Lab.

<!-- 
  - Install Jupyterlab

    ```bash
    mamba create -n jupyter jupyterlab 
    mamba activate opensees
    python -m ipykernel install --user --name opensees --display-name "Python (opensees)"
    ```
-->



                           ┌─┐┌─┐┌─┐┌─┐  ┌──┌─┐┌─┐ ┌──
                           └─┘├─┘└──┘ │ ─┘  └──└───┘
           ───────────────────┘Berkeley, California ──────────────────────
                                   © UC Regents

