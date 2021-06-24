CONDA_ENV = ./env/
CONDA_ACTIVATE = eval "$$(conda shell.bash hook)"; conda activate $(CONDA_ENV); export PYTHONPATH=`pwd`:`pwd`/src:$${PYTHONPATH}

.PHONY: setup

setup: environment.yml
	mamba env create --force --prefix $(CONDA_ENV) --file environment.yml

