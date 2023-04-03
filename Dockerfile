FROM continuumio/miniconda3

RUN apt update

#RUN conda create -n vina python=3
RUN conda install -y -c conda-forge numpy swig boost-cpp sphinx sphinx_rtd_theme

RUN apt update
RUN apt install -y git
RUN git clone https://github.com/ccsb-scripps/AutoDock-Vina
WORKDIR /AutoDock-Vina


RUN apt update
RUN apt install -y build-essential
#RUN apt install -y build-essentials
RUN apt-get install -y manpages-dev
RUN apt-get install -y libboost-all-dev swig

WORKDIR /AutoDock-Vina/build/python
RUN python setup.py build install

RUN pip install vina
RUN conda install -y -c conda-forge pdbfixer
RUN pip install deepchem
RUN pip install minio
