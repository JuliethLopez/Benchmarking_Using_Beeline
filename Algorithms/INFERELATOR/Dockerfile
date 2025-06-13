#basado en ejemplo de arboreto
FROM continuumio/anaconda3:2024.02-1

LABEL maintainer = "Método añadido por Julieth López <julalopezcas@unal.edu.co>"

USER root 

WORKDIR /

COPY docs/ /docs/
COPY docs/environment.yml /docs/environment.yml
#COPY docs/setup.py docs/setup.py

RUN conda env create -f docs/environment.yml

# Use conda run to install inferelator in development mode
#RUN conda run -n inferenv python docs/setup.py develop

# OPTIONAL: pip install version if not using setup.py directly
RUN conda run -n inferenv python -m pip install inferelator

COPY runINFERELATOR.py /

RUN mkdir data/

RUN apt-get update

RUN apt-get install time
