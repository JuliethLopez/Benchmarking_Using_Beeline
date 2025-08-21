#!/bin/bash

# Building docker for the different algorithms 
echo "This may take a while..."

BASEDIR=$(pwd)

# You may remove the -q flag if you want to see the docker build status
cd $BASEDIR/Algorithms/ARBORETO
docker build -q -t arboreto:base .
if ([ $? = 0 ] && [ "$(docker images -q arboreto:base 2> /dev/null)" != "" ]); then
    echo "Docker container for ARBORETO is built and tagged as arboreto:base"
elif [ "$(docker images -q arboreto:base 2> /dev/null)" != "" ]; then
    echo "Docker container failed to build, but an existing image exists at arboreto:base"
else
    echo "Oops! Unable to build Docker container for ARBORETO"
fi


cd $BASEDIR/Algorithms/JUMP3/
docker build -q -t jump3:base .
if ([ $? = 0 ] && [[ "$(docker images -q jump3:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for JUMP3 is built and tagged as  jump3:base"
elif [ "$(docker images -q jump3:base 2> /dev/null)" != "" ]; then
    echo "Docker container failed to build, but an existing image exists at jump3:base"
else
    echo "Oops! Unable to build Docker container for JUMP3"
fi


cd $BASEDIR/Algorithms/PIDC/
docker build -q -t pidc:base .
if ([ $? = 0 ] && [[ "$(docker images -q pidc:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for PIDC is built and tagged as pidc:base"
elif [ "$(docker images -q pidc:base 2> /dev/null)" != "" ]; then
    echo "Docker container failed to build, but an existing image exists at pidc:base"
else
    echo "Oops! Unable to build Docker container for PIDC"
fi


cd $BASEDIR/Algorithms/PNI/
docker build -q -t pni:base .
if ([ $? = 0 ] && [[ "$(docker images -q pni:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for PNI is built and tagged as pni:base"
elif [ "$(docker images -q pni:base 2> /dev/null)" != "" ]; then
    echo "Docker container failed to build, but an existing image exists at pni:base"
else
    echo "Oops! Unable to build Docker container for PNI"
fi


cd $BASEDIR/Algorithms/PPCOR/
docker build -q -t ppcor:base .
if ([ $? = 0 ] && [[ "$(docker images -q ppcor:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for PPCOR is built and tagged as ppcor:base"
elif [ "$(docker images -q ppcor:base 2> /dev/null)" != "" ]; then
    echo "Docker container failed to build, but an existing image exists at ppcor:base"
else
    echo "Oops! Unable to build Docker container for PPCOR"
fi


#cd $BASEDIR/Algorithms/SCSGL/
#docker build -q -t scsgl:base .
#if ([ $? = 0 ] && [[ "$(docker images -q scsgl:base 2> /dev/null)" != "" ]]); then
#    echo "Docker container for SCSGL is built and tagged as scsgl:base"
#elif [ "$(docker images -q scsgl:base 2> /dev/null)" != "" ]; then
#    echo "Docker container failed to build, but an existing image exists at scsgl:base"
#else
#    echo "Oops! Unable to build Docker container for SCSGL"
#fi


cd $BASEDIR/Algorithms/INFERELATOR/
docker build -t inferelator:base .
if ([ $? = 0 ] && [[ "$(docker images -q inferelator:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for INFERELATOR is built and tagged as inferelator:base"
elif [ "$(docker images -q inferelator:base 2> /dev/null)" != "" ]; then
    echo "Docker container failed to build, but an existing image exists at inferelator:base"
else
    echo "Oops! Unable to build Docker container for INFERELATOR"
fi


cd $BASEDIR/Algorithms/TF_hdWGCNA/
docker build -q -t tf_hdwgcna:base .
if ([ $? = 0 ] && [[ "$(docker images -q tf_hdwgcna:base 2> /dev/null)" != "" ]]); then
    echo "Docker container for TF_hdWGCNA is built and tagged as tf_hdwgcna:base"
elif [ "$(docker images -q tf_hdwgcna:base 2> /dev/null)" != "" ]; then
    echo "Docker container failed to build, but an existing image exists at tf_hdwgcna:base"
else
    echo "Oops! Unable to build Docker container for TF_hdWGCNA"
fi