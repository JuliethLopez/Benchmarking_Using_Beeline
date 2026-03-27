#!/bin/bash

data=VSC

# move reference files to split folder
mkdir -p $data"_split/"
for file in ClusterIds.csv PseudoTime.csv parameters.txt refNetwork.csv; do
    echo "$file -> $data_split/"
    cp "$file" set$i/
done

# move first 30 folders in ls, sets of 34 folders
t_sets=9   # number of sets you want
for ((i=1; i<=t_sets; i++)); do
    echo $data"_split/set"$i
    mkdir -p $data"_split/set"$i
    ls -d $data/* | head -n 34 | xargs -I {} mv {} $data"_split/set"$i/
done

# copy files
for file in ClusterIds.csv PseudoTime.csv parameters.txt refNetwork.csv; do
    for i in {1..9}; do
        echo "$file -> set$i"
        cd $data"_split/"
        cp "$file" set$i/
    done
done

# remove all folders called tf_hdwgcna
find . -type d -name "tf_hdwgcna" -exec rm -rf {} +

# Print all folders for a datasets
dataset='GSD'

for ncells in 200 1000 2000; do
    for i in {1..100}; do
        echo $dataset'-'$ncells'-'$i
    done
done