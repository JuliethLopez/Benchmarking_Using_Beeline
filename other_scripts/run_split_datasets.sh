# crate config files
dataset=GSD
for i in {1..9}; do
    cp config-files/config_$dataset.yaml config-files/config_$dataset"_"$i.yaml
done

# run jobs
dataset=GSD
for i in {1..9}; do
    echo config_$dataset"_"$i running
    bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q normal -G jaguar_analysis -M10000 -R 'select[mem>10000] rusage[mem=10000]' python BLRunner.py --config config-files/config_$dataset"_"$i.yaml
done

# move results from sets to original folder
dataset=VSC
for set in {1..9}; do
    echo $dataset"_split/set"$set/*" -> "$dataset
    mv $dataset"_split/set"$set/* $dataset
done

# Look for lacking files
# -f files, -s not empty files
for d in outputs/Curated/mCAD/mCAD-2000-*; do
    [ -f "$d/hdWGCNA/rankedEdges.csv" ] || basename "$d"
done

# delete uncomplete folders
for d in outputs/Curated/VSC/VSC-2000-*; do
    [ -f "$d/INFERELATOR/rankedEdges.csv" ] || rm -r "$d/INFERELATOR"
done