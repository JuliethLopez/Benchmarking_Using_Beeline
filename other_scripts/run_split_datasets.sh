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

# move results from sets to original folder ,PPCOR
dataset=GSD
method=GRNBOOST2
for set in {1..100}; do
    for ncells in {2000}; do
        echo "outputs/Curated_25_04_2026/"$dataset"/"$dataset"-"$ncells"-"$set"/"$method" -> outputs/Curated/"$dataset"/"$dataset"-"$ncells"-"$set"/"
        mv "outputs/Curated_25_04_2026/"$dataset"/"$dataset"-"$ncells"-"$set"/"$method "outputs/Curated/"$dataset"/"$dataset"-"$ncells"-"$set"/"
    done
done


# Look for lacking files
# -f files, -s not empty files
#HSC,VSC,GSD,mCAD
#PPCOR,PIDC,GRNBOOST2,INFERELATOR,hdWGCNA
for d in outputs/Curated/HSC/HSC-*; do
    for method in {PPCOR,PIDC,GRNBOOST2,INFERELATOR,hdWGCNA}; do
        [ -f "$d/"$method"/rankedEdges.csv" ] || basename "$d/"
    done
done

for d in outputs/Curated/mCAD/mCAD-*; do
    for method in INFERELATOR; do
        [ -f "$d/"$method"/not_imputation_not_prior/network.tsv.gz" ] || basename "$d/"
    done
done

for d in outputs/Curated/GSD/GSD-*; do
    f="$d/INFERELATOR/not_imputation_not_prior/network.tsv.gz"
    [ -f "$f" ] || basename "$d"
done

for d in outputs/Curated/GSD/GSD-*/INFERELATOR/not_imputation_not_prior; do
    if ls "$d"/random_seed* 1> /dev/null 2>&1; then
        echo "$d"
    fi
done


# delete uncomplete folders
for d in outputs/Curated/mCAD/mCAD-*; do
    [ -f "$d/INFERELATOR/not_imputation_not_prior/crossvalidation_performance.tsv" ] || rm -r "$d/INFERELATOR/not_imputation_not_prior/crossvalidation_performance.tsv"
done

# delete the whole folder
#& rm -r "$d/hdWGCNA"
for d in outputs/Curated/HSC/HSC-*; do
    rm -r "$d/INFERELATOR" 
done