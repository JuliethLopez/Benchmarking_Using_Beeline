#!/bin/bash

# go to the directory of the script
cd /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline
set -e  # stop if any command fails

echo "Starting BLEvaluator runs Curated..."

echo "Running GSD..."
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config config-files/config_GSD.yaml --time --memory --motifs --sepr
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config ../STREAMLINE/STREAMLINE/config-files/Curated/config_GSD.yaml --jaccard  

echo "Running VSC..."
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config config-files/config_VSC.yaml --time --memory --motifs --sepr
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config ../STREAMLINE/STREAMLINE/config-files/Curated/config_VSC.yaml --jaccard

echo "Running HSC..."
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config config-files/config_HSC.yaml --time --memory --motifs --sepr
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config ../STREAMLINE/STREAMLINE/config-files/Curated/config_HSC.yaml --jaccard

echo "Running mCAD..."
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config config-files/config_mCAD.yaml --time --memory --motifs --sepr
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config ../STREAMLINE/STREAMLINE/config-files/Curated/config_mCAD.yaml --jaccard

echo "All BLEvaluator runs completed successfully for Curated datasets."

echo "Starting BLEvaluator runs Experimental..."

echo "Running ESC..."
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config config-files/config_ESC_sepr.yaml --sepr

echo "Running HEP..."
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config config-files/config_HEP_sepr.yaml --sepr

