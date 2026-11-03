#!/bin/bash

# go to the directory of the script
cd /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline
set -e  # stop if any command fails

echo "Starting BLEvaluator runs..."

echo "Running GSD..."
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config config-files/config_GSD_stability.yaml --spearman

echo "Running VSC..."
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config config-files/config_VSC_stability.yaml --spearman

echo "Running HSC..."
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config config-files/config_HSC_stability.yaml --spearman

echo "Running mCAD..."
bsub -o logs/log.%J.out -e logs/log.%J.err -n 1 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLEvaluator.py --config config-files/config_mCAD_stability.yaml --spearman

echo "All BLEvaluator runs completed successfully."