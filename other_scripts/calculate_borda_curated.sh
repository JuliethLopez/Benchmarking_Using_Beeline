#!/bin/bash

set -e  # stop if any command fails

echo "Starting BLEvaluator runs..."

echo "Running GSD..."
python BLEvaluator.py --config config-files/config_GSD.yaml --borda

echo "Running VSC..."
python BLEvaluator.py --config config-files/config_VSC.yaml --borda

echo "Running HSC..."
python BLEvaluator.py --config config-files/config_HSC.yaml --borda

echo "Running mCAD..."
python BLEvaluator.py --config config-files/config_mCAD.yaml --borda

echo "All BLEvaluator runs completed successfully."

echo "Dividing diferent borda methods in folders"
python other_scripts/borda_make_folders.py