#!/bin/bash

set -e  # stop if any command fails

echo "Starting BLEvaluator runs..."

#echo "Running ESC..."
#python BLEvaluator.py --config config-files/config_ESC.yaml --borda

#echo "Running hHEP..."
#python BLEvaluator.py --config config-files/config_HEP.yaml --borda

echo "All BLEvaluator runs completed successfully"

echo "Dividing diferent borda methods in folders"
python other_scripts/borda_make_folders.py