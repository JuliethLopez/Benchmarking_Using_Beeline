#!/bin/bash
echo "Generating graphics for Curated networks..."

# overview, boxplots
python BLPlotter-Curated-groups.py --config config-files/config_mCAD.yaml -o /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated -v -a -r -e
python BLPlotter-Curated-groups.py --config config-files/config_GSD.yaml -o /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated -v -a -r -e
python BLPlotter-Curated-groups.py --config config-files/config_VSC.yaml -o /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated -v -a -r -e
python BLPlotter-Curated-groups.py --config config-files/config_HSC.yaml -o /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated -v -a -r -e

# general overview per dataset
python BLPlotter-Curated.py --config config-files/config_mCAD.yaml -o /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated -v
python BLPlotter-Curated.py --config config-files/config_GSD.yaml -o /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated -v
python BLPlotter-Curated.py --config config-files/config_VSC.yaml -o /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated -v
python BLPlotter-Curated.py --config config-files/config_HSC.yaml -o /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated -v

# random values for AUROC are 0.5, but for AUPRC and EPr depend on the density, that is calculated when written 0
python BLPlotter-Curated-datasets.py --config config-files/config_curated.yaml -o /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline/outputs_graphics/Curated -v
