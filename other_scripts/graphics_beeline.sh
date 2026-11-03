cd /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/Beeline

# correr en BEELINE
# overview for: undirected, directed, sepr
python BLPlotter-Curated-groups.py --config /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/STREAMLINE/STREAMLINE/config-files/Curated/config_GSD.yaml --auprc --auroc --epr --signed --overview --output outputs_graphics/Curated
python BLPlotter-Curated-groups.py --config /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/STREAMLINE/STREAMLINE/config-files/Curated/config_HSC.yaml --auprc --auroc --epr --signed --overview --output outputs_graphics/Curated
python BLPlotter-Curated-groups.py --config /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/STREAMLINE/STREAMLINE/config-files/Curated/config_mCAD.yaml --auprc --auroc --epr --signed --overview --output outputs_graphics/Curated
python BLPlotter-Curated-groups.py --config /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/STREAMLINE/STREAMLINE/config-files/Curated/config_VSC.yaml --auprc --auroc --epr --signed --overview --output outputs_graphics/Curated
python BLPlotter-Curated-datasets.py --config config-files/config_curated.yaml --overview --output outputs_graphics/Curated

# python BLPlotter-Curated-groups.py --config /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/STREAMLINE/STREAMLINE/config-files/Curated/config_GSD.yaml --overview --output outputs_graphics/Curated
# python BLPlotter-Curated-groups.py --config /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/STREAMLINE/STREAMLINE/config-files/Curated/config_HSC.yaml --overview --output outputs_graphics/Curated
# python BLPlotter-Curated-groups.py --config /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/STREAMLINE/STREAMLINE/config-files/Curated/config_mCAD.yaml --overview --output outputs_graphics/Curated
# python BLPlotter-Curated-groups.py --config /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/STREAMLINE/STREAMLINE/config-files/Curated/config_VSC.yaml --overview --output outputs_graphics/Curated

# tiempo y motivos
Rscript other_scripts/other_plots.R

# experimental
# overview for: epr and sepr
python BLPlotter-Experimental.py --config /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/STREAMLINE/STREAMLINE/config-files/experimental/config_ESC.yaml --epr --overview --output outputs_graphics/experimental
python BLPlotter-Experimental-datasets.py --config /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/STREAMLINE/STREAMLINE/config-files/experimental/config_ESC.yaml --overview --output outputs_graphics/experimental

python BLPlotter-Experimental.py --config /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/STREAMLINE/STREAMLINE/config-files/experimental/config_HEP.yaml --epr --overview --output outputs_graphics/experimental
python BLPlotter-Experimental-datasets.py --config /lustre/scratch127/humgen/projects_v2/jaguar_analysis/analysis/jl38/STREAMLINE/STREAMLINE/config-files/experimental/config_HEP.yaml --overview --output outputs_graphics/experimental