#!/bin/bash
mkdir -p logs

# Curated (2 CPUs)
bsub -o logs/log.%J.out -e logs/log.%J.err -n 2 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLRunner.py --config config-files/config_GSD_hdWGCNA.yaml
bsub -o logs/log.%J.out -e logs/log.%J.err -n 2 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLRunner.py --config config-files/config_mCAD_hdWGCNA.yaml
bsub -o logs/log.%J.out -e logs/log.%J.err -n 2 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLRunner.py --config config-files/config_HSC_hdWGCNA.yaml
bsub -o logs/log.%J.out -e logs/log.%J.err -n 2 -q long -G jaguar_analysis -M20000 -R 'select[mem>20000] rusage[mem=20000]' python BLRunner.py --config config-files/config_VSC_hdWGCNA.yaml

# Experimental (4 CPUs)
bsub -o logs/log.%J.out -e logs/log.%J.err -n 4 -q long -G jaguar_analysis -M40000 -R 'select[mem>40000] rusage[mem=40000]' python BLRunner.py --config config-files/config_ESC_hdWGCNA.yaml
bsub -o logs/log.%J.out -e logs/log.%J.err -n 4 -q long -G jaguar_analysis -M40000 -R 'select[mem>40000] rusage[mem=40000]' python BLRunner.py --config config-files/config_HEP_hdWGCNA.yaml

# Job <321438> is submitted to queue <long>.
# Job <321439> is submitted to queue <long>. klilled, running<481796>
# Job <321440> is submitted to queue <long>.
# Job <321441> is submitted to queue <long>.

# Finished! :)
# Job <321442> is submitted to queue <long>.
# Job <321443> is submitted to queue <long>.