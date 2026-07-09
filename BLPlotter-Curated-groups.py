#!/usr/bin/env python
# coding: utf-8

import os
import sys
from turtle import left
import yaml
import argparse
import itertools
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
import multiprocessing
from pathlib import Path
import concurrent.futures
from itertools import permutations
from collections import defaultdict
from multiprocessing import Pool, cpu_count
from networkx.convert_matrix import from_pandas_adjacency
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.lines as line
import matplotlib.patches as patches

# local imports
import BLEval as ev

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description='Generate plots from evaluation results.')
    parser.add_argument('-c','--config', default='config.yaml', help="Config files.")
    parser.add_argument('-a', '--auprc', action="store_true", default=False, help="AUPRC box plot.")
    parser.add_argument('-r', '--auroc', action="store_true", default=False, help="AUROC box plot.")
    parser.add_argument('-e', '--epr', action="store_true", default=False, help="Early precision box plot.")
    parser.add_argument('-v', '--overview', action="store_true", default=False, help="Overview plot.")
    parser.add_argument('-o', '--output', default='.', help="Output directory.")
    return parser

def parse_arguments():
    opts = get_parser().parse_args()
    return opts

def boxplot(opts, evalConfigs, datasets, randValue, resTypeFile, resTypeName):
    plt.rcParams.update({'font.size': 14})
    DFs = []

    for i, dataset in enumerate(datasets):
        evalConfig = evalConfigs[i]
        DF = pd.read_csv(str(evalConfig.output_settings.base_dir) + '/' \
                           + str(evalConfig.input_settings.datadir).split("inputs")[1] + '/' \
                           + str(evalConfig.output_settings.output_prefix) \
                           + '-' + resTypeFile + '.csv', header = 0, index_col = 0)
        DF = DF.T
        DF.index.name = 'dataset'
        DF.reset_index(inplace = True)  
        DFs.append(DF)  

    f, axes = plt.subplots(len(datasets), 1, figsize=(1.5*max([len(DF.columns) for DF in DFs]), 5*len(datasets)))
     
    for i, dataset in enumerate(datasets):
        DF = DFs[i]

        if "undirected" in resTypeFile:
            # value_vars=['PIDC', 'GENIE3', 'GRNBOOST2', 'PPCOR', 'SCODE', 'GRISLI', 'SINGE', 'SCNS', 'LEAP', 'SINCERITIES', 'GRNVBEM', 'SCRIBE', 'SCSGL'])
            # value_vars=['SINCERITIES', 'SCRIBE', 'SINGE', 'PPCOR', 'PIDC', 'GENIE3', 'LEAP', 'GRNBOOST2', 'GRISLI', 'GRNVBEM', 'SCNS', 'SCODE', 'SCSGL'])
            value_vars=['PIDC', 'GRNBOOST2', 'PPCOR', 'INFERELATOR', 'hdWGCNA', 'BORDA', 'sBORDA', 'mBORDA', 'smBORDA']
        else:
            value_vars=['GRNBOOST2', 'INFERELATOR']
        
        modifiedDF = pd.melt(DF, id_vars = ['dataset'], value_vars = value_vars)
        
        conditions = [
            modifiedDF['dataset'].str.contains('-200-'),
            modifiedDF['dataset'].str.contains('-1000-'),
            modifiedDF['dataset'].str.contains('-2000-')
        ]
        outputs = ['200', '1000', '2000']
        modifiedDF['category'] = pd.Series(np.select(conditions, outputs, 'Other'))
        existing_cats = ['200', '1000', '2000']
        modifiedDF['category'] = modifiedDF['category'].astype('category').cat.set_categories(existing_cats, ordered=True)
              
        if len(datasets) > 1: subax = axes[i]
        else: subax= axes
        
        ax = sns.lineplot(y = randValue, x = range(-1,len(DF.columns)+1), ax = subax)
        ax.lines[0].set_linestyle("--")
        ax.lines[0].set_color("gray")
        
        ax = sns.boxplot(y = 'value', x = 'variable', data = modifiedDF, fliersize = 0, hue = "category", palette = sns.color_palette("Set1"), ax = subax)
        subax.set_xticklabels(subax.get_xticklabels(), rotation=45, ha='right')
        
        handles, labels = ax.get_legend_handles_labels()
        sns.swarmplot(y ='value', x = 'variable', data = modifiedDF, hue = "category", alpha = 0.5, palette = ['k'], dodge = True, ax = subax)
        ax.legend_.remove()
        ax.legend(handles, labels, title = 'Number of cells', bbox_to_anchor=(.5, 1), loc = 'upper center', ncol = 3)
        ax.set_ylim([0.0,1])
        ax.set_ylabel(resTypeName, fontsize = 18)
        ax.set_title(resTypeName +' values for ' + dataset +' Network', fontsize = 18)
        plt.tight_layout()
    
    subax.set_xlabel('Algorithm', fontsize = 18)
    file = opts.output + '/' + '-'.join([str(c.output_settings.output_prefix) for c in evalConfigs]) + '-boxplot-' + resTypeFile + '.pdf'
    print("Boxplot saved to " + file)
    plt.savefig(file, dpi = 300)

def COplot(inputDF, width = 12, height = 7, randValues = [], shape = [], 
            palettes = [], levels = [], rotation = [], switch = [], metric_type = 'ratio'):

    rowNames = inputDF.index
    maxRows = len(rowNames)
    # Total column calculation remains to handle subgroups correctly
    total_cols = sum([len(inputDF[lvl].columns) for lvl in levels])
    pad = 2
    fSize = (width, height)

    f = plt.figure(figsize=fSize)

    # Main plot area setup
    ax = plt.subplot2grid((2 * maxRows + 3, len(levels)), (0, 0), colspan=len(levels), rowspan=2 * maxRows + 2)
    ax.set_yticks(np.arange(0, maxRows + pad))
    ax.set_xticks(np.arange(0, total_cols + pad))
    
    Indices = [""] + list(rowNames) + [""]
    ax.set_yticklabels(Indices, fontsize=18)

    for spine in ax.spines.values():
        spine.set_visible(False)

    # Zebra stripes background
    alt = False
    for rowIdx in range(maxRows):
        rect_color = (0.9, 0.9, 0.9) if alt else (1, 1, 1)
        ptch = patches.Rectangle((0, rowIdx + 0.5), width=total_cols + 1, height=1,  
                                 edgecolor=(1,1,1), facecolor=rect_color, zorder=-1)
        alt = not alt
        ax.add_artist(ptch)

    col_ptr = 0 
    for levlIdx in range(len(levels)):
        lvl = levels[levlIdx]
        # Use 1.0 as default random value for ratios if not provided
        randValue = randValues[levlIdx] if levlIdx < len(randValues) else 1.0
        
        # Access the specific palette for this metric
        current_pal = palettes[levlIdx]
        pal_len = len(current_pal) # CHANGE: Dynamically get palette length

        subcols = inputDF[lvl].columns
        n_sub = len(subcols)

        # Global metric scaling (Metric specific min/max)
        # This allows raw counts to be mapped correctly to the color scale
        v_min_global = inputDF[lvl].min().min()
        v_max_global = inputDF[lvl].max().max()

        # Plot Metric Header
        plt.text(col_ptr + (n_sub / 2) + 0.5, maxRows + 2, 
                 lvl, fontsize=18, fontweight='bold', ha="center", va="center")

        for sc_idx, sc in enumerate(subcols):
            # Subgroup header (e.g., cell count labels)
            plt.text(col_ptr + sc_idx + 1, maxRows + 1, sc,
                     fontsize=14, rotation=rotation[levlIdx],
                     ha="center", va="center")

            for rowIdx in range(maxRows):
                value = inputDF.loc[rowNames[rowIdx], (lvl, sc)]
                x_pos = col_ptr + sc_idx + 1
                y_pos = rowIdx + 1

                if np.isnan(value):
                    plt.text(x_pos, y_pos, '-', fontsize=24, ha="center", va="center", color='#2E4053')
                    continue

                # --- COLOR MAPPING ---
                # Normalize raw value to [0, 1] relative to the metric's global range
                if v_max_global != v_min_global:
                    norm_val = (value - v_min_global) / (v_max_global - v_min_global)
                else:
                    norm_val = 0.5
                
                # Assign color: Dark gray if below random threshold, else map to palette
                if value < randValue:
                    col = '#2E4053'
                else:
                    # Use dynamic palette length - 1 to find the color index
                    palette_idx = int(np.clip(norm_val * (pal_len - 1), 0, pal_len - 1))
                    col = current_pal[palette_idx]

                # --- SHAPE RENDERING ---
                if shape[levlIdx] == 'arrow':
                    # Directional arrows based on ratio vs random
                    txt = '\u2197' if value > randValue else ('\u2198' if value < randValue else '=')
                    txt_col = '#4CAF50' if value > randValue else ('#E53935' if value < randValue else '#2E4053')
                    plt.text(x_pos, y_pos, txt, fontsize=24, ha="center", va="center", color=txt_col)

                elif shape[levlIdx] == 'text':
                    txt = round(value, 2) if value < 10 else int(value)
                    plt.text(x_pos, y_pos, txt, fontsize=11, ha="center", va="center", color='#2E4053')

                else:
                    # Scaling shapes visually based on norm_val
                    patch = None
                    if shape[levlIdx] == 'c':
                        patch = patches.Circle((x_pos, y_pos), radius=np.sqrt(norm_val)/2.5, facecolor=col, edgecolor='k')
                    elif shape[levlIdx] == 's':
                        size = norm_val * 0.8
                        patch = patches.Rectangle((x_pos - size/2, y_pos - size/2), size, size, facecolor=col, edgecolor='k')
                    elif shape[levlIdx] == 'hm':
                        patch = patches.Rectangle((x_pos - 0.45, y_pos - 0.5), 0.9, 1, facecolor=col, edgecolor=col)
                    elif shape[levlIdx] == 'f':
                        patch = patches.Rectangle((x_pos - 0.5, y_pos - 0.5), 1, 1, facecolor=col, edgecolor=col)
                    
                    if patch:
                        ax.add_artist(patch)

                    # Text Overlay for full-cell shapes (hm/f)
                    if shape[levlIdx] in ['hm', 'f']:
                        # switch[levlIdx] determines if high values get white text for contrast
                        textCol = 'white' if switch[levlIdx] and norm_val > 0.6 else 'black'
                        plt.text(x_pos, y_pos, round(value, 1), fontsize=10, ha="center", va="center", color=textCol)

        col_ptr += n_sub 

    # Clean up main axis
    ax.yaxis.set_ticks_position('none') 
    ax.xaxis.set_ticks_position('none') 
    ax.set_xticklabels([])
    
    # --- LEGEND SCALES ---
    for levlIdx in range(len(levels)):
        # not legend for arrow shape
        if shape[levlIdx] == 'arrow':
            continue
        
        # Add color range legend
        lvl = levels[levlIdx]
        v_min = inputDF[lvl].min().min()
        v_max = inputDF[lvl].max().max()
        
        ax_leg = plt.subplot2grid((2*maxRows + 3, len(levels)), (2*maxRows + 2, levlIdx))
        plt.subplots_adjust(left=0.1, right=0.9, bottom=0.15)
        
        # Display the palette as a continuous image
        ax_leg.imshow(np.arange(len(palettes[levlIdx])).reshape(1, len(palettes[levlIdx])),
                  cmap=mpl.colors.ListedColormap(list(palettes[levlIdx])),
                  interpolation="nearest", aspect="auto")        
        
        ax_leg.set_yticks([])
        # The legend labels now represent the raw count/ratio min/max
        ax_leg.set_xticks([0, len(palettes[levlIdx]) - 1])
        ax_leg.set_xticklabels([f'{v_min:.1f}', f'{v_max:.1f}'], fontsize=12)
        ax_leg.set_title(f"Scale: {lvl}", fontsize=10)

    return f


def main():
    opts = parse_arguments()
    config_files = opts.config.split(',')

    # 1. CARGA DE CONFIGURACIÓN
    datasets = []
    evalConfigs = []

    for config_file in config_files:
        with open(config_file, 'r') as conf:
            evalConfig = ev.ConfigParser.parse(conf)
            evalConfigs.append(evalConfig)
            # Obtenemos el nombre del dataset
            datasets.append(str(os.path.basename(evalConfig.input_settings.datadir)))

    # 2. CÁLCULO DE RANDOM PREDICTORS
    randPredictor = {}
    for i, dataset in enumerate(datasets):
        evalConfig = evalConfigs[i]
        true_path = os.path.join(evalConfig.input_settings.datadir, 
                                 evalConfig.input_settings.datasets[0]['name'], 
                                 evalConfig.input_settings.datasets[0]['trueEdges'])
        trueEdgesDF = pd.read_csv(true_path).drop_duplicates()
        trueEdgesDF = trueEdgesDF[trueEdgesDF.Gene1 != trueEdgesDF.Gene2]
        numGenes = len(np.unique(trueEdgesDF[['Gene1', 'Gene2']]))
        randPredictor[dataset] = len(trueEdgesDF) / (numGenes * (numGenes - 1))

    # Generate AUPRC boxplot     
    if (opts.auprc):
        print('\n\nGenerating AUPRC boxplot...')
        boxplot(opts, evalConfigs, datasets, randPredictor[dataset], 'AUPRC-undirected', 'AUPRC')
        boxplot(opts, evalConfigs, datasets, randPredictor[dataset], 'AUPRC-directed', 'AUPRC')

    # Generate AUROC boxplot     
    if (opts.auroc):
        print('\n\nGenerating AUROC boxplot...')
        boxplot(opts, evalConfigs, datasets, 0.5, 'AUROC-undirected', 'AUROC')
        boxplot(opts, evalConfigs, datasets, 0.5, 'AUROC-directed', 'AUROC')
        
    # Generate EPr boxplot     
    if (opts.epr):
        print('\n\nGenerating Early Precision boxplot...')
        boxplot(opts, evalConfigs, datasets, randPredictor[dataset], 'EPr-undirected', 'Early Precision')
        boxplot(opts, evalConfigs, datasets, randPredictor[dataset], 'EPr-directed', 'Early Precision')
    # 3. GENERACIÓN DEL OVERVIEW HEATMAP
    if opts.overview:
        print('\n\nGenerating Grouped Overview Heatmap...')
        
        # --- EXTRACCIÓN DE RUTAS DESDE CONFIG (CORREGIDO) ---
        # Convertimos el objeto PosixPath a string para poder usar split
        datadir_str = str(evalConfigs[0].input_settings.datadir)
        
        # Extraemos la subcarpeta después de "inputs"
        # Esto maneja casos como "inputs/Curated/mCAD" -> "Curated/mCAD"
        if "inputs/" in datadir_str:
            sub_dir = datadir_str.split("inputs/")[1]
        elif "inputs" in datadir_str:
            sub_dir = datadir_str.split("inputs")[1].lstrip('/')
        else:
            sub_dir = datadir_str # Fallback por si la ruta no contiene "inputs"

        out_base = evalConfigs[0].output_settings.base_dir
        db_prefix = evalConfigs[0].output_settings.output_prefix
        
        # Reconstrucción de la ruta: outputs/Curated/mCAD/
        full_output_path = Path(out_base) / sub_dir
        
        # Define metrics to plot (for the plot name) and their suffixes (the real name)
        #resType = ['AUPRC Ratio', 'AUROC Ratio', 'Early Precision Ratio'] #pale3, pale3h,  pale3r, randvalues 0,0.5,0
        #resSuffixes = ['AUPRC-undirected', 'AUROC-undirected', 'EPr-undirected']
        #resType = ['AUPRC Ratio', 'AUROC Ratio', 'Early Precision Ratio'] #pale3, pale3h, pale3r, randvalues 0,0.5,0
        #resSuffixes = ['AUPRC-directed', 'AUROC-directed', 'EPr-directed']
        #resType = ['Early Precision Ratio (Inhibitory)', 'Early Precision Ratio (Activation)'] #pale3, pale3h, randvalues 0,0
        #resSuffixes = ['EPr-Inhibitory', 'EPr-Activation']
        resType = ['Motif-FBL', 'Motif-FFL', 'Motif-MI'] #arrow, arrow, arrow, randvalues 0 #MEJOR USAR OTHER_PLOTS.R o CUrated-sin grupos
        resSuffixes = ['NetworkMotifs-FBL', 'NetworkMotifs-FFL', 'NetworkMotifs-MI']
        #resType = ['Memory', 'Time'] #pale5r, pale5, randvalues 0, 0 #MEJOR USAR OTHER_PLOTS.R
        #resSuffixes = ['Memory', 'Times']
        
        cell_groups = ['200', '1000', '2000']
        
        # Nombres de archivos finales
        resTypeFileName = [f"{db_prefix}-{suffix}" for suffix in resSuffixes]
        
        # Estructura MultiIndex
        multIndTuple = [(res, grp) for res in resType for grp in cell_groups]
        
        first_file = full_output_path / f"{resTypeFileName[0]}.csv"
        
        if not first_file.exists():
            print(f"Error: No se encontró el archivo en {first_file}")
            return
            
        sample_df = pd.read_csv(first_file, index_col=0)
        overviewDF = pd.DataFrame(index=sample_df.index, columns=pd.MultiIndex.from_tuples(multIndTuple))

        # Ingesta de datos
        for i, res in enumerate(resType):
            current_file = full_output_path / f"{resTypeFileName[i]}.csv"
            if current_file.exists():
                raw_res = pd.read_csv(current_file, index_col=0)
                for group in cell_groups:
                    # Buscamos de forma más flexible (que contenga el número)
                    cols = [c for c in raw_res.columns if str(group) in c]
                    if cols:
                        # Calculamos la mediana de las columnas encontradas
                        values = raw_res[cols].median(axis=1)
                        # Dividimos por el random predictor si es una métrica de ratio, sino dejamos el valor crudo
                        if 'Ratio' in res:
                            print(f"Metric: {res}, group: {group}")
                            print(values)
                            print(randPredictor[datasets[0]])
                            overviewDF.loc[:, (res, group)] = values / randPredictor[datasets[0]]
                            print(overviewDF.loc[:, (res, group)])
                        elif 'Stability' in res:
                            # Guidance uses the first value for stability
                            overviewDF.loc[alg, (res, group)] = values.iloc[0]
                        elif res == 'Memory':
                            overviewDF.loc[:, (res, group)] = values / 1048576 # Convertir a GB
                        else:
                            overviewDF.loc[:, (res, group)] = values.median()


                    else:
                        print(f"Advertencia: No se encontraron columnas para el grupo {group} en {res}")

        # Limpieza y Orden
        overviewDF = overviewDF.dropna(how='all')
        #overviewDF = overviewDF.sort_values(by=(resType[0], cell_groups[-1]), ascending=True)
        #undirected
        #overviewDF = overviewDF.reindex(['PIDC', 'GRNBOOST2', 'PPCOR', 'INFERELATOR', 'hdWGCNA', 'BORDA', 'sBORDA', 'mBORDA', 'smBORDA']).iloc[::-1]
        # memory, time
        #overviewDF = overviewDF.reindex(['PIDC', 'GRNBOOST2', 'PPCOR', 'INFERELATOR', 'hdWGCNA']).iloc[::-1]
        #directed, motifs
        overviewDF = overviewDF.reindex(['GRNBOOST2', 'INFERELATOR']).iloc[::-1]
        print("DataFrame para Overview Heatmap (limpio y ordenado):")
        print(overviewDF)
            
        print("DataFrame para Overview Heatmap:")
        #print(overviewDF)
        # Paletas y Plot (Usando la función COplot con nombres a la izquierda y escalas abajo)
        pale1 = sns.cubehelix_palette(11, dark=1, light = 0)
        pale2 = sns.color_palette("plasma_r",11) #color invertido
        pale3 = sns.color_palette("viridis",11)
        pale3h = sns.cubehelix_palette(11, reverse = True)#[2:]
        pale3r = sns.color_palette("magma",11)
        pale4 = sns.color_palette("magma_r",13)[:-2]
        pale5 = sns.cubehelix_palette(rot=-0.3,n_colors = 11)  #color invertido
        pale5r = sns.color_palette('cividis_r', n_colors = 12)[:-1]
        pale_bi = sns.color_palette("coolwarm", 11)

        # ---Parameter values for CoPlot---
        # # Precision metrics
        # rand_values = [0, 0.5, 0.5]
        # shape_list = ['hm', 'hm', 'hm']
        # palettes_list = [pale3, pale3h, pale3r]
        # rotation_list = [0, 0, 0]
        # switch_list = [False, False, False]

        # # memory and time for methods without borda
        # rand_values = [0, 0]
        # shape_list = ['hm', 'hm']
        # palettes_list = [pale5r, pale5]
        # rotation_list = [45, 45]
        # switch_list = [True, True]

        # motifs
        rand_values = [1, 1, 1]
        shape_list = ['arrow', 'arrow', 'arrow']
        palettes_list = [pale1, pale1, pale1]
        rotation_list = [45, 45, 45]
        switch_list = [False, False, False]


        COplot(
            overviewDF, 
            width=16, # Un poco más ancho para acomodar nombres y escalas
            height=len(overviewDF) * 0.7 + 6, 
            randValues=rand_values, # Valores de random predictor para cada métrica. Todos los que calculan el radio toman 1 como random value
            shape=shape_list, # modos: hm = heatmap, c = circle, s = square, rs = rounded square, w = wedge, b = bar, f = full square, text = value text, arrow = up/down/neutral arrows
            palettes=palettes_list,
            levels=resType,  
            rotation=rotation_list, # rotate labels
            switch=switch_list # invert text color in graphic
        )
        
        output_name = f"{opts.output}/{db_prefix}_grouped_overview_{'_'.join(resSuffixes)}.pdf"
        plt.savefig(output_name, bbox_inches='tight')
        print(f"Éxito: Archivo guardado como {output_name}")

if __name__ == '__main__':
    main()