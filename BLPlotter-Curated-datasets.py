#!/usr/bin/env python
# coding: utf-8

import os
import sys
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
    '''
    :return: an argparse ArgumentParser object for parsing command
        line parameters
    '''
    parser = argparse.ArgumentParser(
        description='Generate plots from evaluation results.')

    parser.add_argument('-c','--config', default='config.yaml',
        help="Comma delimited configuration file(s) containing list of datasets "
              "algorithms and output specifications.\n")

    parser.add_argument('-a', '--auprc', action="store_true", default=False,
        help="Generate box plot of AUPRC values for evaluated algorithms.\n")
    
    parser.add_argument('-r', '--auroc', action="store_true", default=False,
        help="Generate box plot of AUROC values for evaluated algorithms.\n")
    
    parser.add_argument('-e', '--epr', action="store_true", default=False,
        help="Generate box plot of early precision values for evaluated algorithms.\n")

    parser.add_argument('-v', '--overview', action="store_true", default=False,
      help="Generate plot of AUPRC and early precision ratios relative to a random predictor.\n")
    
    parser.add_argument('-o', '--output', default='.',
      help="Output directory for generated plots.\n")

    return parser


def parse_arguments():
    '''
    Initialize a parser and use it to parse the command line arguments
    :return: parsed dictionary of command line arguments
    '''
    parser = get_parser()
    opts = parser.parse_args()

    return opts


def boxplot(opts, evalConfigs, datasets, randValue, resTypeFile, resTypeName):
    
    # Set plot variables
    plt.rcParams.update({'font.size': 14})
    
    DFs = []
    for i, dataset in enumerate(datasets):
        evalConfig = evalConfigs[i]
        
        # Read output file containing evaluated values
        DF = pd.read_csv(str(evalConfig.output_settings.base_dir) + '/' \
                           + str(evalConfig.input_settings.datadir).split("inputs")[1] + '/' \
                           + str(evalConfig.output_settings.output_prefix) \
                           + '-' + resTypeFile + '.csv', header = 0, index_col = 0)
        
        DF = DF.T

        #  Add row index as column
        DF.index.name = 'dataset'
        DF.reset_index(inplace = True)  
        
        DFs.append(DF)  

    
    # Initialize plot window
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
            modifiedDF['dataset'].str.contains('-1000-')
            #modifiedDF['dataset'].str.contains('-2000-')
        ]
        outputs = ['200', '1000']
        existing_cats = ['200', '1000', '2000']
        modifiedDF['category'] = pd.Series(np.select(conditions, outputs, '2000')) #'0'
        modifiedDF['category'] = modifiedDF['category'].astype('category').cat.set_categories(existing_cats, ordered=True)
        modifiedDF['category'] = modifiedDF['category'].astype('category').cat.reorder_categories(['200', '1000', '2000'])
              
        if len(datasets) > 1:
            subax = axes[i]
        else:
            subax= axes
        
        # Indicate random predictors value
        ax = sns.lineplot(y = randValue, 
                          x = range(-1,len(DF.columns)+1), ax = subax)
        ax.lines[0].set_linestyle("--")
        ax.lines[0].set_color("gray")
        
        # Plot the evaluated values as a box plot
        ax = sns.boxplot(y = 'value', x = 'variable', data = modifiedDF, fliersize = 0,
                        hue = "category",
                        palette = sns.color_palette("Set1"),
                        ax = subax)
        # Rotate x-axis labels
        subax.set_xticklabels(
            subax.get_xticklabels(),
            rotation=45,
            ha='right'  # aligns labels with ticks more accurately
        )
        # Get the legend from just the box plot
        handles, labels = ax.get_legend_handles_labels()
        sns.swarmplot(y ='value', x = 'variable', data = modifiedDF,
                        hue = "category",
                        alpha = 0.5, palette = ['k'],
                        dodge = True,
                        ax = subax)
        # Remove the old legend
        ax.legend_.remove()
        # Add just the handles/labels from the box plot back
        ax.legend(handles, labels,
                  title = 'Number of cells', bbox_to_anchor=(.5, 1), loc = 'upper center', ncol = len(modifiedDF['category'].cat.categories))
      
        ax.set_ylim([0.0,1])
        ax.set_ylabel(resTypeName, fontsize = 18)
        ax.set_title(resTypeName +' values for ' + dataset +' Network', fontsize = 18)
        plt.tight_layout()
    
    subax.set_xlabel('Algorithm', fontsize = 18)

               
    file = opts.output + '/' \
            + '-'.join([str(c.output_settings.output_prefix) for c in evalConfigs]) \
            + '-boxplot-' + resTypeFile + '.pdf'
    print("Boxplot saved to " + file)
    plt.savefig(file, dpi = 300)

def COplot(inputDF, width = 12, height = 7, randValues = [], shape = [], 
            palettes = [], levels = [], rotation = [], switch = [], metric_type = 'score'):   #    metric_type = 'ratio' or 'score'

    levls = levels
    rowNames = inputDF.index
    maxRows = len(inputDF.index)
    maxCols = len(inputDF.columns)
    pad = 2
    fSize = (width,height)

    f = plt.figure(figsize=fSize)

    ax_main = plt.subplot2grid((2*maxRows + 3, maxCols), (0, 0), colspan = maxCols, rowspan = 2*maxRows + 2)
    ax_main.set_yticks(np.arange(0,maxRows+pad))
    ax_main.set_xticks(np.arange(0,maxCols+pad))

    Indices = [""] + list(rowNames)  +  [""]

    ax_main.set_yticklabels(Indices, fontsize=30, fontweight='bold')
    ax_main.spines['top'].set_visible(False)
    ax_main.spines['right'].set_visible(False)
    ax_main.spines['bottom'].set_visible(False)
    ax_main.spines['left'].set_visible(False)

    # Alternate gray and white row background colors
    alt = False
    for rowIdx in range(len(rowNames)):
        if alt:
            ptch = patches.Rectangle((0,rowIdx+0.5),
                               width = maxCols + 1,
                               height = 1,  
                               edgecolor = (1,1,1),
                               facecolor = (0.9,0.9,0.9),)

        else:
            ptch = patches.Rectangle((0,rowIdx+0.5),
                   width = maxCols,
                   height = 1,  
                   edgecolor = (1,1,1),
                   facecolor = (1,1,1),)
        alt = not(alt)
        # ax.add_artist(ptch)
        ax_main.add_artist(ptch)

    colCnt = 0
    levl_col_widths = [len(inputDF[lvl].columns) for lvl in levls]

    for levlIdx in range(len(levls)):
        randValue = randValues[levlIdx]
        print(f"Plotting level: {levls[levlIdx]} with random predictor value: {randValue}")
        colStart = colCnt
        colNames = inputDF[levels[levlIdx]].columns
        
        # Plot metric column headers
        header_x = colStart + (len(colNames) / 2.0) + 0.5
        plt.text(header_x, maxRows + 2, 
                        levls[levlIdx], fontsize=30, rotation=rotation[levlIdx],
                        fontweight='bold',
                        ha="center", va="center",
                        bbox=dict(boxstyle="round",
                        ec=(1,1,1), fc=(1,1,1)))

        for colIdx in range(len(colNames)):
                colCnt += 1

                # Plot dataset column headers
                plt.text(colStart + colIdx + 1, maxRows + 1, colNames[colIdx],
                         fontsize=25, rotation=rotation[levlIdx],
                         ha="center", va="center",
                         bbox=dict(boxstyle="round",
                         ec=(1,1,1,0),
                         fc=(1,1,1,0)))

                for rowIdx in range(len(rowNames)):

                    value = inputDF.loc[rowNames[rowIdx],levls[levlIdx]][colNames[colIdx]]
                    
                    if shape[levlIdx] == 'arrow':                    

                        valList = ['\u2197','\u2198','\u003d']

                        # Configure arrow direction and color for arrow shape cells
                        if value > 1:
                            txt = valList[0]
                            col = '#4CAF50'
                            fSize = 24

                        elif value < 1:
                            txt = valList[1]
                            col = '#E53935'
                            fSize = 24

                        elif value == 1:
                            txt = valList[2]
                            col = '#2E4053'
                            fSize = 24

                        else:
                            # Missing values?
                            txt = '-'
                            col = '#2E4053'
                            fSize = 24

                        plt.text(colStart+colIdx+1, rowIdx+1, 
                            txt, fontsize= fSize, rotation=0,
                            ha="center", va="center", color = col,
                            bbox=dict(boxstyle="round", 
                            ec=(1,1,1,0), fc=(1,1,1,0)))

                    if shape[levlIdx] == 'text': 
                        fSize = 23

                        if value < 1:
                            txt = round(value,2)
                        else:
                            txt = int(value)

                        # Plot value as text in text cells
                        col = '#2E4053'
                        fSize = 11
                        plt.text(colStart+colIdx+1, rowIdx+1, 
                            txt, fontsize= fSize, rotation=0,
                            ha="center", va="center", color = col,
                            bbox=dict(boxstyle="round", 
                            ec=(1,1,1,0), fc=(1,1,1,0)))              
                    
                    if shape[levlIdx] not in ['text','arrow']:
                        txt = value
                        
                        if np.isnan(value):
                            # Plot '-' text and gray cell color for a missing value
                            txt = '-'
                            col = '#2E4053'
                            fSize = 24

                            plt.text(colStart+colIdx+1, rowIdx+1, 
                            txt, fontsize= fSize, rotation=0,
                            ha="center", va="center", color = col,
                            bbox=dict(boxstyle="round", 
                            ec=(1,1,1,0), fc=(1,1,1,0)))
                            continue
                        elif value < randValue:
                            # Plot a gray cell color for a values worse than random 
                            col = '#2E4053'
                        else:
                            # Plot a cell color based on the range of cell values and column palette
                            value = round(value,3)
                            
                            rangesD = (round(inputDF[levls[levlIdx]][colNames[colIdx]].min(),3), 
                                      round(inputDF[levls[levlIdx]][colNames[colIdx]].max(),3))

                            value = (value-rangesD[0])/(rangesD[1]-rangesD[0])
                            col = palettes[levlIdx][int(np.floor((value)*10))]


                        if shape[levlIdx] == 'c':
                            # Plot a circle sized by the cell value
                            circle1=patches.Circle((colStart+colIdx+1,rowIdx+1),
                                               radius = np.sqrt(value)/2.5,
                                               facecolor=col,
                                               edgecolor = 'k',)

                        elif shape[levlIdx] == 's':
                            # Plot a fixed size square with a black border
                            size = value*0.8
                            circle1=patches.Rectangle((colStart+colIdx+1-(size/2),rowIdx+1-(size)/2),
                                       width = size,
                                       height = size,    
                                       facecolor=col,
                                       edgecolor = 'k',)
                            
                        
                        elif shape[levlIdx] == 'hm':
                            # Plot a fixed size square with no border
                            size = 1
                            circle1=patches.Rectangle((colStart+colIdx+1-(size/2),rowIdx+1-(size)/2),
                                       width = size-0.1,
                                       height = size,    
                                       facecolor = col,
                                       edgecolor = col,)
                            

                        elif shape[levlIdx] == 'rs':
                            # Plot a rounded square sized by the cell value
                            size = size*0.8
                            if size <= 0.15:
                                size = 0.15
                            boxPad = 0.075
                            newVal = size - boxPad*2
                            circle1=patches.FancyBboxPatch((colStart+colIdx+1-(newVal/2),rowIdx+1-(newVal)/2),
                                       width = newVal,
                                       height = newVal,    
                                       facecolor=col,
                                       edgecolor = 'k',
                                       boxstyle=patches.BoxStyle("Round", pad=boxPad))
                                
                        elif shape[levlIdx] == 'w':
                            # Plot a sector of a circle sized by the cell value
                            circle1=patches.Wedge((colStart+colIdx+1,rowIdx+1),
                                       r = 0.4,
                                       theta1 = 0,
                                       theta2 = round(value*360,2),    
                                       facecolor = col,
                                       edgecolor = 'k',)
                        elif shape[levlIdx] == 'b':
                            # Plot a bar with height sized by the cell value
                            circle1=patches.Rectangle((colStart+colIdx+0.6,rowIdx+0.65),
                                       width = 0.75,
                                       height = value,    
                                       facecolor=col,
                                       edgecolor = 'k',)
                            
                        elif shape[levlIdx] == 'f':
                            # Plot fixed size flat color square
                            circle1=patches.Rectangle((colStart+colIdx+1,rowIdx+0.6),
                                       width = 1,
                                       height = 1,    
                                       facecolor = col,
                                       edgecolor = col,)


                        ax_main.add_artist(circle1)
                        
                        textCol = ['black','white']
                        if switch[levlIdx]:
                            textCol[0] = 'white'
                            textCol[1] = 'black'

                        # Plot the cell value as text in the maximum and minimum value cells
                        if value >= 1:
                            plt.text(colStart+colIdx+1, rowIdx+1, 
                            round(txt,1), fontsize= 22, rotation=0, fontweight='bold',
                            ha="center", va="center", color = textCol[0],
                            bbox=dict(boxstyle="round", 
                            ec=(1,1,1,0), fc=(1,1,1,0)))
                        if value <= 0:
                            plt.text(colStart+colIdx+1, rowIdx+1, 
                            round(txt,1), fontsize= 22, rotation=0, fontweight='bold',
                            ha="center", va="center", color = textCol[1],
                            bbox=dict(boxstyle="round", 
                            ec=(1,1,1,0), fc=(1,1,1,0)))
    
    ax_main.yaxis.set_ticks_position('none') 
    ax_main.xaxis.set_ticks_position('none') 
    ax_main.set_xticklabels([])

    # Track the start index position of our column layout dynamically
    leg_start_col = 0
    for levlIdx in range(len(levls)):
        current_span = levl_col_widths[levlIdx]
        
        # SKIP if arrow shape (preserving your condition)
        if shape[levlIdx] == 'arrow':
            leg_start_col += current_span
            continue
        
        # Legend grids must share the exact same global scale as the main plot (`maxCols`) 
        # instead of `len(levls)`. Each metric colorbar spans its exact underlying datasets (`colspan = current_span`)
        ax_leg = plt.subplot2grid((2*maxRows + 3, maxCols), (2*maxRows + 2, leg_start_col), colspan=current_span)
        
        ax_leg.imshow(np.arange(len(palettes[levlIdx])).reshape(1, len(palettes[levlIdx])),
              cmap=mpl.colors.ListedColormap(list(palettes[levlIdx])),
              interpolation="nearest", aspect="auto")        
        ax_leg.yaxis.set_ticks_position('none') 
        ax_leg.xaxis.set_ticks_position('none')
        ax_leg.set_yticklabels([]) 
        ax_leg.set_xticks([0.5, len(palettes[levlIdx]) - 2])
        ax_leg.set_xticklabels(['Low/Poor', 'High/Good'], fontsize=20, fontweight='bold')  

        # Progress the starting grid location for the next legend subplot block
        leg_start_col += current_span

    # Global spacing adjustment
    plt.subplots_adjust(bottom=0.15, hspace=0.6)

    return f

def main():
    opts = parse_arguments()
    
    # 1. Load the config file and extract dataset names
    with open(opts.config, 'r') as conf:
        config_dict = yaml.safe_load(conf)
        conf.seek(0)
        # Use ev.ConfigParser as per guidance to maintain config object consistency
        evalConfig = ev.ConfigParser.parse(conf)

    # CHANGE: Extract dataset names specifically from the list in your YAML
    dataset_list = [d['name'] for d in config_dict['input_settings']['datasets']]
    
    # 2. Calculate randPredictor EXACTLY as in Guidance Code
    # This ensures the 'Ratio' metrics use the correct baseline density
    randPredictor = {}
    print("Calculating network densities for random predictor baseline...")
    for ds_name in dataset_list:
        # Guidance logic: construct path to reference network
        # Assumes directory structure: [input_dir]/[dataset_dir]/[ds_name]/[trueEdges]
        true_edges_path = (
            f"{evalConfig.input_settings.datadir}/{ds_name}/"
            f"{evalConfig.input_settings.datasets[0]['trueEdges']}"
        )
        
        trueEdgesDF = pd.read_csv(true_edges_path, header=0, index_col=None)
        # Remove self-edges and duplicates as per guidance
        trueEdgesDF = trueEdgesDF[trueEdgesDF.Gene1 != trueEdgesDF.Gene2]
        trueEdgesDF.drop_duplicates(inplace=True)
        
        # Compute number of possible directed edges
        numGenes = len(np.unique(trueEdgesDF.loc[:,['Gene1','Gene2']]))
        numPossibleEdges = numGenes * (numGenes - 1)
        randPredictor[ds_name] = trueEdgesDF.shape[0] / numPossibleEdges
        print(f"Random Predictor for {ds_name}: {randPredictor[ds_name]:.4f}")

    # 3. Setup Metrics for Heatmap
    # resType = ['AUPRC Ratio', 'Stability Across Datasets', 'Early Precision Ratio']
    # resTypeFileName = ['AUPRC-undirected', 'Jaccard', 'EPr-undirected']
    resType = ['Stability Across Datasets','Stability Across Runs'] 
    resTypeFileName = ['Jaccard', 'Spearman']
    # resType = ['AUPRC Ratio', 'Stability Across Datasets', 'Early Precision Ratio']
    # resTypeFileName = ['AUPRC-directed', 'Jaccard', 'EPr-directed']
    # resType = ['Motif-FBL', 'Motif-FFL', 'Motif-MI'] #MEJOR USAR OTHER_PLOTS.R o CUrated-sin grupos
    # resTypeFileName = ['NetworkMotifs-FBL', 'NetworkMotifs-FFL', 'NetworkMotifs-MI']
    # resType = ['Memory', 'Time'] #MEJOR USAR OTHER_PLOTS.R
    # resTypeFileName = ['Memory', 'Times']
    
    # Peek at first dataset to initialize algorithm list
    first_ds = dataset_list[0]
    sample_path = (
        f"{evalConfig.output_settings.base_dir}/Curated/{first_ds}/"
        f"{first_ds}-{resTypeFileName[0]}.csv"
    )
    sample_df = pd.read_csv(sample_path, index_col=0)
    algs = sample_df.index

    # Initialize MultiIndex DataFrame
    multIndTuple = [(res, ds) for res in resType for ds in dataset_list]
    overviewDF = pd.DataFrame(index=algs, columns=pd.MultiIndex.from_tuples(multIndTuple))

    # 4. Aggregate Data with EXACT guidance calculation logic
    print("Reading result files and computing overview values...")
    for i, res_label in enumerate(resType):
        suffix = resTypeFileName[i]
        
        for ds in dataset_list:
            # CHANGE: Path logic follows the specific structure in your request
            file_path = f"{evalConfig.output_settings.base_dir}/Curated/{ds}/{ds}-{suffix}.csv"

            if os.path.exists(file_path):
                ResDF = pd.read_csv(file_path, index_col=0).T # Transpose as per guidance
                
                for alg in algs:
                    if alg in ResDF.columns:
                        # CHANGE: Use guidance calculation logic for ratios vs medians
                        if 'Ratio' in res_label:
                            overviewDF.loc[alg, (res_label, ds)] = ResDF[alg].median() / randPredictor[ds]
                        elif 'Stability' in res_label:
                            # Guidance uses the first value for stability
                            overviewDF.loc[alg, (res_label, ds)] = ResDF.iloc[0][alg]
                        elif res_label == 'Memory':
                            overviewDF.loc[alg, (res_label, ds)] = ResDF[alg].median() / 1048576 # Convert to GB
                        else:
                            overviewDF.loc[alg, (res_label, ds)] = ResDF[alg].median()
            else:
                print(f"Warning: File not found {file_path}")

    # 5. Sorting and Plotting
    overviewDF = overviewDF.dropna(how='all')
    #overviewDF = overviewDF.sort_values(by=(resType[0]), ascending=True)
    # undirected 
    overviewDF = overviewDF.reindex(['PIDC', 'GRNBOOST2', 'PPCOR', 'INFERELATOR', 'hdWGCNA', 'BORDA', 'sBORDA', 'mBORDA', 'smBORDA']).iloc[::-1]
    # # directed, motifs
    # overviewDF = overviewDF.reindex(['GRNBOOST2', 'INFERELATOR']).iloc[::-1]
    # # memory, time
    # overviewDF = overviewDF.reindex(['PIDC', 'GRNBOOST2', 'PPCOR', 'INFERELATOR', 'hdWGCNA']).iloc[::-1]
    
    print(overviewDF)

    pale1 = sns.cubehelix_palette(11, dark=1, light = 0)
    pale2 = sns.color_palette("plasma_r",11) #color invertido
    pale3 = sns.color_palette("viridis",11)
    pale3h = sns.cubehelix_palette(11, reverse = True)#[2:]
    pale3r = sns.color_palette("magma",11)
    pale4 = sns.color_palette("magma_r",13)[:-2] #color invertido
    pale5 = sns.cubehelix_palette(rot=-0.3,n_colors = 11) #color invertido
    pale5r = sns.color_palette('cividis_r', n_colors = 12)[:-1] #color invertido
    pale_bi = sns.color_palette("coolwarm", 11)
    
    # ---Parameter values for CoPlot---
    # # Precision metrics
    # rand_values = [0, 0, 0]
    # shape_list = ['hm', 'hm', 'hm']
    # palettes_list = [pale3, pale3h, pale3r]
    # rotation_list = [0, 0, 0]
    # switch_list = [False, False, False]

    # Fidelity
    rand_values = [0, 0]
    shape_list = ['hm', 'hm']
    palettes_list = [pale3, pale3h]
    rotation_list = [0, 0]
    switch_list = [False, False]

    # # memory and time for methods without borda
    # rand_values = [0, 0]
    # shape_list = ['hm', 'hm']
    # palettes_list = [pale5r, pale5]
    # rotation_list = [45, 45]
    # switch_list = [True, True]

    # # motifs
    # rand_values = [1, 1, 1]
    # shape_list = ['arrow', 'arrow', 'arrow']
    # palettes_list = [pale1, pale1, pale1]
    # rotation_list = [45, 45, 45]
    # switch_list = [False, False, False]

    fig = COplot(
        overviewDF.astype(float), 
        width=2 * len(dataset_list) * len(resType), 
        height=len(algs) + 5, 
        randValues=rand_values,
        shape=shape_list, 
        palettes=palettes_list,
        levels=resType,  
        rotation=rotation_list,
        switch=switch_list,
    )

    # Save output
    output_pdf = f"{opts.output}/Heatmap_{'_'.join(resTypeFileName)}.pdf"
    fig.savefig(output_pdf, bbox_inches='tight')
    print(f"Heatmap saved to {output_pdf}")

    fig.savefig(
    f"{opts.output}/Heatmap_{'_'.join(resTypeFileName)}.png",
    dpi=300,
    bbox_inches="tight")

if __name__ == '__main__':
  main()
