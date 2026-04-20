import os
from pathlib import Path

def getMemory(evalObject, dataset):
    """
    Extract peak memory usage (Maximum Resident Set Size) for each algorithm
    in a dataset, using the same structure as getTime().
    
    Returns:
        dict: {algorithm_name: memory_in_kB}
    """

    outDir = str(evalObject.output_settings.base_dir) + \
             str(evalObject.input_settings.datadir).split("inputs")[1] + \
             "/" + dataset["name"] + "/"

    algo_memory = dict()
    algos = evalObject.input_settings.algorithms

    for algo in algos:
        path = outDir+algo[0]+"/time.txt"

        # Case 1: single file
        if Path(path).exists():
            mem = parse_memory_file(path)

        # Case 2: multiple trajectory files
        else:
            import pandas as pd

            PTFile = str(evalObject.input_settings.datadir.joinpath(
                dataset['name'] + '/' + dataset['cellData']
            ))

            PTData = pd.read_csv(PTFile, header=0, index_col=0)
            colNames = PTData.columns

            mem = -1

            for idx in range(len(colNames)):
                path_i = outDir+algo[0]+"/time"+str(idx)+".txt"
                mem_i = parse_memory_file(path_i)
                mem = max(mem, mem_i)

        if mem == -1:
            print("Skipping memory computation for", algo[0], "on dataset", dataset["name"])
            continue

        algo_memory[algo[0]] = mem

    return algo_memory

def parse_memory_file(path):
    """
    Extract Maximum Resident Set Size (kB) from time.txt
    """

    try:
        with open(path, "r") as f:
            for line in f:
                line = line.strip()

                if line.startswith("Maximum resident set size"):
                    return float(line.split(":")[-1].strip().split()[0])

        print("Memory not found in:", path)
        return -1

    except FileNotFoundError:
        print("File not found:", path)
        return -1

    except ValueError:
        print("Parsing error in:", path)
        return -1