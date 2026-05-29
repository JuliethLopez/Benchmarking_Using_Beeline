from optparse import OptionParser
import os
import sys
import pandas as pd
import dask
from dask.distributed import Client, LocalCluster
from arboreto.algo import grnboost2, genie3
from arboreto.utils import load_tf_names

# First approach
# dask.config.set({
#     "distributed.comm.serialization": ["pickle"]
# })

# Second approach
#dask.config.set({"distributed.scheduler.allowed-failures": 10})

def parseArgs(args):
    parser = OptionParser()

    parser.add_option('', '--algo', type = 'str',
                      help='Algorithm to run. Can either by GENIE3 or GRNBoost2')

    parser.add_option('', '--inFile', type='str',
                      help='Path to input tab-separated expression SamplesxGenes file')

    parser.add_option('', '--outFile', type = 'str',
                      help='File where the output network is stored')

    (opts, args) = parser.parse_args(args)

    return opts, args

def main(args):
    opts, args = parseArgs(args)
    inDF = pd.read_csv(opts.inFile, sep = '\t', index_col = 0, header = 0)

    # Configurar el clúster Dask local o conectar a uno existente (aqui se debe poner expose en docker/singularity run)
    cluster = None
    client = None

    scheduler = os.environ.get("DASK_SCHEDULER", None)

    if scheduler is not None:
        try:
            client = Client(scheduler, timeout="5s")
            print("Connected to external cluster:", client)

        except Exception as e1:
            print("Could not connect to external cluster:", e1)
            client = None

    if client is None:
        try:
            # Dfine number of cpu
            n_cpus = int(os.environ.get("LSB_DJOB_NUMPROC", 2))
            print(f"Detected CPUs: {n_cpus}")

            # do not use local for singularity/HPC
            # cluster = LocalCluster(
            #     n_workers=n_cpus,
            #     threads_per_worker=1,
            #     processes=False,
            #     memory_limit="auto",
            #     silence_logs=True,
            #     dashboard_address=None
            # )
            # client = Client(cluster)
            client = Client(n_workers=n_cpus,
                            threads_per_worker=1,
                            processes=False)
            print("Started LocalCluster:", cluster)

        except Exception as e2:
            print("Could not start LocalCluster:", e2)

            # -----------------------------
            # 3. Last resort: in-process execution
            # -----------------------------
            client = Client(processes=False)
            print("Started simple in-process Client:", client)


    if opts.algo == 'GENIE3':
        network = genie3(inDF.to_numpy(), client_or_address = client, gene_names = inDF.columns)
        #network = genie3(inDF.to_numpy(), gene_names = inDF.columns)
        network.to_csv(opts.outFile, index = False, sep = '\t')

    elif opts.algo == 'GRNBoost2':
        network = grnboost2(inDF.to_numpy(), client_or_address = client, gene_names = inDF.columns)
        network.to_csv(opts.outFile, index = False, sep = '\t')

    else:
        print("Wrong algorithm name. Should either be GENIE3 or GRNBoost2.")
                        
if __name__ == "__main__":
    main(sys.argv)
