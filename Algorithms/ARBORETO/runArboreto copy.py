# prueba_hilos.py
import dask
dask.config.set({"distributed.comm.serialization": ["pickle"]})
from dask.distributed import Client, LocalCluster
import numpy as np

cluster = LocalCluster(n_workers=2, threads_per_worker=2, processes=False)
client = Client(cluster)
print("Cluster threads:", cluster)
f = client.submit(lambda x: x.sum(), np.arange(1000000))
print("Result:", f.result())
client.close()
cluster.close()
