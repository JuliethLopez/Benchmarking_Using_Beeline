# Include packages

print("RUNNING PIDC!!!!!")

using NetworkInference
using LightGraphs
using Base.Threads

println("RUNNING PIDC correctly!!!!!")
println("Running with ", Threads.nthreads(), " threads.") # Verify multi-threading

algorithm = PIDCNetworkInference()

dataset_name = string(ARGS[1])

@time genes = get_nodes(dataset_name);

@time network = InferredNetwork(algorithm, genes);

write_network_file(string(ARGS[2]), network);

