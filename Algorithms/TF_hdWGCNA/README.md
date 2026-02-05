*This README.md file was generated on 1/19/2023 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# TF_hdWGCNA: Transcription factor - higly simensional weighted gene coexpression network analysis

This is the instruction on how to integrate TF_hdWGCNA ([[Paper](https://www.cell.com/cell-reports-methods/fulltext/S2667-2375(23)00127-3)] [[GitHub](https://github.com/smorabit/hdWGCNA)]) to BEELINE. Please follow the following steps:

1. **Create TF_hdWGCNA folder:** Create a folder called TF_hdWGCNA under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Create runTF_hdWGCNA.R script:** In the TF_hdWGCNA folder, create an R script runTF_hdWGCNA.R to learn graph from target datasets.

3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. The Dockerfile will run the TF_hdWGCNA algorithm within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for TF_hdWGCNA.

        cd $BASEDIR/Algorithms/TF_hdWGCNA/
        docker build -q -t TF_hdWGCNA:base .
        if ([[ "$(docker images -q TF_hdWGCNA:base 2> /dev/null)" != "" ]]); then
            echo "Docker container for TF_hdWGCNA is built and tagged as TF_hdWGCNA:base"
        else
            echo "Oops! Unable to build Docker container for TF_hdWGCNA"
        fi

5. **Create TF_hdWGCNARunner.py script:** After buliding the Docker image, create a Python script called TF_hdWGCNARunner.py in Beeline/BLRun folder to setup BLRun objects to read inputs and run TF_hdWGCNA inside the Docker image, and also parse the output for evaluation. Specifically, the TF_hdWGCNARunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the input data file (i.e., expression data), and processes them into the format required by the TF_hdWGCNA. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data file (i.e., expression data). It also specifies where the outputs are written. The docker container run TF_hdWGCNA when the parameters are passed. 
   - ``parseOutput()`` : This function reads the TF_hdWGCNA-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add TF_hdWGCNA to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to TF_hdWGCNA. 

    - add "import BLRun.TF_hdWGCNARunner as TF_hdWGCNA"
    - add "'TF_hdWGCNA':TF_hdWGCNA.generateInputs" to InputMapper
    - add "'TF_hdWGCNA':TF_hdWGCNA.run" to AlgorithmMapper
    - add "'TF_hdWGCNA':TF_hdWGCNA.parseOutput" to OutputParser

7. **Add TF_hdWGCNA to config.yaml:** The final step is to add the new algorithm TF_hdWGCNA and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.

        - name: "TF_hdWGCNA"
          params: 
              should_run: [True]