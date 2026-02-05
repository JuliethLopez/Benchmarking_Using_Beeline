*This README.md file was generated on 1/19/2023 by Yiqi Su (yiqisu@vt.edu)*
<!-- remove all comments (like this) before final save  -->

# INFERELATOR: Transcription factor - higly simensional weighted gene coexpression network analysis

This is the instruction on how to integrate INFERELATOR ([[Paper](https://doi.org/10.1093/bioinformatics/btac117)] [[GitHub](https://github.com/flatironinstitute/inferelator)]) to BEELINE. Please follow the following steps:

1. **Create INFERELATOR folder:** Create a folder called INFERELATOR under Beeline/Algorithms for the new method to ensure easy set-up and portability and avoid conflicting libraries/software versions that may arise from the GRN algorithm implmentations.

2. **Create runINFERELATOR.R script:** In the INFERELATOR folder, create an R script runINFERELATOR.R to learn graph from target datasets.

3. **Create a Dockerfile:** Create a "Dockerfile" that contains necessary software specifications and commands listed in a specific order from top to bottom. The Dockerfile will run the INFERELATOR algorithm within the Docker container.

4. **Add the Dockerfile to initialize.sh script:** Once the Dockerfile is ready, add the following lines to 'initialize.sh' script to create Docker image for INFERELATOR.

        cd $BASEDIR/Algorithms/INFERELATOR/
        docker build -q -t INFERELATOR:base .
        if ([[ "$(docker images -q INFERELATOR:base 2> /dev/null)" != "" ]]); then
            echo "Docker container for INFERELATOR is built and tagged as INFERELATOR:base"
        else
            echo "Oops! Unable to build Docker container for INFERELATOR"
        fi

5. **Create INFERELATORRunner.py script:** After buliding the Docker image, create a Python script called INFERELATORRunner.py in Beeline/BLRun folder to setup BLRun objects to read inputs and run INFERELATOR inside the Docker image, and also parse the output for evaluation. Specifically, the INFERELATORRunner.py script contains three functions:

   - ``generateInputs()`` : This function reads the input data file (i.e., expression data), and processes them into the format required by the INFERELATOR. 
   - ``run()`` : This function constructs a "docker run" system command with parameters including the path of input data file (i.e., expression data). It also specifies where the outputs are written. The docker container run INFERELATOR when the parameters are passed. 
   - ``parseOutput()`` : This function reads the INFERELATOR-specific output (i.e., outFile.txt) and formats it into a ranked edgelist comma-separated file (i.e., rankedEdges.csv) with columns Gene1, Gene2, and EdgeWeight. The Gene1 column should contain regulators, the Gene2 column the targets, and EdgeWeight column the absolute value of the weight predicted for edge (regulator,target). The ranked edgelist file will be subsequently used by BLEval object. 

6. **Add INFERELATOR to runner.py:** Next, update runner.py script in Beeline/BLRun folder by adding information related to INFERELATOR. 

    - add "import BLRun.INFERELATORRunner as INFERELATOR"
    - add "'INFERELATOR':INFERELATOR.generateInputs" to InputMapper
    - add "'INFERELATOR':INFERELATOR.run" to AlgorithmMapper
    - add "'INFERELATOR':INFERELATOR.parseOutput" to OutputParser

7. **Add INFERELATOR to config.yaml:** The final step is to add the new algorithm INFERELATOR and any necessary parameters to the config.yaml located in Beeline/config-files folder. Note that currently BEELINE can only handle one parameter set at a time eventhough multiple parameters can be passed onto the single parameter object.

        - name: "INFERELATOR"
          params: 
              should_run: [True]