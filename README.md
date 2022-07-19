# Introgression Detection - simulating data with different parameters and analyse it using different tools

## Instructions for replicating the pipeline

This GitHub repository contains [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines for replicating the analysis of my Master Thesis "A performance comparison of tools to detect introgressed fragments". These pipelines were tested on the LiSC (Life Science Compute Cluster of the University of Vienna) and also on a Linux operating system (Ubuntu 20.04.2 LTS).

First of all, it is necessary to install [Anaconda](https://www.anaconda.com/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html). After that, the following commands can be used to create a virtual environment: 

	conda config --set safety_checks disabled
	conda config --set channel_priority strict
	conda env create -f environment.yml
	conda activate introgression-kruisz

# Download the required tools

To download the required tools that cannot be installed through `conda`, the following commands can be used:

	mkdir ext && cd ext

	# Download SPrime and pipeline
	mkdir SPrime && cd SPrime
	wget https://faculty.washington.edu/browning/sprime.jar
	git clone https://github.com/YingZhou001/sprimepipeline
	chmod a+x sprimepipeline/pub.pipeline.pbs/tools/map_arch_genome/map_arch
	cd ..

	# Download SkovHMM
	mkdir SkovHMM && cd SkovHMM
	git clone https://github.com/LauritsSkov/Introgression-detection
	cd ..


It is also necessary to download the file `ms.tar.gz` from [Hudson Lab](http://home.uchicago.edu/~rhudson1/source/mksamples.html). Next decompress it under the `ext` folder and compile it with the following commands:

	cd msdir
	${CONDA_PREFIX}/bin/gcc -o ms ms.c streec.c rand1.c -lm
	cd ../..

The tool sstar is already installed as a package using the previously created environment.

# Run the pipelines

After installing the tools, it is recommended to test the pipeliens locally beforehand with a dry-run. This can be done, for example, with the following command:

	snakemake -s workflows/1src/simulation_proportion.snake -np

It is also recommended to run the individual pipelines one after the other. For the first part, the data simulations with the different parameters, the following commands can be used (`-c` specifies the number of threads):

	snakemake -s workflows/1src/simulation_proportion.snake -c 1
	snakemake -s workflows/1src/simulation_introtime.snake -c 1
	snakemake -s workflows/1src/simulation_divtime.snake -c 1

For the second part, going through the individual tools SPrime, SkovHMM and sstar, the following commands (in any order) can be used:

	snakemake -s workflows/1src/sprime_proportion.snake -c 1
	snakemake -s workflows/1src/sprime_introtime.snake -c 1
	snakemake -s workflows/1src/sprime_divtime.snake -c 1

	snakemake -s workflows/1src/sstar_proportion.snake -c 1
	snakemake -s workflows/1src/sstar_introtime.snake -c 1
	snakemake -s workflows/1src/sstar_divtime.snake -c 1

For the SkovHMM pipeline it is necessary to create a custom environment, because it requires `Python2`. Therefore '--use-conda' must be added to the command:

	snakemake -s workflows/1src/skovhmm_proportion.snake --use-conda -c 1
	snakemake -s workflows/1src/skovhmm_introtime.snake --use-conda -c 1
	snakemake -s workflows/1src/skovhmm_divtime.snake --use-conda -c 1

If the pipelines are also to be sent as jobs to a cluster, then it is recommended to perform the following steps: 
Create a  [profiles](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) depending on the cluster. 
The profile used in case for the `SLURM Workload Manager` is located in `config/slurm/config.yaml`.
Creating a folder named `logs_slurm` and then submitting the jobs with, for example, the following commands (`-j` specifies the number of threads (cluster)):

	mkdir logs_slurm
	snakemake -s workflows/1src/simulation_proportion.snake --profile config/slurm -j 200

It is also recommended, since the pipelines run for a while, to use the nohup command, which suppresses the HUP-signal and thus allows a program to continue running even if you have logged off the system. This works with the following example command:

	nohup snakemake -s workflows/1src/simulation_proportion.snake --profile config/slurm -j 200 &

There is also a small script to create individual plots. This pipeline can be executed with the following commands:

	snakemake -s workflows/plots/plots_proportion.snake
	snakemake -s workflows/plots/plots_introtime.snake
	snakemake -s workflows/plots/plots_divtime.snake
