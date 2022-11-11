# Simulation datasets

The simulation_lognormal and simulation_normal contain the code, result and figures about datasets from *lognormal* distribution, *normal* distribution, respectively.

## simulation_lognormal

`00rlnorm_mean.R`, `01rlnorm_sd.R` and `02rlnorm_mean_sd.R` are the scripts to generate the datasets and calculate the `FPR` and `FNR` for the 100 samples and 1000 features. `00rlnorm_mean_ss.R`, `01rlnorm_sd_ss.R` and `02rlnorm_mean_sd_ss.R` are the code to generate the datasets and calculate the results of `FPT` and `FNT` for 10 samples and 1000 features. The results are saved the directory of `result`. Then, the scripts of `scripts` directory are used to visualize the results and the figures are saved the `figures` directory.

## simulation_normal

`00rnorm_mean.R`, `01rnorm_sd.R` and `02rnorm_meansd.R` are the scripts to generate the datasets and calculate the `FPR` and `FNR` for the 100 samples and 1000 features. `00rnorm_mean_smallsamples.R`, `01rnorm_sd_smallsamples.R` and `02rnorm_meansd_smallsamples.R` are the scripts to generate the datasets and calculate the `FPR` and `FNR` for the 10 samples and 1000 features. The results also are saved the `result` directory. Then the scripts of `scripts` directory are used to visualize the results and the figures are saved the `figures` directory.

# real_data

The `run_multi_daa.prop.R` was used to perform the differential abundance analysis for each datasets (The input is the files of `dataset` and return the results into the `result` dir). `plot_num.R` and  plot_upset.prop.R were used to visualize the results.

# Requirements

```
MicrobiotaProcess,
metagenomeSeq,
ANCOMBC,
tidyverse,
GUniFrac,
MicrobiomeStat,
matrixStats

LEfSe
```

## requirements R packages

```
Biocpkgs <- c('MicrobiotaProcess', 'metagenomeSeq', 'ANCOMBC', 'tidyverse', 'GUniFrac', 'MicrobiomeStat', 'matrixStats')
for (i in Biocpkgs){
    if (!requireNamespace(i, quietly = TRUE)){
        BiocManager::install(i)
    }
}
```

## LEfSe

The download of `LEfSe` refer to the `https://bitbucket.org/biobakery/biobakery/downloads/`. The `format_input.py` and `run_lefse.py` should be added to your PATH environment. The test can refer the following code.

```
format_input.py --help
run_lefse.py --help
```
