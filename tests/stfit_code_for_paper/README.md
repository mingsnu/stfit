This folder contains all the code to reproduce all the results shown in the paper.

1. `main1.R`: code for reproducing all graphs in both the paper and supplementary material and Table S1 and S2 in supplementary material. (Take about 1 to 2 minutes to run.)
2. `main2.R`: codes for reporducing all Tables in the paper and Table S3 to S10 in the supplementary material, which depends on all experiment results in the landsat_simulation_study folder. (Take a few seconds to run if the experiment results are ready.)
3. `landsat_simulation_study` folder contains multiple subfolders, and within each there is an R file which can be used to reproduce the experiment results saved in the `output` folder. It usually takes several hours to run each experiment (except for the kriging ones, which usually only takes less than 15 minutes.), so to save time reader can just load the pre-saved results in the `output` folders. To run the code yourself: It is recommended to set the working direcotry to the path where R file exists, and all the experiment resutls will be saved in the `output` folder saved in 'rds` format. 


