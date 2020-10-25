This folder contains all the code to reproduce all the results shown in the paper.

**System Info**

    macOS Catalina
    Version 10.15.6
    MacBook Pro (Retina, 15-inch, Mid 2015)
    Processor 2.5 GHz Quad-Core Intel Core i7
    Memory 16 GB 1600 MHz DDR3

** R Session Info**

    R version 4.0.0 (2020-04-24)
    Platform: x86_64-apple-darwin17.0 (64-bit)
    Running under: macOS Catalina 10.15.6
    
    Matrix products: default
    BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
    LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
    
    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    
    attached base packages:
    [1] parallel  stats     graphics  grDevices utils     datasets  methods   base     
    
    other attached packages:
     [1] RColorBrewer_1.1-2  stfit_0.99.8        rasterVis_0.47     
     [4] latticeExtra_0.6-29 lattice_0.20-41     Matrix_1.2-18      
     [7] doParallel_1.0.15   iterators_1.0.12    foreach_1.5.0      
    [10] dplyr_0.8.5         feather_0.3.5       raster_3.1-5       
    [13] sp_1.4-1           
    
    loaded via a namespace (and not attached):
     [1] Rcpp_1.0.4.6      pillar_1.4.4      compiler_4.0.0    tools_4.0.0      
     [5] lifecycle_0.2.0   tibble_3.0.1      viridisLite_0.3.0 pkgconfig_2.0.3  
     [9] png_0.1-7         rlang_0.4.6       cli_2.0.2         rstudioapi_0.11  
    [13] hexbin_1.28.1     vctrs_0.2.4       hms_0.5.3         grid_4.0.0       
    [17] tidyselect_1.0.0  glue_1.4.0        R6_2.4.1          jpeg_0.1-8.1     
    [21] fansi_0.4.1       purrr_0.3.4       magrittr_1.5      codetools_0.2-16 
    [25] ellipsis_0.3.0    assertthat_0.2.1  crayon_1.3.4      zoo_1.8-8        
   

**Instruction**
You can reproduce all the results presented in the paper by running `main1.R` and `main2.R`. The results will be saved in the `output` folder. (Note: please always set the current working direcotry to the code location for all the R script presented in this folder.)

1. `main1.R`: code for reproducing all graphs in both the paper and supplementary material and Table S1 and S2 in supplementary material. (Take about 1 to 2 minutes to run.)
2. `main2.R`: codes for reporducing all Tables in the paper and Table S3 to S10 in the supplementary material, which depends on all experiment results in the `landsat_simulation_study` folder. (Take a few seconds to run if the experiment results are ready.)
3. `landsat_simulation_study` folder contains multiple subfolders, and within each there is an R file which can be used to reproduce the experiment results saved in the `output` folder. It usually takes several hours to run each experiment (except for the kriging ones, which usually only takes less than 15 minutes.), so to save time reader can just load the pre-saved results in the `output` folders. To run the code yourself: It is recommended to set the working direcotry to the path where R file exists, and all the experiment resutls will be saved in the `output` folder saved in 'rds` format. 


