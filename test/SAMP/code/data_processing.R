library(h5)

## reading one file and visualize ===================
file <- h5file("../n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP.004/2017.01.01/SMAP_L3_SM_P_20170101_R14010_001.h5")
list.datasets(file)
long = file["/Soil_Moisture_Retrieval_Data_AM/longitude"][]
long[long==-9999]=NA
lat = file["/Soil_Moisture_Retrieval_Data_AM/latitude"][]
lat[lat==-9999]=NA
mat =file["/Soil_Moisture_Retrieval_Data_AM/soil_moisture"][]
mat[mat == -9999] = NA
h5close(file)
yrange = range(lat,na.rm = TRUE)
xrange = range(long,na.rm = TRUE)
levelplot(raster(mat), xmn=xrange[1], xmx=xrange[2], ymn=yrange[1], ymx=yrange[2])

## reading one month of data and visualize ===============
if(length(list.files("../n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP.004/")) !=365)
  stop("number of day not correct!")
smapdat = matrix(NA,365, 391384)
month = sprintf("%02d", 1:12)
doy = sprintf("%02d", 1:31)
k=0
for(j in 1:length(month)){
  for(i in 1:length(doy)){
    
    fname =  list.files(paste0("../n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP.004/2017.",month[j],".", doy[i], "/"), 
                        pattern = "*.h5$", full.names = TRUE)
    # fname = paste0("../n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP.004/2017.",month[j],".", doy[i], 
    #                "/SMAP_L3_SM_P_2017", month[j], doy[i], "_R14010_001.h5")
    if(length(fname) > 1){
      warning("month ", month[j], "doy", doy[i], " has more than one h5 file.")
      fname = fname[1]
    }
    if(length(fname) == 1){
      k = k + 1
      if(file.exists(fname)){
        cat(k, "\n")
        cat("Loading", fname, "\n")
        file <- h5file(paste0("../n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP.004/2017.01.", doy[i], "/SMAP_L3_SM_P_201701",
                              doy[i], "_R14010_001.h5"))
        #list.datasets(file)
        mat =file["/Soil_Moisture_Retrieval_Data_AM/soil_moisture"][]
        mat[mat == -9999] = NA
        h5close(file)
        smapdat[k,] = c(t(mat)) ## stacking in row order
      }
    }
  }
}
saveRDS(smapdat, "../data/smapdat_rowstack.rds")


