#THEORY OF DATA SYSTEMS - R Reader for SMAP soil moisture data
#Date: April 5, 2018
#What: This R script was built to read data from downloaded from the SMAP website.

#LIBRARIES ==============================================================
#source("https://bioconductor.org/biocLite.R")
#biocLite("rhdf5")
#install.packages("smapr")
#install.packages("xml2")
#n
library(smapr)
library(xml2)
library(raster)
library(rvest)
library(httr)
library(rappdirs)
library(rhdf5)
library(rgdal)
library(MBA)
library(fields)

#USER INPUT ==========================================================
#enter your username and password for "EOSDIS Earthdata Login"
Sys.setenv(ed_un = "mingsnu", ed_pw = "X36ZNkOO")

# searching across a date range
start_date <- as.Date("2017-01-01")
end_date <- as.Date("2017-01-05")

#DATA DOWNLOAD =======================================================
date_sequence <- seq(start_date, end_date, by = 1)
date_sequence_len<-length(date_sequence)
#where the file will be downloaded
direct<-getwd()
soil_moisture<-array(NA,dim=c(391384,date_sequence_len))

for(d in 1:date_sequence_len){
  files<-find_smap(id = "SPL3SMP", dates = date_sequence[d], version = 4)
  downloads <- download_smap(files, directory = direct, overwrite = TRUE)
  extract<-extract_smap(downloads,"/Soil_Moisture_Retrieval_Data_AM/soil_moisture")
  if(d==1){coord<-coordinates(extract)}
  soil_moisture[,d]<-values(extract)
}

#SAVE TXT FILE ======================================================
write.table(soil_moisture,"mySMAPdata.txt", sep="\t")
write.table(coord,"mySMAPcoord.txt", sep="\t")

