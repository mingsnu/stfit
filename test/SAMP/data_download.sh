#!/usr/bin/env bash
# 
# Filename: data_download.sh
# Author: Weicheng Zhu
# Contact: mingsnu@gmail.com
# Created: Sun May 13 21:21:07 2018 (-0500)
# Version: 
# Last-Updated: Mon May 14 00:29:07 2018 (-0500)
#           By: Weicheng Zhu

# Commentary: 
#
#for i in 2017.01.0{1..9} 2017.01.{10..31} 2017.02.0{1..9} 2017.02.{10..29} 2017.03.0{1..9} 2017.03.{10..31};
##for i in 2017.{0{5,7,8},10,12}.{1..31};do  echo $i; done;

## for i in 2017.0{5}.0{1..9} 2017.01.{10..31} 2017.02.0{1..9} 2017.02.{10..29} 2017.03.0{1..9} 2017.03.{10..31};
for i in 2017.{0{5,7,8},10,12}.{1..31} 2017.{0{4,6,9}, 11}.{1..30};
do  echo $i;
    wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off https://n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP.004/$i/;
done;

# wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies --no-check-certificate --auth-no-challenge=on -r --reject "index.html*" -np -e robots=off https://n5eil01u.ecs.nsidc.org/SMAP/SPL3SMP.004/2017.01.03/
