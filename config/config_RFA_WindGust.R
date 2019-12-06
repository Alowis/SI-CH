################################################################################
## file of R environment configuration
################################################################################


###########
## Paths ##
###########

dir.study = getwd()  # exemple: setwd("~/Rstudy_RFA")


if (!exists("dir.functions"))  dir.functions = paste0(dir.study,'/functions')
if (!exists("dir.input"))      dir.input     = paste0(dir.study,"/in")
if (!exists("dir.output"))     dir.output    = paste0(dir.study,"/out")
if (!exists("dir.obs"))        dir.obs       = paste0(dir.input,"/obs")
if (!exists("dir.figures"))    dir.figures   = paste0(dir.output,"/figures")
if (!exists("dir.gust"))       dir.gust      = paste0(dir.output,"/gust")

dir.create(dir.output,showWarnings = F,recursive = T)
dir.create(dir.figures,showWarnings = F,recursive = T)
dir.create(dir.gust,showWarnings = F,recursive = T)

## downloading library
## """"""""""""""""""""""""

library(igraph)
library(chron)
library(spatgraphs)
library(classInt)
library(lmom)
library(lmomRFA)
# library(extRFA)



# source functions directory
# """"""""""""""""""""""""""" 
for (func in list.files(dir.functions)) {
  print(func)
  try(source(paste(dir.functions,func,sep="/")))
}


# try(library(data.table))
# # 
# library(ncdf4)    # lecture netcdf
# try(library(RcppRoll))  # needed for roll_mean function in 0.1_facteur_charge_observe_France.R
# try(library(MASS))
# #library(ks)
# try(library(CDFt))
# try(library(zoo))
# try(library(maps))
# try(library(magrittr))
# #BBORGY 08/20/2014
# try(require(sp))
# try(require(mapdata))
# 
# ## Chargement des fonctions
# ## """"""""""""""""""""""""
# 
#for (func in list.files(dir.functions)) {
# print(func)
# try(source(paste0(dir.functions,"/",func)))
#}
# print("OK")
# source(paste0(dir.study,"/0.2_Calcul_FF_region_AJ.R"))

# 
# ## application of local EVA
# if(!("easypackages" %in% installed.packages()[,1])){install.packages("easypackages")}; library(easypackages)
# packages(c("extRemes", "rworldmap", "colorspace", "fields", "RColorBrewer"))
# source("functions/color_palette_rain_maps.R")
# source("functions/fevdBoot.R")
# source("functions/SpacializedMap.R")
# 
# library(extRemes)  #fevd
# 
# 
# library(Matrix)   # Call the newer version of the Matrix package
# library(ismev)
# library(extRemes)
# library(caTools)
# library(Gmisc)
# library(ggplot2)  # wind rose
# library(RColorBrewer) #winde rose
# # RFA
# library(extRFA)
# library(igraph)
# library(spatgraphs)
# library(chron)
