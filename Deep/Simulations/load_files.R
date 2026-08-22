# load libraries 
library(torch)
library(coro)
library(methods)
library(tidyverse)
library(spatstat)
library(randomForest)
library(xgboost)
library(Matrix)
library(glmnet)
library(tibble)
library(knockoff)

# src from previous paper
src_dir <- "/home/qianzhao_umass_edu/Research/TiltedKnockoff/src/"
source(paste0(src_dir, "sample_obs.R"))
source(paste0(src_dir, "sample_obs_casecontrol.R"))
source(paste0(src_dir, "knockoff_selection.R"))

# deep tilted knockoff functions 
path_deep_tilted <- "/home/qianzhao_umass_edu/Research/TiltedKnockoff/Deep/DeepTiltedKnockoff/"

source(paste0(path_deep_tilted, "filmnet.R"))
source(paste0(path_deep_tilted, "machine_V6.R"))
source(paste0(path_deep_tilted, "stat_V3.R"))
source(paste0(path_deep_tilted, "utils_V2.R"))
source(paste0(path_deep_tilted, "mmd_eval_V2.R"))





