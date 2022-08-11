# PlanetScope-gap-filling
## --------------- Description ----------------------------------------------
# "autoOCBGF_all" object-class based gap-filling (OCBGF) method used for gap-filling the entire time series.
# Please note that this version is not used for water body yet.
#
# Please refer to :
#
# Jing Wang, Calvin KF Lee, Xiaolin Zhu, Ruyin Cao, Yating Gu, Shengbiao Wu, Jin Wu*,
# "A new object-class based gap-filling method for PlanetScope satellite image time series", 
# Remote sensing of Enviroment, Accepted in 16/06/2022
# 
## --------------- Usuage----------------------------------------------------
#     
# Main Function : autoOCBGF_all
# 
# This main function includes three steps : 
# Step 1 : Read PlanetScope time-series images one by one from the folder.
#          The time-series images should have same spatial extent (row * column * band).
# Step 2 : Conduct OCBGF with four sub-steps.
#          Sub-step1 : Pixel-level quality control
#          Sub-step2 : Object segmentation and classification 
#          Sub-step3 : Scenario-specific gap-filling
#          Sub-step4 : Post-image-processing using guided filter
# Step 3 : Write each gap-filled result in original folder of each PlanetScope time-series image.
#
# Input arguments :
# input_dir : The folder of PlanetScope time-series files, each file is named as "planet_order_*", e.g. "planet_order_20180108",
#             each file includes the to-be-gapfilled PlanetScope image "*_cloudseriesnoqc.tif" 
# size_imgrow, size_imgcol : The row (column) number of each PlanetScope image. The row and column numbers should be consistent for the entire time series. 
# filt_win : The size of the structure element for pixel quality control; e.g. 5 
# prc_cloud : The temporal percentile threshold for pixel quality control; e.g. 1
# klist : The range of the number of optimal classes 
# 
# Please note that due to the MATLAB memory limitation, OCBGF divided PlanetScope image time series into discrete cubic blocks, 
#            for example, a block with a spatial extent of 1111 * 1111 PlanetScope pixels (3333m * 3333m). 
#            
# Output arguments :
# planet_all: The gap-filled results; Write gap-filled results to the folders of the specified image(s) by "*_interp"
# planet_valid: Pixel quality index (1: valid; 2: gap-filled with Single-reference-image Scenario; 3: gap-filled with Two-reference-images Scenario;4: gap-filled with others); As you like, you can also write planet_valid by "*_index"
# adj_temp: Adjcent dates (DOY;day-of-year) used for gap-filling each missing pixel; As you like, you can also write adj_temp by "*_temp"
# 
# Example Data in Dropbox: https://www.dropbox.com/sh/r6x82vsf0cg2j93/AAAEsRV9DFif4Tjb_xM1uwrfa?dl=0
# -------------------------------------------------------
# Author: Jing Wang (w1012j@163.com)
# Last Date: 19/06/2022
