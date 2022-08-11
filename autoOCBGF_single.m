%% --------------- Description ----------------------------------------------
% "autoOCBGF_single" object-class based gap-filling (‘OCBGF’) method used for gap-filling specified image(s) in a time series.
% Please note that this version is not used for water body yet.
%
% Please refer to :
%
% Jing Wang, Calvin KF Lee, Xiaolin Zhu, Ruyin Cao, Yating Gu, Shengbiao Wu, Jin Wu*,
% "A new object-class based gap-filling method for PlanetScope satellite image time series", 
% Remote sensing of Enviroment, Accepted in 16/06/2022
% 
%% --------------- Usuage----------------------------------------------------
%     
% Function : autoOCBGF_single
% 
% This main function includes three steps : 
% Step 1 : Read PlanetScope time-series images one by one from the folder.
%          The time-series images should have same spatial extent (row * column * band).
% Step 2 : Conduct OCBGF with four sub-steps.
%          Sub-step1 : Pixel-level quality control
%          Sub-step2 : Object segmentation and classification 
%          Sub-step3 : Scenario-specific gap-filling
%          Sub-step4 : Post-image-processing using guided filter
% Step 3 : Write each gap-filled result in original folder of each PlanetScope time-series image.
%
% Input arguments :
% input_dir : The folder of PlanetScope time-series files, each file is named as "planet_order_*", e.g. "planet_order_20180108",
%             each file includes the to-be-gapfilled PlanetScope image "*_cloudseriesnoqc.tif" 
% ind_single : The image(s) no. of specified target image(s)  
% size_imgrow, size_imgcol : The row (column) number of each PlanetScope image. The row and column numbers should be consistent for the entire time series. 
% filt_win : The size of the structure element for pixel quality control; e.g. 5 
% prc_cloud : The temporal percentile threshold for pixel quality control; e.g. 1
% klist : The range of the number of optimal classes 
% 
% Please note that due to the MATLAB memory limitation, OCBGF divided PlanetScope image time series into discrete cubic blocks, 
%            for example, a block with a spatial extent of 1111 * 1111 PlanetScope pixels (3333m * 3333m). 
%            
% Output arguments :
% planet_all: The gap-filled results; Write gap-filled results to the folders of the specified image(s) by "*_interp"
% planet_valid: Pixel quality index (1: valid; 2: gap-filled with Single-reference-image Scenario; 3: gap-filled with Two-reference-images Scenario;4: gap-filled with others); As you like, you can also write planet_valid by "*_index"
% adj_temp: Adjcent dates (DOY;day-of-year) used for gap-filling each missing pixel; As you like, you can also write adj_temp by "*_temp"
% 
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 19/06/2022
% -------------------------------------------------------
%% --------------------------------------------------------------------------
function autoOCBGF_single
clear
clc
tic
%% parameter setting; an example for Iowa-cropland site
input_dir='E:\OCBGF_exampledata\';% file name of the time series
dir_planet=dir([input_dir,'planet_order_*']);% file dir of the time series
ind_single=1:length(dir_planet);%[38];% specify the target image(s) no. in the time series
filt_win=9;% the size of the structure element
prc_cloud=1;% the percentage threshold
klist=8:10;% the range of the number of optimal classes (K range)
size_imgrow=3334;% row number of image
size_imgcol=3334;% col number of image
%%
num_optcls=max(klist);% the default number of classes
write_dir=input_dir;
size_day=length(dir_planet);
rownum_block=round(size_imgrow/1111);% the number of row blocks
colnum_block=round(size_imgcol/1111);% the number of col blocks
for block=1:rownum_block*colnum_block% for each block
    % spatial extent of the block
    [block_row,block_col]=ind2sub([rownum_block,colnum_block],block);
    bksize_rows=ceil(size_imgrow/rownum_block);
    bksize_cols=ceil(size_imgcol/colnum_block);
    block_rows=max((block_row-1)*bksize_rows+1-50,1);
    block_rowe=min(block_row*bksize_rows+50,size_imgrow);
    block_cols=max((block_col-1)*bksize_cols+1-50,1);
    block_cole=min(block_col*bksize_cols+50,size_imgcol);
    
    %% Step1: read image time series
    [planet_all,date,sub_info,geoTags]=read_series_block(dir_planet,block_rows,block_rowe,block_cols,block_cole,'cloudseriesnoqc');
    [size1,size2,size3,~]=size(planet_all);
    
    % mask water body
    nan_band=planet_all(:,:,1,:)<=0 | isnan(planet_all(:,:,1,:));
    mask_backgrd=sum(nan_band,4)==size_day;% water body
    planet_all(repmat(nan_band,[1,1,size3,1]))=nan;
    nan_band=nan_band & ~mask_backgrd;
    planet_all=reshape(planet_all,[size1*size2,size3,size_day]);
    
    %% Step2: Conduct OCBGF; Sub-step1: Pixel-level quality control
    planet_valid=pixel_qualitycontrol(planet_all,nan_band,filt_win,prc_cloud);
    planet_all(repmat(~planet_valid,[1,size3,1]))=nan;
    % percentage of background (exclude water body)
    mask_backgrd=sum(isnan(planet_all(:,1,:)),3)==size_day;
    planet_valid(repmat(mask_backgrd,[1,1,size_day]))=nan;
    planet_all(repmat(mask_backgrd,[1,size3,size_day]))=1;
    % percentage of whole image contaminated by clouds
    nan_prc=sum(isnan(planet_all(:,1,:)),1)./(size1*size2);
    [num_nan,ord_valid]=sort(nan_prc(:));
    nan_zero=find(nan_prc<=0);
    % derive clearest image(s) and class map(s) among the time series
    [imgall_zerosele,date_zerosele,cls_zerosele,num_optcls]=derive_cls_zerosele(nan_zero,ord_valid,planet_all,size1,size2,size3,date,num_optcls,klist);    
    adj_temp=zeros([size1*size2,3,size_day]);% adjcent dates used to gap fill missing information
    prc_thre=0.99-sum(mask_backgrd(:))/(size1*size2);
    % gap-fill process for each image with gaps
    for ii=1:length(ind_single)
        doy=find(ord_valid==ind_single(ii))
        if num_nan(doy)<=0% clear image
            continue;
        elseif num_nan(doy)>=prc_thre % NAN for entire image
            cur_date=ord_valid(doy); % current order of this image
            [targ_img,adj_temp]=gapfilling_weight(planet_all,cur_date,adj_temp,size1,size2,date);
            planet_all(:,:,cur_date)=targ_img;
        else
            targ_img=planet_all(:,:,ord_valid(doy));% target image for interpolation
            num_gap=num_nan(doy);
            
            delt_date=date-date(ord_valid(doy));
            bef_dates=[find(delt_date<0,1,'last'):-1:1,size_day:-1:find(delt_date==0)];
            aft_dates=[find(delt_date>0,1,'first'):size_day,1:find(delt_date==0)];
            %% Sub-step2:Object segmentation and classification 
            adj_ord1=find(nan_prc(bef_dates(1:min(5,size_day-1)))<0.05);
            adj_ord2=find(nan_prc(aft_dates(1:min(5,size_day-1)))<0.05);
            if isempty(adj_ord1)
                [~,adj_ord1]=sort(nan_prc(bef_dates(1:min(5,size_day-1))));
            end
            if isempty(adj_ord2)
                [~,adj_ord2]=sort(nan_prc(aft_dates(1:min(5,size_day-1))));
            end
            adj_max=planet_all(:,:,[bef_dates(adj_ord1(1)),aft_dates(adj_ord2(1))]);
            [~,ind_zerosele]=min(abs(date_zerosele-date(ord_valid(doy))));
            img_zerosele=imgall_zerosele(:,:,ind_zerosele);
            
            [updt_cls,segm]=search_closepix4(targ_img,adj_max,img_zerosele,cls_zerosele,size1,size2,num_optcls,klist);
            %% Sub-step3: Scenario-specific gap-filling  
            % Scenario 1 Single-reference-image            
            [targ_img,adj_temp,planet_valid,num_gap]=gapfilling_single(planet_all,targ_img,adj_temp,planet_valid,num_gap,updt_cls,segm,ord_valid(doy),date);
            % Scenario 2 Two-reference-image for remain missing data
            [targ_img,adj_temp,planet_valid]=gapfilling_double(planet_all,targ_img,adj_temp,planet_valid,num_gap,nan_prc,prc_thre,updt_cls,segm,ord_valid(doy),date);
            % backup scenario not falling into scenario 1 or 2
            [targ_img,adj_temp,planet_valid]=gapfilling_backup(planet_all,targ_img,adj_temp,planet_valid,num_gap,nan_prc,prc_thre,updt_cls,segm,ord_valid(doy),date);
            planet_all(:,:,ord_valid(doy))=targ_img;
            clear targ_img targ_img1 targ_img2 weig_img
        end
    end
    planet_all(planet_all==1)=nan;
    planet_all=reshape(planet_all,[size1,size2,size3,size_day]);
%     planet_valid=reshape(planet_valid,[size1,size2,size_day]);
%     adj_temp=reshape(adj_temp,[size1,size2,3,size_day]);
    %% Sub-step4: Post-image-processing using guided filter
    filter_ref=reshape(imgall_zerosele(:,:,1),size1,size2,size3);
    for doy=1:size_day
        filter_targ=planet_all(:,:,:,doy);
        planet_all(:,:,:,doy)= imguidedfilter(filter_targ,filter_ref);
    end
    %% Step3: write the interpolation results
    write_interpPlanet(planet_all(:,:,:,ind_single),block,write_dir,dir_planet(ind_single),sub_info,geoTags,'_interp');
%     write_interpindexPlanet(planet_valid(:,:,ind_single),block,write_dir,dir_planet(ind_single),sub_info,geoTags,'_index');
%     write_interptempPlanet(adj_temp(:,:,:,ind_single),block,write_dir,dir_planet(ind_single),sub_info,geoTags,'_temp');
    clear planet_new planet_all planet_all1 planet_all2 planet_valid planet_valid1 planet_valid2 adj_temp adj_temp1 adj_temp2
end
Planet_mosaic_v3(dir_planet(ind_single),'_mosaicinterp2');
toc
end