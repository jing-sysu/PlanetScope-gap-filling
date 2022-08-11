%% --------------------------------------------------------------------------
% JingWang wrote it

% Function: The second step of Planet data processing, namely mosaic Planet images within the same ROI
% Usuage:
% Input parameters: 
% file_dir & ROI region (with upleft and lowright latitude and longtitude)
% Output:
% mosaic image   

%% --------------------------------------------------------------------------
function Planet_mosaic_v3(planet_dir,str_suffix)
%image mosaic
for i=1:length(planet_dir)
    cloudy_dir=dir(strcat(planet_dir(i).folder,'\',planet_dir(i).name,'\*_cloudseriesnoqc.tif'));%sub dir
    [cloudy_image,mosaic_info]=geotiffread(strcat(cloudy_dir(1).folder,'\',cloudy_dir(1).name));
    [size1,size2,dim]=size(cloudy_image);

    sub_dir=dir(strcat(cloudy_dir(1).folder,'\temp\*_interp_block_*.tif'));%sub dir
    [~,idx] = sort([sub_dir.datenum]);
    sub_dir = sub_dir(idx);
    if isempty(sub_dir)
        disp(['No data for ',planet_dir(i).name]);
        continue;%if no SR data,it will be not usefull
    end
    %mosaic map info
    [planet_image1,sub_info1]=geotiffread(strcat(sub_dir(1).folder,'\',sub_dir(1).name));
    info = geotiffinfo(strcat(sub_dir(1).folder,'\',sub_dir(1).name));
    geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
    mosaic_image=zeros(size1,size2,dim,2);
    [row_1_lu,col_1_lu]=map2pix(mosaic_info,sub_info1.XWorldLimits(1),sub_info1.YWorldLimits(2));
    [row_1_rd,col_1_rd]=map2pix(mosaic_info,sub_info1.XWorldLimits(2),sub_info1.YWorldLimits(1));
    
    row_1_lu=ceil(row_1_lu);
    col_1_lu=ceil(col_1_lu);
    row_1_rd=floor(row_1_rd);
    col_1_rd=floor(col_1_rd);%change to integer
    
    mosaic_image(row_1_lu:row_1_rd,col_1_lu:col_1_rd,:,1)=planet_image1;
    for sub_i=2:length(sub_dir)
        [planet_image,sub_info]=geotiffread(strcat(sub_dir(sub_i).folder,'\',sub_dir(sub_i).name));%planet_image1,sub_info1
        %change to row,col coordinate
        [row_1_lu,col_1_lu]=map2pix(mosaic_info,sub_info.XWorldLimits(1),sub_info.YWorldLimits(2));
        [row_1_rd,col_1_rd]=map2pix(mosaic_info,sub_info.XWorldLimits(2),sub_info.YWorldLimits(1));
        
        row_1_lu=ceil(row_1_lu);
        col_1_lu=ceil(col_1_lu);
        row_1_rd=floor(row_1_rd);
        col_1_rd=floor(col_1_rd);%change to integer
        
        mosaic_image(row_1_lu:row_1_rd,col_1_lu:col_1_rd,:,2)=planet_image;
        
        %     histogram matching
        mosaic_image(mosaic_image==0)=nan;
        ref_img=squeeze(mosaic_image(:,:,:,1));
        targ_img=squeeze(mosaic_image(:,:,:,2));
        for band=1:dim
            overlap=ref_img(:,:,band)>0 & targ_img(:,:,band)>0;
            if sum(overlap(:)>0)
                input_band=targ_img(:,:,band);
                [mu_input,sigma_input] = normfit(input_band(overlap==1));
                ref_band=ref_img(:,:,band);
                [mu_ref,sigma_ref] = normfit(ref_band(overlap==1));
                param_a=sigma_ref/sigma_input;
                param_b=mu_ref-param_a*mu_input;
                targ_img(:,:,band)=input_band.*param_a+param_b;
            end
        end
        mosaic_image(:,:,:,2)=targ_img;
        mosaic_image(:,:,:,1)=mean(mosaic_image,4,'omitnan');
        mosaic_image(:,:,:,2)=zeros(mosaic_info.RasterSize(1),mosaic_info.RasterSize(2),dim);
    end
    %     mosaic image with average
    mosaic_image_comb=mosaic_image(:,:,:,1);% average in overlap
    mosaic_image_comb(isnan(mosaic_image_comb))=0;
    if sum(mosaic_image_comb(:,:,1)>0,[1,2])/(mosaic_info.RasterSize(1)*mosaic_info.RasterSize(2))<0.2
        disp(['No enough data for ',planet_dir(i).name]);
    end
    geotiffwrite(strcat(cloudy_dir(1).folder,'\',sub_dir(1).name(1:8),str_suffix,'.tif'),mosaic_image_comb,mosaic_info,'GeoKeyDirectoryTag',geoTags);
    clear mosaic_band mosaic_image
end