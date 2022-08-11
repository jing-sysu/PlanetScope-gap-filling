%% --------------------------------------------------------------------------
% Jing Wang wrote it

% Function: read block image time series
% Usuage:
% Input parameters:
% planet_file is the input file dir of image time series
% block_rowstart,block_rowend,block_colstart,block_colend are spatial extent of the block
% str_suffix is suffix of file name
% Output:
% planet_all is image time series block
% date is date info
% sub_info,geoTags are map info

%% --------------------------------------------------------------------------
function [planet_all,date,sub_info,geoTags]=read_series_block(planet_file,block_rowstart,block_rowend,block_colstart,block_colend,str_suffix)
size_day=length(planet_file);
planet_all=zeros(block_rowend-block_rowstart+1,block_colend-block_colstart+1,4,size_day,'single');
date=zeros(size_day,1);
for i=1:size_day
    sub_dir=dir([planet_file(i).folder,'\',planet_file(i).name,'\*',str_suffix,'.tif']);%sub dir
    if isempty(sub_dir)==1% if no enough information
        continue;
    end
    % read image
    [planet_image,geo_info]=geotiffread(strcat(planet_file(i).folder,'\',planet_file(i).name,'\',sub_dir(1).name));%
    planet_loc=planet_image(block_rowstart:block_rowend,block_colstart:block_colend,:);%;%
    % geoinformation
    info = geotiffinfo(strcat(planet_file(i).folder,'\',planet_file(i).name,'\',sub_dir(1).name));
    geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
    % date information, DOY
    date(i)=datenum(sub_dir(1).name(1:8),'yyyymmdd')-datenum([sub_dir(1).name(1:4),'0101'],'yyyymmdd')+1;
    
    planet_all(:,:,:,i)=planet_loc;
    clear planet_image planet_loc
end
[size1,size2,~,~]=size(planet_all);
% write geoinformation
sub_info=geo_info;
sub_info.XWorldLimits(1)=geo_info.XWorldLimits(1)+(block_colstart-1)*geo_info.CellExtentInWorldX;
sub_info.XWorldLimits(2)=geo_info.XWorldLimits(1)+block_colend*geo_info.CellExtentInWorldX;
sub_info.YWorldLimits(2)=geo_info.YWorldLimits(2)-(block_rowstart-1)*geo_info.CellExtentInWorldY;
sub_info.YWorldLimits(1)=geo_info.YWorldLimits(2)-block_rowend*geo_info.CellExtentInWorldY;%correct the bug, not need to take attention for input_ROI
sub_info.RasterSize= [size1,size2];