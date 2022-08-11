%% --------------------------------------------------------------------------
% Jing Wang wrote it

% Function: write pixel quality index to disk
% Usuage:
% Input parameters:
% Planet_new is the pixel quality index
% b is the block no.
% write_dir, planet_dir are the output file dir and file name of pixel quality index
% geo_info, geoTags are map info
% str_suffix is suffix of file name
% Output:
% Write the pixel quality index to disk

%% --------------------------------------------------------------------------
function write_interpindexPlanet(Planet_new,b,write_dir,planet_dir,geo_info,geoTags,str_suffix)
color_map=[1 0 0;0 1 0;0 0 1];%red,green,blue
[~,~,size_day]=size(Planet_new);
for d=1:size_day
    data=Planet_new(:,:,d);
    Planet_interp_index_color = ind2rgb(data,color_map);%convert to RGB image
    folder=strcat(write_dir,planet_dir(d).name,'\temp\');
    if ~exist(folder,'dir')
        mkdir(folder);
    end
    geotiffwrite(strcat(folder,planet_dir(d).name(14:21),str_suffix,'_block_',num2str(b),'.tif'),Planet_interp_index_color,geo_info,'GeoKeyDirectoryTag',geoTags);
end