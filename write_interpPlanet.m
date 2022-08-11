%% --------------------------------------------------------------------------
% Jing Wang wrote it

% Function: write Planet interpolation results to disk
% Usuage:
% Input parameters:
% Planet_new is the interpolation results
% b is the block no.
% write_dir, planet_dir are the output file dir and file name of interpolation results
% geo_info, geoTags are map info
% str_suffix is suffix of file name
% Output:
% Write interpolation results to disk

%% --------------------------------------------------------------------------
function write_interpPlanet(Planet_new,b,write_dir,planet_dir,geo_info,geoTags,str_suffix)
[~,~,~,size_day]=size(Planet_new);
for d=1:size_day
    Planet_temp=Planet_new(:,:,:,d);
    folder=strcat(write_dir,planet_dir(d).name,'\temp\');
    if ~exist(folder,'dir')
        mkdir(folder);
    end
    geotiffwrite(strcat(folder,planet_dir(d).name(14:21),str_suffix,'_block_',num2str(b),'.tif'),Planet_temp,geo_info,'GeoKeyDirectoryTag',geoTags);
end
