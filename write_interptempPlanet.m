%% --------------------------------------------------------------------------
% Jing Wang wrote it

% Function: write adjcent dates used for gap-filling to disk
% Usuage:
% Input parameters:
% Planet_new is the adjcent dates used for gap-filling
% b is the block no.
% write_dir, planet_dir are the output file dir and file name of adjcent dates used for gap-filling
% geo_info, geoTags are map info
% str_suffix is suffix of file name
% Output:
% Write adjcent dates used for gap-filling to disk

%% --------------------------------------------------------------------------
function write_interptempPlanet(Planet_new,b,write_dir,planet_dir,geo_info,geoTags,str_suffix)
[~,~,~,size_day]=size(Planet_new);
for d=1:size_day
    Planet_temp=Planet_new(:,:,:,d);
    folder=strcat(write_dir,planet_dir(d).name,'\temp\');
    if ~exist(folder,'dir')
        mkdir(folder);
    end
    geotiffwrite(strcat(folder,planet_dir(d).name(14:21),str_suffix,'_block_',num2str(b),'.tif'),Planet_temp,geo_info,'GeoKeyDirectoryTag',geoTags);
end