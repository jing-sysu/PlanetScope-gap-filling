function planet_valid=pixel_qualitycontrol(planet_all,nan_band,filt_win,prc_cloud)%mask_backgrd,
[size1,size2,~,size_day]=size(nan_band);
% cloud surrounding pixel
filt=strel('square',filt_win);
planet_nanfilt=imdilate(nan_band,filt);
planet_nanfilt=reshape(planet_nanfilt,[size1*size2,1,size_day]);
% time series outliers
planet_new=planet_all;
planet_new=planet_new-mean(planet_new,1,'omitnan')+mean(planet_new,[1,3],'omitnan'); % centralization
prc_threshold=prctile(planet_new,[prc_cloud/2,100-prc_cloud/2],3);
planet_outlier=planet_new<prc_threshold(:,:,1) | planet_new>prc_threshold(:,:,2);
planet_outlier=sum(planet_outlier,2)>1;
clear planet_new
% detect cloud pixels with planet_valid
for doy=1:size_day
    planet_mis=planet_outlier(:,:,doy);
    mis_pixmph=morphlg(planet_mis,size1,size2,7);% & ~mask_backgrd(:)
    planet_outlier(:,:,doy)=mis_pixmph;
end
% valid index
planet_valid=(1-planet_nanfilt).*(1-planet_outlier);
clear planet_nanfilt planet_outlier 