function [targ_img,adj_temp,planet_valid,num_gap]=gapfilling_single(planet_all,targ_img,adj_temp,planet_valid,num_gap,updt_cls,segm,ord_doy,date)
scenerio_threshold=0.8;% threshold for scenerio determination
[size1,size2]=size(updt_cls);
[~,size3,size_day]=size(planet_all);
% adjcent dates 
delt_date=date-date(ord_doy);
[~,adj_ord]=sort(delt_date,'ComparisonMethod','abs');
i=1;
% iteration for gap-filling
while num_gap>0 && i<=min(10,size_day-1)
    i=1+i;
    % missing pixels in the target image
    mis_pixel=isnan(targ_img(:,1));
    % corresponding valid pixels in adjacent dates
    adj_vid=(planet_valid(:,:,adj_ord(i))>0 & planet_valid(:,:,adj_ord(i))<=2);
    target_vid=(planet_valid(:,:,ord_doy)>0 & planet_valid(:,:,ord_doy)<=2);
    % corresponding missing and unmissing valid pixels in adjacent dates
    adj_mis=adj_vid & mis_pixel;
    adjtgt_unmis=adj_vid & target_vid & mis_pixel==0;
    if sum(adj_mis(:))<min(100,sum(mis_pixel)) || sum(adjtgt_unmis(:))<100 % no corresponding valid pixel
        continue;
    end
    % adjcent image for gap-filling the target image with missing pixels
    adj_img=planet_all(:,:,adj_ord(i));
    % coresponding class of missing and unmissing pixel in adjacent image
    num_clsmis=intersect(updt_cls(adj_mis),updt_cls(adjtgt_unmis));
    % for each class build the gap-filling model and apply it to missing pixels
    for z=1:length(num_clsmis)
        cls_clrvid=updt_cls(:)==num_clsmis(z) & adjtgt_unmis; % clear pixels corresponding to the given class 
        cls_misvld=updt_cls(:)==num_clsmis(z) & adj_mis;
        seg_clr=segm(:).*cls_clrvid;% clear pixels corresponding to the given aegmentation 
        seg_mis=segm(:).*cls_misvld;
        subclass_mis=unique(seg_mis(seg_mis>0));
        subclass_unmis=unique(seg_clr(seg_clr>0));
        % centroid of segmentation
        centers=regionprops(reshape(seg_clr,size1,size2),'centroid');
        centerpos = cat(1,centers.Centroid);
        % for each segmentation build the gap-filling model and apply it to missing pixels
        for sub_z=1:length(subclass_mis)
            sub_clrvid=cls_clrvid & seg_clr==subclass_mis(sub_z);
            sub_misvld=cls_misvld & seg_mis==subclass_mis(sub_z);
            meancorr_clear=0;
            numcorr_clear=0;
            % if belong to the same segmentation object and class
            if sum(sub_clrvid)>10
                corr_clear=corr(targ_img(sub_clrvid,:), adj_img(sub_clrvid,:),'Rows','complete');
                meancorr_clear=mean(diag(corr_clear));
                numcorr_clear=sum(diag(corr_clear)>=scenerio_threshold-0.05);
            end
            if meancorr_clear>scenerio_threshold && numcorr_clear==size3
                targ_img=polyfit_adj(targ_img,adj_img,sub_clrvid,sub_misvld);
            else
                % search the closest segmentation objects belonging to the same class
                [disrms_ord,cum_area]=search_closestobj(sub_misvld,cls_clrvid,seg_clr,adj_img,subclass_unmis,centerpos,size1,size2);
                ite_times=0;
                sumcls_clrvid=0;
                ind_clr=zeros(size(cum_area));
                ind_clr(1)=1;
                while (meancorr_clear<scenerio_threshold || numcorr_clear<size3) && sumcls_clrvid<cum_area(end)
                    sumcls_clrvid=1000*(10^ite_times);
                    ite_times=ite_times+1;
                    ind_clr=cum_area<sumcls_clrvid | ind_clr;
                    sub_clrvid=ismember(seg_clr,subclass_unmis(disrms_ord(ind_clr)));
                    if sum(sub_clrvid)>100
                    corr_clear=corr(targ_img(sub_clrvid,:), adj_img(sub_clrvid,:),'Rows','complete');
                    meancorr_clear=mean(diag(corr_clear));
                    numcorr_clear=sum(diag(corr_clear)>=scenerio_threshold-0.05);
                    end
                end
                if meancorr_clear<scenerio_threshold || numcorr_clear<size3
                    continue;
                end                
                targ_img=polyfit_adj(targ_img,adj_img,sub_clrvid,sub_misvld);
            end
            adj_temp(sub_misvld,1,ord_doy)=date(ord_doy);
            adj_temp(sub_misvld,2,ord_doy)=date(adj_ord(i));
            adj_temp(sub_misvld,3,ord_doy)=date(adj_ord(i));
            planet_valid(sub_misvld,1,ord_doy)=2; %planet_valid(cls_misvld,1,adj_ord(i));
        end
    end
    % update the start condition
    num_gap=(size1*size2-sum(targ_img(:,1)>0))/(size1*size2);
end

function targ_clc=polyfit_adj(targ_clc,adj_clc,best_clrvid,sub_misvld)
band=size(targ_clc,2);
for z=1:band
train_x=adj_clc(best_clrvid,z);
train_y=targ_clc(best_clrvid,z);
idxy = isnan(train_y) | isnan(train_x);
[pram,scale,mu]=polyfit(train_x(~idxy),train_y(~idxy),1);
targ_clc(sub_misvld,z)=(adj_clc(sub_misvld,z)-mu(1))./mu(2).*pram(1)+pram(2);
end

function [disrms_ord,cum_area]=search_closestobj(sub_misvld,cls_clrvid,seg_clr,adj_img,subclass_unmis,centerpos,size1,size2)
% sort by distance to missing pixels
[row_mis,col_mis] = find(reshape(sub_misvld,size1,size2)==1);
misCenter=mean([row_mis,col_mis],1);
dis_misclr=pdist2([centerpos(subclass_unmis,2),centerpos(subclass_unmis,1)],misCenter);
[~,dis_ord]=sort(dis_misclr);
areas=hist(seg_clr(seg_clr>0),subclass_unmis);% pixel no. for each subclass_unmis
% sort by RMSE to missing pixels
rms_clear=mean(abs(adj_img(cls_clrvid,:)-mean(adj_img(sub_misvld,:),1)),2,'omitnan');
rms_img=zeros(size1,size2);
rms_img(cls_clrvid)=rms_clear;
rms_mean=regionprops(reshape(seg_clr,size1,size2),reshape(rms_img,size1,size2),'MeanIntensity');
rms_mean = cat(1,rms_mean.MeanIntensity);
[~,rms_ord]=sort(rms_mean(subclass_unmis));
% search the closest subclass_mis by both distance and RMSE
[~,IA,~]=unique(dis_ord);[~,IB,~]=unique(rms_ord);
[~,disrms_ord]=sort(IA+IB);
cum_area=cumsum(areas(disrms_ord));