function [targ_img,adj_temp,planet_valid]=gapfilling_double(planet_all,targ_img,adj_temp,planet_valid,num_gap,nan_prc,prc_thre,updt_cls,segm,ord_doy,date)
[size1,size2]=size(updt_cls);
[~,size3,size_day]=size(planet_all);
% calculate closest dates
delt_date=date-date(ord_doy);
bef_dates=[find(delt_date<0,1,'last'):-1:1,size_day:-1:find(delt_date>0,1,'first')];
aft_dates=[find(delt_date>0,1,'first'):size_day,1:find(delt_date<0,1,'last')];
bef_delt=delt_date(bef_dates);
bef_delt(bef_delt>0)=366-bef_delt(bef_delt>0);
bef_delt=abs(bef_delt);
aft_delt=delt_date(aft_dates);
aft_delt(aft_delt<0)=366+aft_delt(aft_delt<0);
aft_delt=abs(aft_delt);
befaft_delt=bef_delt.*ones(size(bef_delt,1))+aft_delt'.*ones(size(aft_delt,1));
half_befaft=fliplr(triu(fliplr(befaft_delt)));
[delt_val,delt_ord] = sort(half_befaft(:));
[aft_ord,bef_ord]=ind2sub(size(befaft_delt),delt_ord(delt_val>0));
% iteration of interpolation
k=0;
while num_gap>0 && k<length(bef_ord)
    k=k+1;i=aft_ord(k);j=bef_ord(k);
    mis_pixel=isnan(targ_img(:,1));
    % corresponding valid pixels in adjacent dates
    aft_vid=planet_valid(:,:,aft_dates(i))>0 & planet_valid(:,:,aft_dates(i))<=2;
    bef_vid=planet_valid(:,:,bef_dates(j))>0 & planet_valid(:,:,bef_dates(j))<=2;
    % corresponding missing and unmissing pixels in adjacent dates
    befaft_mis=bef_vid & aft_vid & mis_pixel;
    befaft_unmis=bef_vid & aft_vid & mis_pixel==0;
    if sum(befaft_mis(:))==0 || sum(befaft_unmis(:))==0 || nan_prc(aft_dates(i))>=prc_thre || nan_prc(bef_dates(j))>=prc_thre % no corresponding valid pixel
        continue;
    end
    
    % adjcent image
    aft_img=planet_all(:,:,aft_dates(i));
    bef_img=planet_all(:,:,bef_dates(j));
    % coresponding class of missing pixel in adjacent image
    num_clsmis=intersect(updt_cls(befaft_mis),updt_cls(befaft_unmis));
    for z=1:length(num_clsmis)
        cls_clrvid=updt_cls(:)==num_clsmis(z) & befaft_unmis; %  clear pixels corresponding to the given classmin(valid_threshold(adj_ord(i)),0.9)
        cls_misvld=updt_cls(:)==num_clsmis(z) & befaft_mis;
        seg_clr=segm(:).*cls_clrvid;
        seg_mis=segm(:).*cls_misvld;
        subclass_mis=unique(seg_mis(seg_mis>0));
        subclass_unmis=unique(seg_clr(seg_clr>0));
        % centroid of segmentation
        centers=regionprops(reshape(seg_clr,size1,size2),'centroid');
        centerpos = cat(1,centers.Centroid);
        % for each subclass_mis build the relationship and apply to miss pixels
        for sub_z=1:length(subclass_mis)
            sub_clrvid=cls_clrvid & seg_clr==subclass_mis(sub_z);
            sub_misvld=cls_misvld & seg_mis==subclass_mis(sub_z);
            if sum(sub_clrvid)>100% if belong to the same segmentation object and class
                [targ_img,sub_misvld]=regrss_bef_aft(targ_img,bef_img,aft_img,sub_clrvid,sub_misvld,size1,size2);
            else
                % search the closest segmentation objects belonging to the same class
                [disrms_ord,cum_area]=search_closestobj(sub_misvld,cls_clrvid,seg_clr,aft_img,bef_img,subclass_unmis,centerpos,size1,size2);
                ind_clr=zeros(size(cum_area));
                ind_clr(1)=1;
                ind_clr=cum_area<1000 | ind_clr;
                sub_clrvid=ismember(seg_clr,subclass_unmis(disrms_ord(ind_clr)));
                [targ_img,sub_misvld]=regrss_bef_aft(targ_img,bef_img,aft_img,sub_clrvid,sub_misvld,size1,size2);
            end
            adj_temp(sub_misvld,1,ord_doy)=date(ord_doy);
            adj_temp(sub_misvld,2,ord_doy)=date(bef_dates(j));
            adj_temp(sub_misvld,3,ord_doy)=date(aft_dates(i));
            planet_valid(sub_misvld,1,ord_doy)=3;%planet_valid1(cls_misvld,1,aft_dates(i))*0.8;
        end
    end
    % update the start condition
    num_gap=(size1*size2-sum(targ_img(:,1)>0))/(size1*size2);
    if num_gap<0.005
        cls_misvld=isnan(targ_img(:,1));
        targ_img_mat=reshape(targ_img,size1,size2,size3);
        targ_img_mat=fillmissing(targ_img_mat,'nearest');
        targ_img=reshape(targ_img_mat,[size1*size2,size3]);
        adj_temp(cls_misvld,1,ord_doy)=date(ord_doy);
        adj_temp(cls_misvld,2,ord_doy)=date(ord_doy);
        adj_temp(cls_misvld,3,ord_doy)=date(ord_doy);
        planet_valid(cls_misvld,1,ord_doy)=4;%planet_valid1(cls_misvld,1,aft_dates(i))*0.8;
        num_gap=(size1*size2-sum(targ_img(:,1)>0))/(size1*size2);
    end
end

function [targ_img,sub_misvld]=regrss_bef_aft(targ_img,bef_img,aft_img,sub_clrvid,sub_misvld,size1,size2)
band=size(targ_img,2);
delt_mis=bef_img(sub_misvld,:)-aft_img(sub_misvld,:);
delt_unmis=bef_img(sub_clrvid,:)-aft_img(sub_clrvid,:);
invld=sum(delt_mis>mean(delt_unmis)+2*std(delt_unmis) | delt_mis<mean(delt_unmis)-2*std(delt_unmis),2)<=1;
sub_misvld(sub_misvld>0)=invld;
for z=1:band
% bef_img1=imfilter(reshape(bef_img,size1,size2),fspecial('gaussian',[15 15],5));
% bef_difimg=reshape(bef_img,size1,size2)-bef_img1;
% bef_difimg=bef_difimg(:);
% bef_difimg(isnan(bef_difimg))=0;
% aft_img1=imfilter(reshape(aft_img,size1,size2),fspecial('gaussian',[15 15],5));
% aft_difimg=reshape(aft_img,size1,size2)-aft_img1;
% aft_difimg=aft_difimg(:);
% aft_difimg(isnan(aft_difimg))=0;
train_x=bef_img(sub_clrvid,z)-aft_img(sub_clrvid,z);
train_y=targ_img(sub_clrvid,z)-aft_img(sub_clrvid,z);
train_xx=[ones(size(train_x,1),1),train_x];
[pram,~,~,~,stats]=regress(train_y,train_xx);
if stats(3)<0.05 || pram(1)==0
    bef_misvld=bef_img(sub_misvld,z);
    aft_misvld=aft_img(sub_misvld,z);
    test_x=bef_misvld-aft_misvld;
    test_xx=[ones(size(test_x,1),1),test_x];
    targ_img(sub_misvld,z)=test_xx*pram+aft_img(sub_misvld,z);%+pram(2)*bef_difimg(sub_misvld)+(1-pram(2))*aft_difimg(sub_misvld)
else
    train_x1=bef_img(sub_clrvid,z)-targ_img(sub_clrvid,z);
    train_x2=aft_img(sub_clrvid,z)-targ_img(sub_clrvid,z);
    test_x1=bef_img(sub_misvld,z)-mean(train_x1);
    test_x2=aft_img(sub_misvld,z)-mean(train_x2);
    targ_img(sub_misvld,z)=mean(abs(train_x2))./(mean(abs(train_x1))+mean(abs(train_x2))).*test_x1+mean(abs(train_x1))./(mean(abs(train_x1))+mean(abs(train_x2))).*test_x2;
end
end

function [disrms_ord,cum_area]=search_closestobj(sub_misvld,cls_clrvid,seg_clr,aft_img,bef_img,subclass_unmis,centerpos,size1,size2)
% sort by distance to missing pixels
[row_mis,col_mis] = find(reshape(sub_misvld,size1,size2)==1);
misCenter=mean([row_mis,col_mis],1);
dis_misclr=pdist2([centerpos(subclass_unmis,2),centerpos(subclass_unmis,1)],misCenter);
[~,dis_ord]=sort(dis_misclr);
areas=hist(seg_clr(seg_clr>0),subclass_unmis);
% sort by RMSE to missing pixels
rms_clear=mean((aft_img(cls_clrvid,:)-mean(aft_img(sub_misvld,:),1)).^2,2,'omitnan')+mean((bef_img(cls_clrvid,:)-mean(bef_img(sub_misvld,:),1)).^2,2,'omitnan');
rms_img=zeros(size1,size2);
rms_img(cls_clrvid)=rms_clear;
rms_mean=regionprops(reshape(seg_clr,size1,size2),reshape(rms_img,size1,size2),'MeanIntensity');
rms_mean = cat(1,rms_mean.MeanIntensity);
[~,rms_ord]=sort(rms_mean(subclass_unmis));
% search the closest subclass_mis by both distance and RMSE
[~,IA,~]=unique(dis_ord);[~,IB,~]=unique(rms_ord);
[~,disrms_ord]=sort(IA+IB);
cum_area=cumsum(areas(disrms_ord));