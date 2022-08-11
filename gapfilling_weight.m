function [targ_img,adj_temp]=gapfilling_weight(planet_all,cur_date,adj_temp,size1,size2,date)%cur_date: current order of this image
nan_prc=sum(isnan(planet_all(:,1,:)),1)./(size1*size2);
nonnan=find(nan_prc==0); % other not entire NAN image but has been gap-filled
if isempty(nonnan)==1
    [num_nan1,ord_valid1]=sort(nan_prc(:));
    nonnan=[ord_valid1(1),ord_valid1(2)];
end
delt_date=date(nonnan)-date(cur_date); % difference of adjacent dates
aft_date=nonnan(find(delt_date>0,1,'first')); % the image order before the current image
bef_date=nonnan(find(delt_date<0,1,'last')); % the image order after the current image
if isempty(aft_date)==1 % when the current image has no after date image
    aft_date=nonnan(1);
    weight_aft=(date(cur_date)-date(bef_date))/(date(aft_date)+365-date(bef_date));
    weight_bef=(date(aft_date)+365-date(cur_date))/(date(aft_date)+365-date(bef_date));
elseif isempty(bef_date)==1 % when the current image has no before date image
    bef_date=nonnan(end);
    weight_aft=(date(cur_date)+365-date(bef_date))/(date(aft_date)+365-date(bef_date));
    weight_bef=(date(aft_date)-date(cur_date))/(date(aft_date)+365-date(bef_date));
else
    weight_aft=(date(cur_date)-date(bef_date))/(date(aft_date)-date(bef_date));
    weight_bef=(date(aft_date)-date(cur_date))/(date(aft_date)-date(bef_date));
end
% weoght average of adjacent images
aft_img=planet_all(:,:,aft_date);
bef_img=planet_all(:,:,bef_date);
weig_img(:,:,1)=aft_img.*weight_aft;
weig_img(:,:,2)=bef_img.*weight_bef;
targ_img=sum(weig_img,3,'omitnan');
adj_temp(:,1,cur_date)=date(cur_date);
adj_temp(:,2,cur_date)=date(bef_date);
adj_temp(:,3,cur_date)=date(aft_date);