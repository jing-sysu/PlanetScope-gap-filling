function [img_zerosele,date_zerosele,cls_zerosele,num_optcls]=derive_cls_zerosele(nan_zero,ord_valid,planet_all,size1,size2,size3,date,num_optcls,klist)%     klist=3:6;%5:10;%the number of clusters you want to try
se=strel('square',3);%Morphological structuring element
if isempty(nan_zero)==1
    nan_zero=ord_valid(1);
    img_zerosele=planet_all(:,:,nan_zero);
    myfunc = @(X,K)(kmeans(X, K));
    eva = evalclusters(img_zerosele,myfunc,'CalinskiHarabasz','klist',klist);
    if isnan(eva.OptimalK)==0
        num_optcls=eva.OptimalK;
    end
    cls_zerosele=kmeans(img_zerosele,num_optcls,'MaxIter',1000);%
    MAX_nan=sum(isnan(cls_zerosele))+1;
    [~,reord]=sort(abs(date-date(nan_zero)));
    planet_new=planet_all(:,:,reord);
    while sum(isnan(cls_zerosele))>0 && sum(isnan(cls_zerosele))<MAX_nan
        [other_nan,other_ord]=sort(sum(isnan(planet_new(isnan(cls_zerosele),1,:)),1));
        cls_othersele=kmeans(planet_new(:,:,other_ord(1)),num_optcls,'MaxIter',1000);%;
        [img_zerosele,cls_zerosele]=fill_closepix(img_zerosele,cls_zerosele,planet_new(:,:,other_ord(1)),cls_othersele,size1,size2);
        % img_zerosele=fill_closepix(img_zerosele,planet_new(:,:,other_ord(1)),cls_othersele,size1,size2,size3,num_optcls);
        % if sum(isnan(cls_zerosele))<10
        %     cls_zerosele(img_zerosele(:,1)~=1)=fillmissing(cls_zerosele(img_zerosele(:,1)~=1),'previous');
        % end
        MAX_nan=max(other_nan(:));
    end
    date_zerosele=date(nan_zero);
else
    idx = randperm(length(nan_zero));
    nan_zerosele=nan_zero(idx(1:min(5,length(nan_zero))));
    img_zerosele=reshape(planet_all(:,:,nan_zerosele),size1*size2,min(5,length(nan_zero))*size3);
    myfunc = @(X,K)(kmeans(X, K));
    eva = evalclusters(img_zerosele,myfunc,'CalinskiHarabasz','klist',klist);
    if isnan(eva.OptimalK)==0
        num_optcls=eva.OptimalK;
    end
    cls_zerosele=kmeans(img_zerosele,num_optcls,'MaxIter',1000);%

    img_zerosele=planet_all(:,:,nan_zero);
    date_zerosele=date(nan_zero);
end
cls_zerosele=reshape(cls_zerosele,size1,size2);
img_open=imopen(cls_zerosele,se);
cls_zerosele=imclose(img_open,se);%remove tiny points
cls_zerosele=cls_zerosele(:);

function [curimg_vect,curcls_vect]=fill_closepix(curimg_vect,curcls_vect,aftimg_vect,aftcls_vect,size1,size2)
mis_mat=reshape(isnan(curcls_vect),size1,size2);
if sum(mis_mat(:))>0
    se = strel('square', 31);
    mis_dilmat=imdilate(mis_mat,se);
    num_cls=unique(aftcls_vect(mis_mat(:) & aftcls_vect>0));
    aftcls_dilvet=aftcls_vect(mis_dilmat & ~mis_mat);
    curcls_dilvet=curcls_vect(mis_dilmat & ~mis_mat);
    for i=1:length(num_cls)
        curcls_i=curcls_dilvet(aftcls_dilvet==num_cls(i));
        curcls_iput=mode(curcls_i,1);
        curcls_vect(mis_mat(:)==1 & aftcls_vect==num_cls(i))=curcls_iput;
        
        curimg_dilvet=curimg_vect(mis_dilmat(:) & ~mis_mat(:) & curcls_vect==curcls_iput,:);
        aftimg_dilvet=aftimg_vect(mis_dilmat(:) & ~mis_mat(:) & aftcls_vect==num_cls(i),:);
        test_y=mean(curimg_dilvet,1,'omitnan')-mean(aftimg_dilvet,1,'omitnan')+aftimg_vect(mis_mat(:)==1 & aftcls_vect==num_cls(i),:);
        curimg_vect(mis_mat(:)==1 & aftcls_vect==num_cls(i),:)=test_y;
    end
end
