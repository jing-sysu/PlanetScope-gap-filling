function [gapfill_cls,dilLabl_complt]=search_closepix4(curimg_vect,adjimgs_vect,fullimg_vect,fullcls_vect,size1,size2,num_optcls,klist)
% tic
scale=1/4;
se=strel('square',3);%Morphological structuring element
% Step1: classify for adjacent images
[~,size_band,size_dim]=size(adjimgs_vect);
adjimgs_mat=reshape(adjimgs_vect,[size1,size2,size_band*size_dim]);
adjimgs_resmat=imresize(adjimgs_mat,scale);
[newsize1,newsize2,~]=size(adjimgs_resmat);
adjimgs_resvect=reshape(adjimgs_resmat,[newsize1*newsize2,size_band*size_dim]);
myfunc = @(X,K)(kmeans(X, K));
eva = evalclusters(adjimgs_resvect,myfunc,'CalinskiHarabasz','klist',klist);
if isnan(eva.OptimalK)==0
    num_optcls=eva.OptimalK;
end
adjcls_vect=kmeans(adjimgs_resvect,num_optcls,'MaxIter',1000);%
fullcls_mat=reshape(fullcls_vect,[size1,size2]);
fullcls_resmat=imresize(fullcls_mat,[newsize1,newsize2],'nearest');
fullimg_mat=reshape(fullimg_vect,[size1,size2,size_band]);
fullimg_resmat=imresize(fullimg_mat,scale);
fullimg_resvet=reshape(fullimg_resmat,[newsize1*newsize2,size_band]);
if sum(adjcls_vect>0)/(newsize1*newsize2)>0.8%0.6%
%     adjcls_vect=search_closepix(adjcls_vect,fullcls_resmat(:),newsize1,newsize2);
    [adjimgs_resvect,adjcls_vect]=fill_closepix(adjimgs_resvect,adjcls_vect,fullimg_resvet,fullcls_resmat(:),newsize1,newsize2);
else
    adjcls_vect=fullcls_resmat(:);
    adjimgs_resvect=fullimg_resvet;
end
adjcls_mat=reshape(adjcls_vect,[newsize1,newsize2]);

% Step2: classify for current image by integrating temporal information from adjacent images
curimg_mat=reshape(curimg_vect,[size1,size2,size_band]);
curimg_resmat=imresize(curimg_mat,scale);
curimg_resvet=reshape(curimg_resmat,[newsize1*newsize2,size_band]);
curimg_resvet=[curimg_resvet,adjimgs_resvect];
curcls_resvect=kmeans(curimg_resvet,num_optcls,'MaxIter',1000);%
curcls_resmat=reshape(curcls_resvect,newsize1,newsize2);

% Step3: image spatial segementation by using complete adjacent images
adjimgs_resmat=reshape(adjimgs_resvect,[newsize1,newsize2,size(adjimgs_resvect,2)]);
adjgrey_complt=mean(adjimgs_resmat,3,'omitnan');
wind=fspecial('laplacian',0);
adjlapl_complt=imfilter(adjgrey_complt,wind,'replicate');
adjfilt_complt=adjgrey_complt-adjlapl_complt;% image enhance using Laplacian filter
[adjedge_complt,threshold] = edge(adjfilt_complt,'Canny');% edge extract
% fudgeFactor = 0.5;
% adjedge_complt = edge(adjgrey_complt,'Canny',threshold * fudgeFactor);
adjedge_complt(1,:)=1;
adjedge_complt(end,:)=1;
adjedge_complt(:,1)=1;
adjedge_complt(:,end)=1;
se90 = strel('line',3,90);
se0 = strel('line',3,0);
se45 = strel('line',3,45);
sedisk = strel('disk',1);
diledge_complt = imdilate(adjedge_complt,[se90 se0 se45 sedisk]);% line connection
dilwatsd_complt = watershed(diledge_complt);% water shed segmentation
dilLabl_complt = bwlabel(dilwatsd_complt);% segmentation objects
dilLabl_complt(dilLabl_complt==0)=nan;
dilLabl_complt = fillmissing(dilLabl_complt,'nearest');
% imshow(dilLabl_complt,[])
% dilLabl_complt = label2rgb(dilwatsd_complt);
% imshow(dilLabl_complt)

if sum(curcls_resvect>0)/(newsize1*newsize2)<0.95%0.2
    gapfill_cls=adjcls_mat;
else
    % Step4: genenrate the class image for missing gaps using adjacent class image as base and current class image and spatial segmentation image as adjust
    gapfill_cls=curcls_resmat;
    mis_mat=isnan(gapfill_cls);
    mis_num1=sum(mis_mat(:))+1;
    mis_num2=sum(mis_mat(:));
    while mis_num1>0 && mis_num2<mis_num1
        mis_num1=mis_num2;
        num_cls=unique(adjcls_vect(mis_mat(:)));
        for i=1:length(num_cls)
            adjcls_misi=adjcls_mat==num_cls(i) & mis_mat;
            Lablcls_i=dilLabl_complt(adjcls_misi);
            [C,~,~] = unique(Lablcls_i);% in case of cracks on results
            for j=1:length(C)
                adjcls_unmisij=adjcls_mat==num_cls(i) & ~mis_mat & dilLabl_complt==C(j);
                if sum(adjcls_unmisij,'all')>0
                    [Cj,~,Nj] =unique(curcls_resvect(adjcls_unmisij));
                    adjcls_misij=adjcls_mat==num_cls(i) & mis_mat & dilLabl_complt==C(j);
                    if length(Cj)>1
                        a_counts = accumarray(Nj,1)/sum(adjcls_unmisij(:));
                        [a_num,a_ind]=sort(a_counts,'descend');
                        [row_Cj,col_Cj] = find(adjcls_misij==1);
                        CjCenter=mean([row_Cj,col_Cj],1);
                        [row_2,col_2]=find(curcls_resmat.*adjcls_unmisij==Cj(a_ind(2)));
                        mean_second=mean(pdist2([row_2,col_2],CjCenter),'all','omitnan');
                        [row_1,col_1]=find(curcls_resmat.*adjcls_unmisij==Cj(a_ind(1)));
                        mean_first=mean(pdist2([row_1,col_1],CjCenter),'all','omitnan');
                        if a_num(2)>0.4 && mean_second<mean_first
                            curcls_iput=Cj(a_ind(2));
                        else
                            curcls_iput=Cj(a_ind(1));
                        end
                    else
                        curcls_iput=Cj(1);
                    end
                    gapfill_cls(adjcls_misij)=curcls_iput;
                end
            end
        end
        mis_mat=isnan(gapfill_cls);
        mis_num2=sum(mis_mat(:));
    end
    if mis_num2>0
        num_cls=unique(adjcls_vect(mis_mat(:)));
        for i=1:length(num_cls)
            adjcls_misi=adjcls_mat==num_cls(i) & mis_mat;
            if sum(adjcls_mat==num_cls(i) & ~mis_mat,'all')>0
                curcls_iput=mode(gapfill_cls(adjcls_mat==num_cls(i) & ~mis_mat),1);
                gapfill_cls(adjcls_misi)=curcls_iput;
            else
                num_clsall=unique(gapfill_cls(gapfill_cls>0));
                curcls_iput=num_clsall(end)+1;
                gapfill_cls(adjcls_misi)=curcls_iput;
            end
        end
    end
    % num_cls=unique(adjcls_vect(mis_mat(:)));
    % for i=1:length(num_cls)
    %     adjcls_i=adjcls_mat==num_cls(i) & mis_mat;
    %     Lablcls_i=dilLabl_complt(adjcls_i);
    %     [C,ia,ic] = unique(Lablcls_i);
    %     a_counts = accumarray(ic,1)/sum(adjcls_i(:));
    %     if max(a_counts(:))<0.7 && sum(adjcls_i(:))/(newsize1*newsize2)>0.02 % similarity<0.7
    %         for j=1:length(C)
    %             if sum(~mis_mat & dilLabl_complt==C(j),'all')>0 && a_counts(j)>0.1
    %                 curcls_iput=mode(curcls_resvect(~mis_mat & dilLabl_complt==C(j)),1);
    %                 adjcls_iput=mode(adjcls_vect(curcls_resmat==curcls_iput & dilLabl_complt~=C(j),:),1);
    %                 if adjcls_iput~=num_cls(i)
    %                     gapfill_cls(dilLabl_complt==C(j) & adjcls_mat==num_cls(i))=adjcls_iput;
    %                 end
    %             end
    %         end
    %     end
    % end
end
gapfill_cls=imresize(gapfill_cls,[size1,size2],'nearest');
img_open=imopen(gapfill_cls,se);
gapfill_cls=imclose(img_open,se);%remove tiny points
dilLabl_complt=imresize(dilLabl_complt,[size1,size2],'nearest');
% toc

function [curimg_vect,curcls_vect]=fill_closepix(curimg_vect,curcls_vect,aftimg_vect,aftcls_vect,size1,size2)
mis_mat=reshape(isnan(curcls_vect),size1,size2);
if sum(mis_mat(:))>0
    se = strel('square', 7);
    mis_dilmat=imdilate(mis_mat,se);
    num_cls=unique(aftcls_vect(mis_mat(:) & aftcls_vect>0));
    aftcls_dilvet=aftcls_vect(mis_dilmat & ~mis_mat);
    curcls_dilvet=curcls_vect(mis_dilmat & ~mis_mat);
    for i=1:length(num_cls)
        curcls_i=curcls_dilvet(aftcls_dilvet==num_cls(i));
        if isempty(curcls_i)
            curcls_i=curcls_dilvet;
        end
        curcls_iput=mode(curcls_i,1);
        curcls_vect(mis_mat(:)==1 & aftcls_vect==num_cls(i))=curcls_iput;
        
        curimg_dilvet=curimg_vect(mis_dilmat(:) & ~mis_mat(:) & curcls_vect==curcls_iput,:);
        aftimg_dilvet=aftimg_vect(mis_dilmat(:) & ~mis_mat(:) & aftcls_vect==num_cls(i),:);
        if isempty(aftimg_dilvet)
            aftcls_i=mode(aftcls_dilvet(curcls_dilvet==curcls_iput));
            aftimg_dilvet=aftimg_vect(mis_dilmat(:) & ~mis_mat(:) & aftcls_vect==aftcls_i,:);
        end
        test_y=mean(curimg_dilvet,1,'omitnan')-repmat(mean(aftimg_dilvet,1,'omitnan'),[1,size(curimg_dilvet,2)./size(aftimg_dilvet,2)])+repmat(aftimg_vect(mis_mat(:)==1 & aftcls_vect==num_cls(i),:),[1,size(curimg_dilvet,2)./size(aftimg_dilvet,2)]);
        curimg_vect(mis_mat(:)==1 & aftcls_vect==num_cls(i),:)=test_y;
    end
end