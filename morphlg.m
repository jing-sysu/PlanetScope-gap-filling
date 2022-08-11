function new_img=morphlg(img,size1,size2,wind)
img=reshape(img,[size1,size2]);

se=strel('square',wind);%Morphological structuring element
img_open=imopen(img,se);
img_close=imclose(img_open,se);%remove tiny points
img_bwopen=bwareaopen(img_close,250);% cloudoc=bwareaopen(cloudoc,250);
img_dilat=imdilate(img_bwopen,se);
new_img=img_dilat(:);