%% Fine-tuning on the cloud and shadow masks using morphological processing 
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%%
function [realoc,mcmask]=morphlg_cr(csmask,img,haze,shade)
% img_band=img(:,:,1);
se1=strel('square',7);%[0 1 0;1 1 1;0 1 0];%%Morphological structuring element% se1=strel('square',5);
cloudo=imopen(csmask,se1);
se2=strel('square',7);
cloudoc=imclose(cloudo,se2);%remove tiny points
cloudoc=bwareaopen(cloudoc,25);% cloudoc=bwareaopen(cloudoc,250);
cloudoc=imdilate(cloudoc,se1);
% cloudoc = edgegrow(img,cloudoc,haze,shade);
rmask=cloudoc==0 & img==0;%real image mask after cloud detection
realo=imopen(rmask,se1);% se2=strel('square',7);%
realoc=imclose(realo,se2);%remove tiny points
realoc = bwareaopen(realoc,36);% realoc = bwareaopen(realoc,360);
% realoc=seieve(planet_curr,realoc,haze(:,:,:,ceil(i/interval)),shade(:,:,:,ceil(i/interval)));
mcmask=~realoc & img==0;%final cloud mask by PCA