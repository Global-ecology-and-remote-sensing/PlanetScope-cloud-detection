%% --------------- Description ----------------------------------------------
% AUTOmatic SpatioTemporal Integration approach for Automatic Cloud and Shadow Screening (‘STI-ACSS’) 
%
% This code is used for automatic cloud and cloud shadow detection in
% PlanetScope time-series images. This version is not used for snow/ice/water/mountains
% detection yet.
%
% Please refer to :
%
% Jing Wang, Dedi Yang, Shuli Chen, Xiaolin Zhu, Shengbiao Wu, Marc Bogonovich, Zhengfei Guo, Zhe Zhu, Jin Wu*,
% "Automatic cloud and cloud shadow detection in tropical areas for PlanetScope satellite images", 
% Remote sensing of Enviroment, Volume 264, 2021, 112604, ISSN 0034-4257, https://doi.org/10.1016/j.rse.2021.112604. Accepted in 11/07/2021
% 
%% --------------- Usuage----------------------------------------------------
%     
% Function : autoSTIACSS
% 
% This main function includes three steps : 
% Step 1 : read PlanetScope time-series images one by one from the folder.
%          The time-series images should have same spatial extent (row * column * band).
% Step 2 : conduct STI-ACSS.
% Step 3 : write each cloud and cloud shadow mask in original folder of each PlanetScope time-series image.
%
% Input arguments :
% file_dir : the folder of PlanetScope time-series files, each file is named as "planet_order_*", e.g. "planet_order_20180108",
%            each file includes the to-be-detected PlanetScope image "*_mosaicnoqc.tif" and standard PlanetScope subfiles downloaded from https://www.planet.com/, 
%            each subfile includes "*_AnalyticMS_SR.tif","*_metadata.json","*_AnalyticMS_metadata.xml", and "*_udm/udm2.tif".
% mask_dir : the folder of water mask, which should have the same image
%            size with PlanetScope image. In this water mask: 1 for land surface, 0 for water body. 
% interval : Due to the MATLAB memory limitation, we read PlanetScope files with every n (i.e. interval) images. 
%            For example, an image time series generally includes <=30 images, 
%            if there are 80 images in a time series, the interval will be 3.
% upscale :  Due to the MATLAB memory limitation, we upscale PlanetScope
%            image with a upscale factor. For example, for a image with a spatial extent of 3334 * 3334 PlanetScope pixels (10km * 10km), the upscale is set as 1.5.
% bck_size : STI-ACSS divided PlanetScope image time series into discrete cubic blocks, the block size if ~500m * 500m.
%            for example, for a image with a spatial extent of 3334 * 3334 PlanetScope pixels (10km * 10km) and upscale the image with a upscale factor of 1.5, 
%            the block size is 100 (~500m/3m/1.5).
% 
% Output arguments :
% cloud and cloud shadow mask: 0 : shadow
%                              64 : background
%                              128 : clear 
%                              255 : cloud
% 
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%%                    
function autoSTIACSS
tic
file_dir ='E:\Wangjing\STI-ACSS\exampledata\2019\PSScene4Band\';% floder name of PlanetScope time series flles
planet_dir = dir([file_dir,'planet_order_*']);%planet data dir
dir_original = dir([file_dir,planet_dir(1).name,'\*_mosaicnoqc.tif']);%Planet data
planet_temp = imread([file_dir,planet_dir(1).name,'\',dir_original.name]);% read one image of this Planet time series
[size_d1,size_d2,size3]=size(planet_temp);
% mask_dir='E:\Wangjing\Download\BCI_enlarge\class\Planet\Planet_watermask_20181004_MP.tif';%file name of water mask
% Planet_mask=imread(mask_dir);
Planet_mask=ones(size_d1,size_d2);
interval=1;%time series interval
upscale=1.5;
bck_size=100;%determine the blocknum1,every 500m spatial extent
size_p1=round(size_d1/upscale/bck_size)*bck_size;% row of resized image size
size_p2=round(size_d2/upscale/bck_size)*bck_size;% col of resized image size

for inte=1:interval
    % Step 1: read PlanetScope time-series images one by one
    inte_len=floor(length(planet_dir)/interval)+(mod(length(planet_dir),interval)>=inte);
    planet_original=zeros(size_p1,size_p2,size3,inte_len);
    planet_date=strings(inte_len,1);
    date=zeros(inte_len,1);
    dirct_x=zeros(inte_len,1);
    dirct_y=zeros(inte_len,1);
    for i=inte:interval:length(planet_dir)
        dir_original = dir([file_dir,planet_dir(i).name,'\*_mosaicnoqc.tif']);%Planet data
        planet_temp= imread([file_dir,planet_dir(i).name,'\',dir_original.name]);% read Planet
        planet_temp= double(planet_temp).*double(repmat(Planet_mask,[1,1,size3]));
        planet_original(:,:,:,ceil(i/interval))=imresize(planet_temp,[size_p1,size_p2],'nearest');
        planet_date(ceil(i/interval))=dir_original.name(1:8);
        date(ceil(i/interval))=datenum(dir_original.name(1:8),'yyyymmdd')-datenum([dir_original.name(1:4),'0101'],'yyyymmdd')+1;
        [dirct_x(ceil(i/interval),:),dirct_y(ceil(i/interval),:)]=calc_proj([file_dir,planet_dir(i).name,'\']);
    end
    planet_date=char(planet_date);
    planet_bg = sum(planet_original,3)==0; %background
    planet_original(planet_original<=0)=nan;
    % Step 2: STI-ACSS process
    [shade,cloud]=kernel_cloudremoval(planet_original,bck_size,date);
    % Step 3: write PlanetScope time-series images one by one
    for i=inte:interval:length(planet_dir)%
        bg_cuur=planet_bg(:,:,:,ceil(i/interval));
        shade_cuur=shade(:,:,:,ceil(i/interval));
        cloud_cuur=cloud(:,:,:,ceil(i/interval));
        % write the cloud and cloud shadow mask of inter-step 1 (Task A-C)
        csmask_cuur = ones(size_p1,size_p2,'uint8')*128;% clear  
        csmask_cuur(bg_cuur)=64;% background         
        csmask_cuur(cloud_cuur)=255;% cloud
        csmask_cuur(shade_cuur)=0;% shadow
        imwrite(csmask_cuur,[file_dir,planet_dir(i).name,'\',planet_date(ceil(i/interval),:),'_csinterstep1.tif'],'tif');
        
        [real_smask,s_mask]=morphlg_cr(shade_cuur,bg_cuur,cloud(:,:,:,ceil(i/interval)),shade(:,:,:,ceil(i/interval)));%shade
        [real_cmask,c_mask]=morphlg_cr(cloud_cuur,bg_cuur,cloud(:,:,:,ceil(i/interval)),shade(:,:,:,ceil(i/interval)));%cloud
        % write the cloud and cloud shadow mask of inter-step 2 (Task A-D)
%         csmask_cuur2 = ones(size_p1,size_p2,'uint8')*128;% clear  
%         csmask_cuur2(csmask_cuur==64)=64;%background
%         csmask_cuur2(c_mask>0)=255;% cloud
%         csmask_cuur2(s_mask>0)=0;% shadow
%         imwrite(csmask_cuur2,[file_dir,planet_dir(i).name,'\',planet_date(ceil(i/interval),:),'_csinterstep2.tif'],'tif');
        
        [new_cmask,new_smask]=match_cloudshade_res(c_mask,s_mask,bg_cuur,dirct_x(ceil(i/interval)),dirct_y(ceil(i/interval))); 
        new_smask=imresize(new_smask,[size_d1,size_d2],'nearest');
        new_cmask=imresize(new_cmask,[size_d1,size_d2],'nearest');
        bg_cuur=imresize(bg_cuur,[size_d1,size_d2],'nearest');
        new_smask=bwareaopen(new_smask,200);
        new_cmask=bwareaopen(new_cmask,200);
        % write the final cloud and cloud shadow mask (Task A-E)
        cs_mask = ones(size_d1,size_d2,'uint8')*128;% clear 
        cs_mask(bg_cuur)=64;% background         
        cs_mask(new_cmask>0)=255;% cloud
        cs_mask(new_smask>0)=0;% shadow
        imwrite(cs_mask,[file_dir,planet_dir(i).name,'\',planet_date(ceil(i/interval),:),'_noqccloudseries1mask.tif'],'tif');
        i
    end
    clear planet_original planet_itr planet_normli planet_mean_band cloud_mask planet_outlier cs_final
end
toc

