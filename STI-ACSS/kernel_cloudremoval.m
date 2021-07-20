%% integration of single-image-based cloud and cloud shadow detection and Multi-temporal-based cloud and cloud shadow detection (task A-C)
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%%
function [shade,cloud]=kernel_cloudremoval(planet_original,size_bck,date)
[size_p1,size_p2,~,size4]=size(planet_original);
%% single-image-based cloud detection method; adaptive threshold based on HOT cloud index
% automatically calculate HOT parameters, slope and intercept, using ATSA method
HOT=HOT_slope(planet_original);
% HOT=planet_original(:,:,1,:)-0.5*planet_original(:,:,3,:); 
% adaptive threshold method
HOT_thre=HOT_thresh(HOT); 
%% HOT sensitive test
% HOT_thre=repmat(mean(HOT(:),'omitnan'),[1,size(date,1)]);
% HOT_thre=prctile(HOT,10,[1,2]);
% HOT_thre=prctile(HOT,30,[1,2]);
% HOT_thre=prctile(HOT,50,[1,2]);
% HOT_thre=prctile(HOT,70,[1,2]);
% HOT_thre=prctile(HOT,90,[1,2]);
% plot_thre(date,HOT_thre,HOT);
% cloud=HOT>reshape(HOT_thre,[1,1,1,size(date,1)]);
%% cloud detection 
cloud=HOT>HOT_thre;
mean_vis=planet_original(:,:,3,:)>prctile(planet_original(:,:,3,:),95,'all') & planet_original(:,:,4,:)>prctile(planet_original(:,:,4,:),95,'all');% brightness
% mean_vis=mean(planet_original(:,:,1:3,:),3,'omitnan')>3000;% brightness
cloud=cloud | mean_vis;% in case of some clouds not detected by HOT
%% single-image-based shadow detection method; floodfill approach based on a shadow index
sqrtB3=sqrt(planet_original(:,:,3,:)./mean(planet_original(:,:,3,:),[1,2],'omitnan')).*sqrt(planet_original(:,:,4,:)./mean(planet_original(:,:,4,:),[1,2],'omitnan'));% shadow index; pixels with low red and nir band are detected as shadow
sqrtB3_dif=zeros([size_p1,size_p2,1,size4]);
% flood fill and calculate the difference after and before floodfill
for z=1:size4
    sqrtB3_date=sqrtB3(:,:,1,z);
    sqrtB3_dif(:,:,1,z)=imfill_date(sqrtB3_date,cloud(:,:,1,z));
end
shade=sqrtB3_dif;%
% this case for detection cloud shadows in mountains areas where shadow coverage are much larger than cloud coverage due to mountain shadow effect
% shade=devi_shadowcover(planet_original,cloud);
% shade=shade & sqrtB3_dif;
clear HOT mean_vis sqrtB3_dif sqrtB3
%% multi-temporal-based cloud and shadow detection method;
cloud_mask=zeros(size_p1,size_p2,1,size4);
planet_outlier=zeros(size_p1,size_p2,1,size4);
planet_itr=planet_original;
% initial iteration setting
iter_bef=1000;
iter_aft=std(planet_itr,0,[1,2,4],'omitnan')./mean(planet_itr,[1,2,4],'omitnan');
itr=0;
thre_perc=5; % outlier threshold
while itr<10 && sum(iter_bef-iter_aft>0.01)>0 %iteration no more than 5 and std difference between i-1 and i less than a threshold
    % Centralization
    bond=prctile(planet_itr,[5,95],[1,2]);
    planet_mean_band=mean(planet_itr.*(planet_itr<bond(2,1,:,:) & planet_itr>bond(1,1,:,:)),[1,2],'omitnan');
%     planet_mean_band=mean(planet_itr,[1,2],'omitnan');
    planet_normli=planet_itr-planet_mean_band;
    % Blocking
    Imeans_date=sepblockfun(planet_normli,[size_bck,size_bck,1,size4],@nanmean);
    block_std=sqrt(sepblockfun(planet_normli.^2,[size_bck,size_bck,1,size4],@nanmean) - Imeans_date.^2);
    % select potential blocks with cloud and shadow
    block_ind=sum(block_std>mean(block_std,[1,2],'omitnan'),3,'omitnan')>0;
    
    block_date_ind=imdilate(block_ind,[0 1 0;1 1 1;0 1 0]);
    % detect outliers for each block
    for block1=1:size_p1/size_bck%50
        for block2=1:size_p2/size_bck%50
            planet_block=planet_normli((block1-1)*size_bck+1:block1*size_bck,(block2-1)*size_bck+1:block2*size_bck,:,:);
            planet_outlier((block1-1)*size_bck+1:block1*size_bck,(block2-1)*size_bck+1:block2*size_bck,1,:)=cloud_detect(planet_block,itr,thre_perc);%,block_ind_add(block1,block2)
        end
    end
    block_res_ind=imresize(block_date_ind,[size_p1,size_p2],'nearest');
    planet_outlier=planet_outlier.*repmat(block_res_ind,[1,1,1,size4]);%block_ind_add;
    % generate cloud and shadow mask
    cloud_mask=cloud_mask | planet_outlier;
    % prepare iteration
    planet_original(repmat(cloud_mask,[1,1,4,1]))=nan;
    planet_itr=planet_original;
    itr=itr+1
    iter_bef=iter_aft;
    iter_aft=std(planet_itr,0,[1,2,4],'omitnan')./mean(planet_itr,[1,2,4],'omitnan');
    clear bond planet_mean_band planet_normli
end
clear planet_itr planet_original
%% generate initial cloud and shadow masks
shade= cloud_mask & shade;
cloud= cloud_mask & cloud;