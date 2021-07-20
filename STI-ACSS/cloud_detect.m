%% Multi-temporal-based cloud and cloud shadow detection using a pair of outlier threshold (i.e. thre_perc)
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%%
function cloud_mask=cloud_detect(planet_original,itr,thre_perc)
[size_p1,size_p2,size3,size4]=size(planet_original);
ind_conv=0;%itr*0.5;

planet_percentile = prctile(planet_original,[thre_perc/2-ind_conv 100-thre_perc/2+ind_conv],[1,2,4]);
Imeans_block=mean(planet_original,[1,2,4],'omitnan');
Istd_block=std(planet_original,[],[1,2,4],'omitnan');
Id_cloudy=sum(planet_original>Imeans_block+2*Istd_block,'all')>0;
Id_shade=sum(planet_original<Imeans_block-2*Istd_block,'all')>0;
if Id_cloudy && Id_shade
    planet_outlier=sum(planet_original<repmat(planet_percentile(1,1,:),[size_p1,size_p2,1,size4]) | planet_original>repmat(planet_percentile(2,1,:),[size_p1,size_p2,1,size4]),3,'omitnan');
elseif Id_cloudy
    planet_outlier=sum(planet_original>repmat(planet_percentile(2,1,:),[size_p1,size_p2,1,size4]),3,'omitnan');
elseif Id_shade
    planet_outlier=sum(planet_original<repmat(planet_percentile(1,1,:),[size_p1,size_p2,1,size4]),3,'omitnan');
else
    planet_outlier=zeros(size_p1,size_p2,1,size4);
end
cloud_mask=planet_outlier>0;