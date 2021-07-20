%% shadow coverage estimation using cloud coverage
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%%
function shade=devi_shadowcover(planet_original,cloud)
[size1,size2,~,size4]=size(planet_original);
shade=zeros(size1,size2,1,size4);
cloud_cover=sum(cloud,[1,2])./sum(planet_original(:,:,1,:)>0,[1,2]);
shadow_cover=min(1-cloud_cover,cloud_cover);
band_mean1=mean(planet_original,[1,2],'omitnan');
planet_mean1=mean(planet_original,[1,2,4],'omitnan');
thick_cloud=sum(band_mean1>planet_mean1,3)>0;
shadow_cover(thick_cloud & shadow_cover<0.1)=0.2;
shadow_cover(shadow_cover<0.1)=0.1;
for z=1:size4
    shadow_threshold1=prctile(planet_original(:,:,4,z),shadow_cover(z)*100,[1,2]);
%     shadow_threshold2=prctile(sqrtB3(:,:,1,z),shadow_cover(z)*100,[1,2]);
    shade(:,:,1,z)=planet_original(:,:,4,z)<shadow_threshold1;% & sqrtB3(:,:,1,z)<shadow_threshold2;
end