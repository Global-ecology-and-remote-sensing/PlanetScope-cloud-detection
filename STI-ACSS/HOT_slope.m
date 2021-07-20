%% Calculate the HOT index following Zhu and Helmer (2018)
% More details please refer to:
%       Zhu, X.L., & Helmer, E.H. (2018). An automatic method for screening clouds and cloud shadows in optical satellite image time series in cloudy regions. 
%       Remote Sensing of Environment, 214, 135-153.
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%%
function HOT=HOT_slope(planet_original)%,dehaze]
[size1,size2,size3,size4]=size(planet_original);
slope=zeros(1,1,1,size4);
intercept=zeros(1,1,1,size4);
coef=zeros(1,1,1,size4);
% dehaze=zeros(size1,size2,size3,size4);
nbins=50;
dn_max=10000;
rmax=0.15*dn_max;
rmin0=0.01*dn_max;
for z=1:size4
    blue_temp=planet_original(:,:,1,z)<=rmax & planet_original(:,:,1,z)>rmin0;
    blue_temp=imerode(blue_temp,[0 1 0;1 1 1;0 1 0]);% seive tiny points
    planet_temp=reshape(planet_original(:,:,:,z),[size1*size2,size3]);
    clear_val=planet_temp(blue_temp==1,[1,3]);% extract blue bands and red bands
    valstart=min(clear_val(:,1));
    valend=max(clear_val(:,1));% start and end value of blue band
    
    interval=(valend-valstart)/nbins;
    clear_bin=zeros(nbins,2);% 50 bins
    if size(clear_val,1)>500
        for i=1:nbins
            inte_start=valstart+(i-1)*interval;
            inte_end=valstart+i*interval;% start and end value of bins
            id_inte=clear_val(:,1)<inte_end & clear_val(:,1)>=inte_start;
            if sum(id_inte)>=20
                clear_inte=clear_val(id_inte,:);
                ind_clear=clear_inte(:,2)>mean(clear_inte(:,2))+3*std(clear_inte(:,2));
                clear_inte(ind_clear,:)=nan;
                top_num=min(20,ceil(0.01*sum(clear_inte(:,2)>0)));
                [~,id_max]=maxk(clear_inte(:,2),top_num);% maxest 20 values of red bands in this bin
                clear_max=clear_inte(id_max,:);% selected clear pixels in binds
                clear_bin(i,:)=mean(clear_max,1,'omitnan');
                clear id_max clear_max
            end
        end
    end
    if sum(clear_bin(:,1)>0)>=0.5*nbins
        band1=clear_bin(:,1);
        band3=clear_bin(:,2);
        ind_valid=band1>0 & band3>0;
%         plot(band1(ind_valid),band3(ind_valid));
        
        coef_mat=corr(band1(ind_valid),band3(ind_valid));
        coef(1,1,1,z)=coef_mat(1);
        p = polyfit(band1(ind_valid),band3(ind_valid),1);
        slope(1,1,1,z)=p(1);
        intercept(1,1,1,z)=p(2);
        if p(1)<1.5
            slope(1,1,1,z)=0;
            intercept(1,1,1,z)=0;
        end
    else
        slope(1,1,1,z)=0;
        intercept(1,1,1,z)=0;
    end
end
if sum(slope==0)<size4
    intercept(slope==0)=mean(intercept(slope>0));
    slope(slope==0)=mean(slope(slope>0));
else
    slope(slope==0)=2;
end
HOT=(planet_original(:,:,1,:).*repmat(slope,[size1,size2,1,1])-planet_original(:,:,3,:)+repmat(intercept,[size1,size2,1,1]))./sqrt(1+repmat(slope,[size1,size2,1,1]).^2);
