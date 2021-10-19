%% generate mosaicked PlanetScope time-series images without using default quality control layer
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 15/10/2021
% -------------------------------------------------------
% Input arguments : 
% file_dir: the folder of PlanetScope time-series files, each file is named as "planet_order_*", e.g. "planet_order_20180108",
%           each file includes the clipped images named as "*_clipnoqc.tif"
% ROI region: either using the latitude and longtitude of the upleft and lowright points 
%             or using the latitude and longtitude of the center point and the radius (m) of the ROI region.
%
% Output arguments :
% Write the mosaicked image in the original folder for each PlanetScope time-series image.         
                    
%% --------------------------------------------------------------------------
function Planet_noqc_mosaic % Step2
file_dir='G:\wangjing\PlanetFusionData\Planet_deciduous\Congo\2018\PSScene4Band\';% replace it with your file dir
planet_dir=dir(strcat(file_dir,'planet_order_*'));%dir
sub_dir1=dir(strcat(file_dir,planet_dir(1).name,'\*_clipnoqc.tif'));%sub dir
[planet_image1,sub_info1]=geotiffread(strcat(file_dir,planet_dir(1).name,'\',sub_dir1(1).name));
geoinfo = geotiffinfo(strcat(file_dir,planet_dir(1).name,'\',sub_dir1(1).name));
geoTags = geoinfo.GeoTIFFTags.GeoKeyDirectoryTag;
[~,~,dim]=size(planet_image1);
% define the ROI region for clipping
% way 1:
center=[7534797 319605]; %Alice spring
input_ROI=[center(1)+5001,center(2)-5001,center(1)-5001,center(2)+5001];%
% way 2:
% input_ROI=[9687287,720655,9682286,730657];% k67 site
% mosaic information
mosaic_info=sub_info1;
mosaic_info.XWorldLimits=round([input_ROI(2),input_ROI(4)]);
mosaic_info.YWorldLimits=round([input_ROI(3),input_ROI(1)]);
RasterExtentInWorldX=mosaic_info.XWorldLimits(2)-mosaic_info.XWorldLimits(1);
RasterExtentInWorldY=mosaic_info.YWorldLimits(2)-mosaic_info.YWorldLimits(1);
mosaic_info.RasterSize = [RasterExtentInWorldY/round(sub_info1.CellExtentInWorldY),RasterExtentInWorldX/round(sub_info1.CellExtentInWorldX)];
for i=1:length(planet_dir)
    sub_dir=dir(strcat(file_dir,planet_dir(i).name,'\*_clipnoqc.tif'));%sub dir
    mosaic_image=zeros(mosaic_info.RasterSize(1),mosaic_info.RasterSize(2),dim,length(sub_dir));
    if isempty(sub_dir)
        disp(['No data for ',planet_dir(i).name]);
        continue;%if no SR data,it will not be used
    end
    % read clipped images for mosaicking 
    for sub_i=1:length(sub_dir) 
        [planet_image,sub_info]=geotiffread(strcat(file_dir,planet_dir(i).name,'\',sub_dir(sub_i).name));
        info = geotiffinfo(strcat(file_dir,planet_dir(i).name,'\',sub_dir(sub_i).name));
        if info.Zone==geoinfo.Zone                     
            %change to row,col coordinate
            [row_1_lu,col_1_lu]=map2pix(mosaic_info,sub_info.XWorldLimits(1),sub_info.YWorldLimits(2));
            [row_1_rd,col_1_rd]=map2pix(mosaic_info,sub_info.XWorldLimits(2),sub_info.YWorldLimits(1));
            
            row_1_lu=ceil(row_1_lu);
            col_1_lu=ceil(col_1_lu);
            row_1_rd=floor(row_1_rd);
            col_1_rd=floor(col_1_rd);%change to integer
            
            mosaic_image(row_1_lu:row_1_rd,col_1_lu:col_1_rd,:,sub_i)=planet_image;
        else % in case it has different zone with the first clipped images
            [XX,YY]=meshgrid(1:size(planet_image,2),1:size(planet_image,1));            
            Y=sub_info.YWorldLimits(2)-3*(YY-1);
            X=sub_info.XWorldLimits(1)+3*(XX-1);
            [lat,lon] = projinv(info,X,Y);
            
            [input_X,input_Y] = projfwd(geoinfo,lat,lon);
            %change to row,col coordinate
            [row,col]=map2pix(mosaic_info,input_X,input_Y);
            row=round(row);
            col=round(col);
            ind_valid=row>=1 & row<=mosaic_info.RasterSize(1) & col>=1 &col<=mosaic_info.RasterSize(2);
            row_new=row(ind_valid);
            col_new=col(ind_valid);
            planet_temp=reshape(planet_image,[size(planet_image,1)*size(planet_image,2),dim]);
            ind_new=sub2ind(mosaic_info.RasterSize,row_new,col_new);
            mosaic_temp=zeros(mosaic_info.RasterSize(1)*mosaic_info.RasterSize(2),dim);
            mosaic_temp(ind_new,:)=planet_temp(ind_valid,:);
            mosaic_image(:,:,:,sub_i)=reshape(mosaic_temp,[mosaic_info.RasterSize(1),mosaic_info.RasterSize(2),dim]);
        end
    end
%     histogram matching
    mosaic_new=zeros(mosaic_info.RasterSize(1),mosaic_info.RasterSize(2),dim,2);
    mosaic_new(:,:,:,1)=mosaic_image(:,:,:,1);
    for sub_i=2:length(sub_dir) 
        mosaic_new(:,:,:,2)=mosaic_image(:,:,:,sub_i);        
        mosaic_new(mosaic_new==0)=nan;
        ref_img=squeeze(mosaic_new(:,:,:,1));
        targ_img=squeeze(mosaic_new(:,:,:,2));
        for band=1:dim
            overlap=ref_img(:,:,band)>0 & targ_img(:,:,band)>0;
            if sum(overlap(:)>0)
                input_band=targ_img(:,:,band);
                [mu_input,sigma_input] = normfit(input_band(overlap==1));
                ref_band=ref_img(:,:,band);
                [mu_ref,sigma_ref] = normfit(ref_band(overlap==1));
                param_a=sigma_ref/sigma_input;
                param_b=mu_ref-param_a*mu_input;
                targ_img(:,:,band)=input_band.*param_a+param_b;
            end
        end
        mosaic_new(:,:,:,2)=targ_img;
        mosaic_new(:,:,:,1)=mean(mosaic_new,4,'omitnan');
        mosaic_new(:,:,:,2)=zeros(mosaic_info.RasterSize(1),mosaic_info.RasterSize(2),dim);
     end
%     mosaic image with average
    mosaic_image_comb=mosaic_new(:,:,:,1);% average in overlap
    mosaic_image_comb(isnan(mosaic_image_comb))=0;
    if sum(mosaic_image_comb(:,:,1)>0,[1,2])/(mosaic_info.RasterSize(1)*mosaic_info.RasterSize(2))<0.2
        disp(['No enough data for ',planet_dir(i).name]);
    end   
    geotiffwrite(strcat(file_dir,planet_dir(i).name,'\',sub_dir(1).name(1:8),'_mosaicnoqc.tif'),mosaic_image_comb,mosaic_info,'GeoKeyDirectoryTag',geoTags);    
    clear mosaic_image_comb mosaic_band mosaic_image
end