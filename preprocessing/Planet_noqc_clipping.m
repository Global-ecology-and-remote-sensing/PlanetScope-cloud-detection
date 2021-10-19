%% generate clipped PlanetScope time-series images without using default quality control layer
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 15/10/2021
% -------------------------------------------------------
% Input arguments : 
% file_dir: the folder of PlanetScope time-series files, each file is named as "planet_order_*", e.g. "planet_order_20180108",
%           each file includes standard PlanetScope subfiles downloaded from https://www.planet.com/, 
%           each subfile includes "*_AnalyticMS_SR.tif".
% ROI region: either using the latitude and longtitude of the upleft and lowright points 
%             or using the latitude and longtitude of the center point and the radius (m) of the ROI region.
%
% Output arguments :
% Write the clipped image in the original folder for each PlanetScope time-series image.         
                    
%% --------------------------------------------------------------------------
function Planet_noqc_clipping % Step1
file_dir='G:\wangjing\PlanetFusionData\Planet_deciduous\Congo\2018\PSScene4Band\';% replace it with your file dir
planet_dir=dir(strcat(file_dir,'planet_order_*'));% file dir
% define the ROI region for clipping
% way 1:
center=[7534797 319605]; %Alice spring
input_ROI=[center(1)+5001,center(2)-5001,center(1)-5001,center(2)+5001];%
% way 2:
% input_ROI=[9687287,720655,9682286,730657];% k67 site
for i=1:length(planet_dir)
    sub_dir1=dir([file_dir,planet_dir(i).name,'\201*']);%sub dir
    isub = [sub_dir1(:).isdir];
    sub_dir = sub_dir1(isub);
    for sub_i=1:length(sub_dir) 
        image_dir=dir(strcat(file_dir,planet_dir(i).name,'\',sub_dir(sub_i).name,'\','*_AnalyticMS_SR.tif'));        
        if isempty(image_dir)
            disp(['No AnalyticMS_SR.tif for ',planet_dir(i).name,'\',sub_dir(sub_i).name]);
            continue;%if no SR data,it will not be used
        end
        [planet_image,geo_info]=geotiffread(strcat(file_dir,planet_dir(i).name,'\',sub_dir(sub_i).name,'\',image_dir.name));%planet surfacere flectance image
        info = geotiffinfo(strcat(file_dir,planet_dir(i).name,'\',sub_dir(sub_i).name,'\',image_dir.name));
        geoTags = info.GeoTIFFTags.GeoKeyDirectoryTag;
        [size1,size2,dim]=size(planet_image);
        % clip with the identified ROI region
        ROI_lu_lat=min(input_ROI(1),geo_info.YWorldLimits(2));
        ROI_lu_lon=max(input_ROI(2),geo_info.XWorldLimits(1));
        ROI_rd_lat=max(input_ROI(3),geo_info.YWorldLimits(1));
        ROI_rd_lon=min(input_ROI(4),geo_info.XWorldLimits(2));
        
        [ROI_lu_row,ROI_lu_col]=map2pix(geo_info,ROI_lu_lon,ROI_lu_lat);% change to row,col coordinate
        [ROI_rd_row,ROI_rd_col]=map2pix(geo_info,ROI_rd_lon,ROI_rd_lat);% change to row,col coordinate
        ROI_lu_row=ceil(ROI_lu_row);
        ROI_lu_col=ceil(ROI_lu_col);
        ROI_rd_row=floor(ROI_rd_row);
        ROI_rd_col=floor(ROI_rd_col);
        if ROI_lu_row<1 || ROI_lu_row>size1 || ROI_lu_col<1 || ROI_lu_col>size2 || ROI_rd_row<1 || ROI_rd_row>size1 || ROI_rd_col<1 || ROI_rd_col>size2
            continue;
        end% without overlap with your ROI
        crop_image=planet_image(ROI_lu_row:ROI_rd_row,ROI_lu_col:ROI_rd_col,:);% crop image
        % geoinfomation of the clipped image
        sub_info=geo_info;
        sub_info.RasterSize = size(crop_image);
        sub_info.XWorldLimits(1)=geo_info.XWorldLimits(1)+(ROI_lu_col-1)*geo_info.CellExtentInWorldX;
        sub_info.XWorldLimits(2)=geo_info.XWorldLimits(1)+ROI_rd_col*geo_info.CellExtentInWorldX;
        sub_info.YWorldLimits(2)=geo_info.YWorldLimits(2)-(ROI_lu_row-1)*geo_info.CellExtentInWorldY;
        sub_info.YWorldLimits(1)=geo_info.YWorldLimits(2)-ROI_rd_row*geo_info.CellExtentInWorldY;%correct the bug, not need to take attention for input_ROI
        % write the clipped image
        geotiffwrite(strcat(file_dir,planet_dir(i).name,'\',sub_dir(sub_i).name,'_clipnoqc.tif'),crop_image,sub_info,'GeoKeyDirectoryTag',geoTags);
        clear planet_image quality_control quality_control_dim
    end
end
 