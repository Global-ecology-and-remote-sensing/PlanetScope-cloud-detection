%% calculate the projection direction from metadata file
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%% 
function [dirct_x,dirct_y]=calc_proj(file_path)
subdir  = dir(file_path);

for i = 1 : length( subdir )
    if( isequal( subdir( i ).name, '.' )==0 &&...
            isequal( subdir( i ).name, '..')==0 &&...
            subdir( i ).isdir)
        subdir_ind=i;
        break;
    end
end
xml_dir=dir([file_path,subdir(subdir_ind).name,'\*_metadata*.xml']);
fid_in=fopen(fullfile(file_path,subdir(subdir_ind).name,xml_dir(1).name),'r');
geo_char=fscanf(fid_in,'%c',inf);
fclose(fid_in);
geo_char=geo_char';
geo_str=strread(geo_char,'%s');

% Read in Solar Azimuth & Elevation angle (degrees)
char_AzA=char(geo_str(strmatch('<opt:illuminationAzimuthAngle',geo_str)+1));
azi=str2double(char_AzA(strfind(char_AzA,'"deg">')+6:strfind(char_AzA,'</')-1));
char_EA=char(geo_str(strmatch('<opt:illuminationElevationAngle',geo_str)+1));
zen=90-str2double(char_EA(strfind(char_EA,'"deg">')+6:strfind(char_EA,'</')-1));

dirct_x=-tand(zen)*sind(azi);
dirct_y=-tand(zen)*cosd(azi);
