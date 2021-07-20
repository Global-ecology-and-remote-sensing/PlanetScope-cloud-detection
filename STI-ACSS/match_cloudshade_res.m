%% object-based cloud and cloud shadow matching
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%%
function [new_cmask,new_smask]=match_cloudshade_res(c_mask,s_mask,bcgd,varargin)
% avoid the tiny object affecting the matching results
c_mask_center = bwareaopen(c_mask,80);% seieve tiny objects
s_mask_center = bwareaopen(s_mask,80);% seieve tiny objects
if sum(s_mask_center(:))==0 || sum(c_mask_center(:))==0 % if on cloud or shadow objects after seieving, keep original cloud or shadow masks
    s_mask_center=s_mask;
    c_mask_center=c_mask;
end
% rough matching
interval1=6;
dirct_x=sign(varargin{1});
dirct_y=-sign(varargin{2});
tanyx=abs(varargin{2}/varargin{1});
startx1=min(600*dirct_x,0);
endx1=max(600*dirct_x,0);
% calculate the edge of shadow mask; as the edge of shadow affect the accuracy of cloud and shadow matching
%% Two methods calculate the shadow edge; first, image shrink to small size
% bcgd_resz=imresize(bcgd,3/5);
% [resz_1,resz_2]=size(bcgd_resz);
% [size_1,size_2]=size(c_mask_center);
% del_edge=ones(size_1,size_2);
% left=floor((size_1-resz_1)/2)+1;
% right=size_1-floor((size_1-resz_1)/2);
% up=floor((size_2-resz_2)/2)+1;
% down=size_2-floor((size_2-resz_2)/2);
% del_edge(left:right,up:down)=bcgd_resz;
%% Two methods calculate the shadow edge; second, image move along the predicted direction
bcgd_edge = imtranslate(bcgd,[round(dirct_x*100), round(tanyx*dirct_y*100)],'FillValues',1);
% Traverse all possible distance along the projected direction
overlap_num=zeros((endx1-startx1)/interval1+1,1);
for i=startx1:interval1:endx1
    j=dirct_y*round(tanyx*abs(i));
    temp_mask = imtranslate(c_mask_center ,[i, j]);
    overlap_mask=temp_mask & (s_mask_center & (~bcgd_edge));
    overlap_num((i-startx1)/interval1+1)=sum(overlap_mask(:));
end

% select the maximum overlapping area
[~,move_ind]=max(overlap_num(:));
[move_i,move_j]=ind2sub(size(overlap_num),move_ind);
% calculate the moved cloud mask
move_cmask = imtranslate(c_mask_center,[startx1+(move_i-1)*interval1, dirct_y*round(tanyx*abs(startx1+(move_i-1)*interval1))]);
move_cmask = move_cmask & (~c_mask_center); % exclude overlap with original cloud mask
% calculate the overlapped moved cloud and shadow mask and matching area
ovlp_mask = move_cmask & s_mask_center;
max_num=sum(ovlp_mask(:));
% calculate the cloud percentage
c_scale=sum(c_mask(:))/sum(bcgd==0,'all');
% calculate the ratio of cloud to shadow percentage
cs_scale=sum(c_mask(:))/sum(s_mask(:));
% 7 senarioes grouped into 3 senatioes, namely clear, thin cloud and matched
% case 1 exclude all; case 2 keep all; case 3 matched
if max_num/sum(move_cmask(:)+eps)<0.2 || max_num/sum(s_mask_center(:)+eps)<0.2 % Overall matching similarity lower than threshold, cloud and shadow layers maybe not well be matched
%     disp('Cloud and shadow layers maybe not well be matched!');
    if move_i==1 || move_i>=size(overlap_num,1)-2 || move_j==1 || move_j>=size(overlap_num,2)-2 || move_j==(0-startx1)/interval1+1 || move_i==(0-starty1)/interval1+1% wrong matching to the edge of searching range and original position (0,0)
        if c_scale<0.011 && max_num/(sum(move_cmask(:))+eps)>=0.28 && max_num/sum(s_mask_center(:)+eps)>0.01
            disp('case3');
            status=3; % matched
        elseif c_scale<0.011  %% few cloud; less than 1%
            disp('case1');
            status=1; % clear; exclude all
        else
            disp('case2');
            status=2; % thin cloud; keep all
        end
    elseif c_scale<0.011 && max_num/(sum(move_cmask(:))+eps)<0.3% few cloud; less than 1%
        disp('case1');
        status=1; % clear; exclude all
    elseif c_scale<0.011 && max_num/(sum(move_cmask(:))+eps)>=0.3
        disp('case3');
        status=3; % matched
    elseif c_scale>0.09 && cs_scale>2 % many thin clouds and cloud much larger than shadow
        disp('case2');
        status=2;% thin cloud; keep all
    else
        disp('case3');
        status=3;% matched
    end
else
    if c_scale>0.05 && cs_scale>2 && max_num/(sum(move_cmask(:))+eps)<0.5 % many thin clouds and cloud much larger than shadow and mathced cloud percentage is low
        disp('case2');
        status=2;%many thin clouds
    else
        disp('case3');
        status=3;%matched
    end
end
switch status
    case 1 % clear
        new_cmask=zeros(size(c_mask));
        new_smask=zeros(size(s_mask));
    case 2 % many clouds
        new_cmask=c_mask;
        new_smask=s_mask;
    case 3 % matched
        clear overlap_num
        % fine matching
        interval2=2;
        startx2=startx1+(move_i-1)*interval1-interval1;
        endx2=startx1+(move_i-1)*interval1+interval1;
        overlap_num=zeros((endx2-startx2)/interval2+1,1);
        for i=startx2:interval2:endx2
            j=dirct_y*round(tanyx*abs(i));
            temp_mask = imtranslate(c_mask_center ,[i, j]);
            overlap_mask=temp_mask & (s_mask_center & (~bcgd_edge));
            overlap_num((i-startx2)/interval2+1)=sum(overlap_mask(:));           
        end
        [~,move_ind]=max(overlap_num(:));
        [move_ii,move_jj]=ind2sub(size(overlap_num),move_ind);
        move_disx=startx2+(move_ii-1)*interval2;
        move_disy=dirct_y*round(tanyx*abs(startx2+(move_ii-1)*interval2));
        % exclude unmatched cloud and shadow objects
        [new_smask,new_cmask]=gene_newmask2(c_mask,s_mask,move_disx,move_disy,bcgd);
end