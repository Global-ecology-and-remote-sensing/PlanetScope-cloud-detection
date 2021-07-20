%% object-based cloud and cloud shadow matching
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%%
function [new_smask,new_cmask]=gene_newmask2(c_mask,s_mask,move_disx,move_disy,bcgd)
% translate cloud mask to shadow position, derive the moved cloud mask
mc_mask=imtranslate(double(c_mask),[move_disx, move_disy],'FillValues',255);% Transfering logical to double to separate background and target in imtranslate
mc_mask=mc_mask==1;
% calculate the moved backgroud including water and image edge
mbg_mask=imtranslate(bcgd,[move_disx, move_disy],'FillValues',1);% Transfering logical to double to separate background and target in imtranslate
mbg_mask=mbg_mask==1;
% moved cloud mask was separated into 2 parts; 1. overlap between moved cloud mask and original cloud mask
ovlp_cmc=mc_mask & c_mask;
% moved cloud mask was separated into 2 parts; 2. pixel can not be cloud and shadow at the same time
mc_mask=mc_mask & (~c_mask);
% calulate the matched cloud and shadow objects
[~,L_mcmask,~]=bwboundaries(mc_mask,'noholes');
mc_num=unique(L_mcmask);
for i=2:length(mc_num)% i=1 is the background
    mc_area=L_mcmask==mc_num(i); % for each moved cloud object, calculate the area
    ovlp_area=(s_mask | bcgd) & mc_area ; % calculate the overlap area between moved cloud object & shadow object %| mbg_mask)
    if sum(ovlp_area(:))/sum(mc_area(:))<0.2 || sum(mc_area(:))<=30 % if lower than the similarity threshold or tiny object
        L_mcmask(L_mcmask==mc_num(i))=0; % exclude the corresponding object
    end
end
% calculate the new shadow mask
[B_smask,L_smask,~]=bwboundaries(s_mask,'noholes');
% plot the shadow mask; similarly could be used to plot the cloud mask
% figure;imshow(label2rgb(L_smask, @jet, [.5 .5 .5]))
% hold on
% for k = 1:length(B_smask)
%    boundary = B_smask{k};
%    plot(boundary(:,2), boundary(:,1), 'w', 'LineWidth', 2)
% end
s_num=unique(L_smask(:));
new_smask=s_mask;
for i=2:length(s_num)
    L_ovlp=L_smask==s_num(i) & (L_mcmask>0 | mbg_mask); % calculate the overlap between shadow mask and matched mask object
    if sum(L_ovlp(:))<30 % If no overlap, exclude from the original shadow mask
        new_smask(L_smask==s_num(i))=0;
    end
end
% calculate the new cloud mask
mmc_mask=imtranslate(L_mcmask>0,[-move_disx, -move_disy],'FillValues',1);% move the matched mask back to cloud mask
movlpcmc_mask=imtranslate(ovlp_cmc,[-move_disx, -move_disy],'FillValues',1);% move the overlaped cloud mask back to cloud mask
mmc_mask=(mmc_mask | movlpcmc_mask | bcgd);% integrated mask used to determine if detected cloud object is right
[~,L_cmask,~]=bwboundaries(c_mask,'noholes');
c_num=unique(L_cmask(:));
new_cmask=c_mask;
for i=2:length(c_num)
    L_ovlp=L_cmask==c_num(i) & (mmc_mask);% calculate the overlap between cloud mask and matched mask object
    if sum(L_ovlp(:))<30% If no overlap, exclude from the original cloud mask
        new_cmask(L_cmask==c_num(i))=0;
    end
end
% figure;imshow(new_smask)
% figure;imshow(new_cmask)