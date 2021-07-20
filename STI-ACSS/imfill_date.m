%% Single-image-based cloud shadow detection using flood-fill approach
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%%
function sqrtB3_shade=imfill_date(sqrtB3_date,cloud)
[size1,size2]=size(sqrtB3_date);
sqrtB3_clear=sqrtB3_date.*(~cloud);
sqrtB3_clear(sqrtB3_clear==0)=nan;
if isnan(prctile(sqrtB3_clear(:),90))==0
    fillvalue=prctile(sqrtB3_clear(:),90);
else
    fillvalue=prctile(sqrtB3_date(:),90);
end
% expand image size 1/4~1/5 to reduce edge effect
expd1=round(size1/4);
expd2=round(size2/4);
sqrtB3_expd = f_expand(sqrtB3_date,expd1*2,expd2*2);
% out of land border
sqrtB3_expd(isnan(sqrtB3_expd))=fillvalue;%max(sqrtB3_expd,[],'all','omitnan');
% imfill image
sqrtB3_filled=imfill(sqrtB3_expd);
% calculate difference between image and filled image
sqrtB3_expddif=sqrtB3_filled-sqrtB3_expd;
sqrtB3_dif=sqrtB3_expddif(expd1+1:size1+expd1,expd2+1:size2+expd2);
sqrtB3_shade=sqrtB3_date<prctile(sqrtB3_date(sqrtB3_dif>0.1),50) | sqrtB3_dif>0.1;