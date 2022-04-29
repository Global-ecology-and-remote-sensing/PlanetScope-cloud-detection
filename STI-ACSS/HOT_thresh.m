%% Adaptive approach to calculate the HOT threshold 
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%%
function threshold=HOT_thresh(HOT_img)
[~,~,~,size4]=size(HOT_img);
threshold=zeros(1,1,1,size4);
interval=10;% interval of threshold interation
startpoint=round(prctile(HOT_img(:),2.5)/interval)*interval;% start point of threshold interation
if std(HOT_img(:),[],'omitnan')/(mean(HOT_img(:),'omitnan')+eps)<1
    endpoint=round(prctile(HOT_img(:),99)/interval)*interval;% end point of threshold interation
else
    endpoint=round(prctile(HOT_img(:),97.5)/interval)*interval;% 
end
x=startpoint:interval:endpoint;
x=reshape(x,[1,1,length(x)]);
% Ymax=prctile(HOT_img,97.5,[1,2]);
% Ymin=prctile(HOT_img,2.5,[1,2]);
% k=(endpoint-startpoint)/(Ymax-Ymin);
% HOT_img=startpoint+k.*(HOT_img-Ymin);
std_dis=zeros(1,1,1,size4);
for i=1:size4
    rep_img=repmat(HOT_img(:,:,1,i),[1,1,length(x)]);
    thre_img=rep_img>x;
    thre_num=sum(thre_img,[1,2]);
%     figure;plot(squeeze(x)',thre_num(:))
%     xlabel('Threshold')
%     ylabel('Numbers of cloud pixels');%ylabel('Shade threshold');
%     set(gca,'FontSize',15,'FontWeight','bold','Linewidth',2);
%     dx = gradient(thre_num(:));
%     dxx = gradient(dx(:));
%     figure;plot(dxx(:))
%     [~,ind_max]=max(dxx);
    [~, ind_max,std_dis(1,1,1,i)]=thre_trangle(thre_num(:));
    threshold(1,1,1,i)=x(1)+interval*(ind_max-1);
%     text(400,0,num2str(threshold(1,1,1,i)))
%     saveas(gcf,['HOT',num2str(i),'.fig'],'fig')
    clear rep_img thre_img
end
% nan_ind=std_dis./max(std_dis(:))<0.2 & threshold>prctile(threshold(:),80);%;10 %normalize the distance to distance/x asis value, namely sin value
% threshold(nan_ind)=mean(threshold(~nan_ind),'omitnan');% faild threshold will be set the minimum value of all thresholds in this time series
clear rep_img thre_img
