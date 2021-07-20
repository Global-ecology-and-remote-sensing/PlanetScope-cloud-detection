%% calculate the turning point of the candidate thresholds-number of cloud pixels higher than candidate thresholds in 2-D plot
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%%
function [val, idx, stdDist]=thre_trangle(data)
if data(1)-data(end)>0 % calculate the cloud threshold
    slope=(data-data(end))./(size(data,1):-1:1)';
    [~,ind_max]=max(slope);% max slope make the algorithom more robust
    x = [ind_max size(data,1)];
    y = [data(ind_max) data(end)];
    offset1=ind_max-1;
    offset2=0;
    perpDist=cal_dist(data,x,y,offset1,offset2);
    % there are several scenarios to find the turning points
    % First two scenarios, when the ind_max is close to start point, but there are two turing points
    dif_Dist=diff(perpDist(:));
    turn_point=find(dif_Dist(1:end-1)>0 & dif_Dist(2:end) < 0);
    if length(turn_point)>1 && perpDist(turn_point(1))/perpDist(turn_point(2))>0.5 % there are two turning points, and the difference between these two distance is not very large
        idx=turn_point(1)+1;
        val=perpDist(idx);
    else % there is only one turning points, or the difference between these two distance is very large, turn to the point with largest distance
        [val, idx] = max(perpDist);
    end
    stdDist=std(perpDist(perpDist>0));
    
    % Second two scenarios, when the ind_max is far away from start point, need to find another turning point
    iter=0;
    while ind_max>20 && iter<3
        iter=iter+1;
        subslope=(data(1:ind_max)-data(ind_max))./(ind_max:-1:1)';
        [~,ind_submax]=max(subslope);
        x = [ind_submax ind_max];
        y = [data(ind_submax) data(ind_max)];
        offset1=ind_submax-1;
        offset2=ind_max-size(data,1);
        perpDist2=cal_dist(data,x,y,offset1,offset2);
        [val2, idx2] = max(perpDist2);
        stdDist2=std(perpDist2(perpDist2>0));
        if val2>val || (ind_max-ind_submax)/abs(offset2)>0.8 % if the distance large than fisrt two senarios
            val=val2;
            idx=idx2;
            stdDist=stdDist2;
            break;
        end
        ind_max=ind_submax;
    end
else % calculate the shadow threshold
    slope=(data-data(1))./(1:size(data,1))';
    [~,ind_max]=max(slope);
    x = [1 ind_max];
    y = [data(1) data(ind_max)];
    offset1=0;
    offset2=ind_max-size(data,1);
    perpDist=cal_dist(data,x,y,offset1,offset2);
    [val, idx] = max(perpDist);
    stdDist=std(perpDist(perpDist>0));
end
end
%% calculate the distance from points on the curve to the line
function perpDist=cal_dist(data,x,y,offset1,offset2)
% The slope of the line connecting the two endpoints
m = ( y(2) - y(1) )/( x(2) - x(1) );
pm= - 1 / m;
% Point on the curve (xc,yc), point on the line (xl,yl)
perpDist = zeros(size(data,1),1);
for i = offset1+1:size(data,1)+offset2
    xc = i ; yc = data(i);
    yl = ( (m * xc) + (m^2 * yc) - (m * x(1)) + y(1) )/(1+ m^2);
    xl = xc - m*(yl - yc);
    % distance^2
    d2 = (xl - xc)^2 + (yl - yc)^2;
    perpDist(i) = d2;
end
end