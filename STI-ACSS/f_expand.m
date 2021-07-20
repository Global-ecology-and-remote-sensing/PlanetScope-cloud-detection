%% expand the matrix spatial extent by mirror symmetry
% -------------------------------------------------------
% Author: Jing Wang (w1012j@163.com)
% Last Date: 17/07/2021
% -------------------------------------------------------
%%
function  Mp = f_expand(M,b1,b2)
% reflective expanding a matrix
% M: mxn
% b1, b2 must be ODD number;
[m,n,~] = size(M);
d1 = floor(b1/2);  d2 = floor(b2/2);
Up = M(d1:-1:1,:,:,:,:);
Down = M(m:-1:m-d1+1,:,:,:,:);
Mp= [Up;M;Down];
Left = Mp(:,d2:-1:1,:,:,:);
Right = Mp(:,n:-1:n-d2+1,:,:,:);
Mp = [Left,Mp,Right];
