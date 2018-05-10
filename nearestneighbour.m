clear
clc
A=dlmread(‘cc.txt');
[n,m] =size(A);
x=zeros(1,3);
d=0;
e=0;
Sum=0;
Ave=0;
pix=0.91; % pixel size in nm
fid = fopen(‘nnd.txt','w+');
%tic
for i=1: n
    i;
    x = A(i,:)';
    idx = nearestneighbour(x, A', 'NumberOfNeighbours', 2); % 2:nearest,3:second nearest, etc.
    d=norm(x-A(idx(1,2),:)’); % 2:nearest,3:second nearest, etc.
    e=d*pix;
    fprintf(fid, '%d\n',e);
    Sum= Sum+ d;
end
Ave=Sum/n*pix
fclose(fid);
%toc