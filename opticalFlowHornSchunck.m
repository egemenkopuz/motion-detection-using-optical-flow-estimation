function [Vx,Vy] = opticalFlowHornSchunck(frame0, frame1, alpha, iter)
%% SETUP
[ydim,xdim,ch] = size(frame0);
if (ch == 3)
    frame0 = rgb2gray(frame0);
end
[ydim_next,xdim_next,ch] = size(frame1);
if (ch == 3)
    frame1 = rgb2gray(frame1);
end
%% ERROR CHECK
if (ydim ~= ydim_next || xdim ~= xdim_next)
    throw(MException('Inconsistent frames','frame dimensions must be same.'));
end
if (iter <= 0)
    throw(MException('Inconsistent parameter','iter > 0.'));
end
%% PROCESS

frame0 = double(frame0);
frame1 = double(frame1);

Vx = zeros(size(frame0(:,:,1)));
Vy = zeros(size(frame1(:,:,1)));

%Horn-Schunck Spatiotemporal derivative estimation
Ix = conv2(frame0, 0.25* [-1 1; -1 1],'same') + conv2(frame1, 0.25*[-1 1; -1 1],'same');
Iy = conv2(frame0, 0.25*[-1 -1; 1 1], 'same') + conv2(frame1, 0.25*[-1 -1; 1 1], 'same');
It = conv2(frame0, 0.25*ones(2),'same') + conv2(frame1, -0.25*ones(2),'same');

%Horn-Schunk Kernel Mask
kernel = [1/12 1/6 1/12; 1/6 0 1/6; 1/12 1/6 1/12];

for i=1:iter
    AvgU = conv2(Vx,kernel,'same');
    AvgV = conv2(Vy,kernel,'same');
    
    Vx = AvgU - (Ix.*((Ix.*AvgU)+(Iy.*AvgV) + It))./(alpha^2+Ix.^2+Iy.^2);
    Vy = AvgV - (Iy.*((Ix.*AvgU)+(Iy.*AvgV) + It))./(alpha^2+Ix.^2+Iy.^2);
end