function [Vx,Vy] = opticalFlowLucasKanade(frame0, frame1, t, w)
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
if (w <= 1)
    throw(MException('Inconsistent parameter','w must be > 1.'));
end
if (t <= 0)
    throw(MException('Inconsistent parameter','t must be > 0.'));
end
%% PROCESS

frame0 = double(frame0);
frame1 = double(frame1);

% spatial gradients 
[Ix, Iy] = imgradientxy(frame1,'Prewitt');
% temporal gradient
It = frame1 - frame0;
% velocity vectors
Vx = zeros(ydim,xdim);
Vy = zeros(ydim,xdim);

w = round(w/2);
for i=w+1:w:xdim-w-1
   for j=w+1:w:ydim-w-1 
       valx = Ix(j-w:j+w,i-w:i+w);
       valy = Iy(j-w:j+w,i-w:i+w);
       valt = It(j-w:j+w,i-w:i+w);
       
       g1 = sum(sum(valx.*valx));
       g2 = sum(sum(valx.*valy));
       g3 = sum(sum(valy.*valy));
       G = [g1,g2;g2,g3];
       
       b1 = sum(sum(valx.*valt));
       b2 = sum(sum(valy.*valt));
       b =[b1;b2];
       
       eigs = eig(G);
       if (min(eigs) >= t)
           u = -inv(G)*b;
           Vx(j,i) = u(1);
           Vy(j,i) = u(2);
       else
           Vx(j,i) = 0;
           Vy(j,i) = 0;
       end
   end
end