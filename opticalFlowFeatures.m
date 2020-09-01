function [Vx,Vy,visiblePoints] = opticalFlowFeatures(pointTracker, frame0, frame1, magnitudeScalar, firstTime)
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
if (magnitudeScalar < 0)
    throw(MException('Inconsistent parameter','magnitudeScalar >= 0.'));
end
%% PROCESS

Vx = zeros(ydim,xdim);
Vy = zeros(ydim,xdim);

% detect points
%oldPoints = detectMinEigenFeatures(frame0);
oldPoints = detectHarrisFeatures(frame0);
% coordinates
points = oldPoints.Location;

if (firstTime == true)
    % first time init
    initialize(pointTracker, points, frame0);
else
    setPoints(pointTracker,points);
end

% check next locations of the points with Lucas-Kanade
[points, isFound] = step(pointTracker, frame1);
visiblePoints = points(isFound, :);
oldInliers = oldPoints(isFound, :);

coordinates = oldInliers.Location;
for i=1:1:size(visiblePoints,1)
    X = coordinates(i,1);
    Y = coordinates(i,2);
    
    newX = visiblePoints(i,1);
    newY = visiblePoints(i,2);
    
    valx = newX - X;
    valy = newY - Y;
    
    Vx(round(Y),round(X)) = valx * magnitudeScalar;
    Vy(round(Y),round(X)) = valy * magnitudeScalar;
end
end