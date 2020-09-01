%% Input (Settings and video)
METRIC_MODE = false;         % to measure performance (disables technical plots)
PLOT_ONLY_TRACKING = false;        % plots only visual tracking

filename = 'video_example.mp4';  
frameSkipRate = 3;          % 2 for metric mode, otherwise 6 or 12
frameCompareInterval = 1;   % 1 or 2

% Implementation Method, sparse(features) or dense(spatial gradients)
implementationMethod = "dense" ;    % "dense" or "sparse"

% for dense method
opticalMethod = "Horn-Schunck";    % "Lucas-Kanade" or "Horn-Schunck"

% for sparse method
magnitudeScalar = 0.01;

% for Lucas-Kanade:
threshold = 50000;
windowSize = 15;

% for Horn-Schunck:
alpha = 50;
iter = 10;

% Blob Analysis
minBlobArea = 200;
maxBlobArea = 15000;

%% Setup
if(~isdeployed)
  cd(fileparts(which(mfilename)));
end
try
video = VideoReader(filename);
catch exception
    throw(exception);
end
nframes = video.NumFrames;
frame0 = read(video,1);
[ydim,xdim,~] = size(frame0);
[X,Y] = meshgrid(1:1:xdim,1:1:ydim);

if (implementationMethod == "sparse")
   pointTracker = vision.PointTracker('MaxBidirectionalError', 2); 
   pointTrackerInit = true;
end

% Morphological operation parameters
seDense = strel('disk',10); % for dense methods
seSparse = strel('disk',20); % for sparse methods
seHeat = strel('disk',25); % for heat map

% Blob Analysis
BlobAnalysis = vision.BlobAnalysis('MinimumBlobArea',minBlobArea,'MaximumBlobArea',maxBlobArea);

%% Main Loop
tic; % for performance measures
iterCount = 0;

for i=1+frameCompareInterval:frameSkipRate:nframes
    
    if (i > 220)
       breakpoint = 5; 
    end
    % process
    frame0 = read(video,i-frameCompareInterval);
    frame1 = read(video,i);
    
    % colored frame
    frame_rgb = frame1;
    
    % from colored to grayscale along with gaussian filtering
    frame0 = imgaussfilt(rgb2gray(frame0),1);
    frame1 = imgaussfilt(rgb2gray(frame1),1);
    
    if (implementationMethod == "sparse")
       % Sparse method
       [Vx,Vy,Points] = opticalFlowFeatures(pointTracker,frame0,frame1,magnitudeScalar,pointTrackerInit);
        pointTrackerInit = false;
    else
       % Dense method
       if (opticalMethod == "Lucas-Kanade")
           % parameters:(previous frame, current frame, treshold, window size)
           [Vx,Vy] = opticalFlowLucasKanade(frame0,frame1,threshold,windowSize);
       else
           % parameters:(previous frame, current frame, alpha, iteration count)
           [Vx,Vy] = opticalFlowHornSchunck(frame0,frame1,alpha,iter);
       end
    end
    
    % magnitudes
    Vm = sqrt(Vx.^2 + Vy.^2);
    % to binary image
    binarized = imbinarize(Vm);
    % morphological closing (dilation)
    if (implementationMethod == "sparse")
        BW = imclose(binarized,seSparse);
    else
        BW = imclose(binarized,seDense);
    end
    
    % blob analysis
    [area,centroid,bbox] = step(BlobAnalysis,BW);
    
    % visiualization
    if METRIC_MODE == false
        if PLOT_ONLY_TRACKING == true
        figure(1);
        shape = insertShape(frame_rgb,'rectangle',bbox,'linewidth',2,'color','y');
        imagesc(shape);
        title('Visual Tracking')
        else
        figure(1);

        subplot(2,2,1)
        imagesc(frame_rgb); 
        hold on;
        quiver(X,Y,Vx,Vy,23,'y');
        hold off;
        title('Optical Flow')
        

        axH = subplot(2,2,2);
        imshow(BW)
        title('Morphed')
        colormap (axH,gray);

        axM = subplot(2,2,3);
        imagesc(imclose(Vm,seHeat),[0,1.0]);
        title('Magnitude')
        colorbar;
        colormap(axM,jet);

        subplot(2,2,4)
        shape = insertShape(frame_rgb,'rectangle',bbox,'linewidth',2,'color','y');
        imagesc(shape);
        title('Visual Tracking')
        end
    end
    
    % TO PLOT CORNERS
%   if (implementationMethod == "sparse")
%       figure(2);
% 
%       subplot(1,1,1);
%       imshow(frame_rgb);
%       hold on;
%       plot(Points(:,1),Points(:,2),'r*');
%       hold off;
%   end
    % TO PLOT UNMORPHED
%   if (implementationMethod == "dense")
%       figure(2);
% 
%       subplot(1,1,1);
%       imshow(BW);
%   end 
    
    iterCount = iterCount + 1;
end

toc;
X = ['Average elapsed time per frame is ', num2str(toc/iterCount), ' seconds'];
disp(X);
release(BlobAnalysis);