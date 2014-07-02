% Ross Jones
% 10.8.13
% Images Testing


%% Image Import
clc

countsTotal = 0;
xTotal = 0;

for i = 1:2:59
    ii = num2str(i);
    fileFolder = ['C:\Users\Ross\Google Drive\Capstone Images\87 - 10.7.13 - FHA 1 (EB)\Images\Image-' ii '.tif'];
    [image, colormap] = imread(fileFolder);
    image = rgb2gray(image);
    [counts, x] = imhist(image);
    countsTotal = countsTotal + counts;
    xTotal = xTotal + x;
end

figure(), stem(xTotal, countsTotal)


%% Process a sequence of images
clc

% This example starts by reading a series of images from a directory into the 
% MATLAB® workspace, storing the images in an m * n * p array. The example 
% then passes the entire array to the stdfilt function and performs standard 
% deviation filtering on each image in the sequence. Note that, to use stdfilt 
% with an image sequence, you must use the nhood argument, specifying a 2-D neighborhood.

% Create an array of filenames that make up the image sequence
fileFolder = 'C:\Users\Ross\Google Drive\Capstone Images\87 - 10.7.13 - FHA 1 (EB)\Images';
% fileFolder = fullfile(matlabroot, 'toolbox', 'images', 'imdemos');
dirOutput = dir(fullfile(fileFolder, 'Image-*.tif'));
fileNames = {dirOutput.name}';
numFrames = numel(fileNames);

image = rgb2gray(imread(fullfile(fileFolder, fileNames{1})));

% Preallocate the array
sequence = zeros([size(image) numFrames], class(image));

% Create image sequence array
for p = 1:numFrames
    sequence(:, :, p) = rgb2gray(imread(fullfile(fileFolder, fileNames{p}))); 
end

% Process sequence
sequenceNew = stdfilt(sequence, ones(3));

% View results
figure()
for k = 1:numFrames
      imshow(sequence(:, :, k));
      title(sprintf('Original Image # %d',k));
      pause(1);
      imshow(sequenceNew(:, :, k), []);         % [] is for high-low, blank means uses image high/low
      title(sprintf('Processed Image # %d',k));
      pause(1);
end


%% Multithresh
close all
clear all
clc

% returns thresh containing N threshold values using Otsu's method. thresh is a 
% 1xN vector which can be used to convert image A into an image with N + 1 
% discrete levels using imquantize.

% thresh = multithresh(A,N);

% returns metric, a measure of the effectiveness of the computed thresholds. 
% metric is in the range [0 1] and a higher value indicates greater 
% effectiveness of the thresholds in separating the input image into N + 1 
% classes based on Otsu's objective criterion.

% [thresh,metric] = multithresh(___) 

% Example

% Import image
filePath = 'C:\Users\Ross\Google Drive\Capstone Images\Binh\';
file = 'Array Image-1.tif';
image = rgb2gray(imread(strcat(filePath, file)));
figure(), imshow(image), title('Original Image')

% Implement multithresh
[thresh, metric] = multithresh(image, 3);
fprintf('Metric = %.3f', metric)
segImage = imquantize(image, thresh);               % apply the thresholds to obtain segmented image
RGB = label2rgb(segImage);                          % convert to color image
figure(), imshow(RGB), title('Thresholded Image')


%% Cell Objoect Identification Example
close all
clear all
clc

% Import image
image = imread('C:\Users\Ross\Google Drive\Capstone Images\87 - 10.7.13 - DDA 1 (YOYO)\Images\Image-31.tif');
image = rgb2gray(image);
figure(), imshow(image), title('Original Image')

% Detect entire cell
[junk, threshold] = edge(image, 'sobel');
fudgeFactor = .5;
BWs = edge(image, 'sobel', threshold * fudgeFactor);
figure(), imshow(BWs), title('Binary Gradient Mask')

% Dilate the image
strucElem90 = strel('line', 3, 90);
strucElem0  = strel('line', 2, 0);
BWs_dil = imdilate(BWs, [strucElem90, strucElem0]);
figure(), imshow(BWs_dil), title('Dilated Gradient Mask')

% Fil interior gaps
BWd_fill = imfill(BWs_dil, 'holes');
figure(), imshow(BWd_fill), title('Binary Image with Filled Holes')

% Remove connected objects on border
% BW_no_bord = imclearborder(BWd_fill, 4);
% figure(), imshow(BW_no_bord), title('Cleared Border Image')

% Smoothen the object
strucElemDiam = strel('diamond', 1);
BW_final = imerode(BWd_fill, strucElemDiam);
BW_final = imerode(BW_final, strucElemDiam);
figure(), imshow(BW_final), title('Segmented Image')

% Outline the object
BW_outline = bwperim(BW_final);
segout = image;
segout(BW_outline) = 255;
figure(), imshow(segout), title('Outlined Original Image')


%% Harmonic Cut & Regularized Centroid Transform
% Yang & Parvin (2003)
close all;
clear all;
clc

% Import image
image = imread('C:\Users\Ross\Google Drive\Capstone Images\87 - 11.21.13 - DDA C1 (EB)\Images\Image-8.tif');
image = rgb2gray(image);
figure(), imshow(image), title('Original Image')

[H, W] = size(image);

tic
% Create offset matrices
C  = zeros(H, W); C(2:H-1, 2:W-1)  = image(2:H-1, 2:W-1);
Cl = zeros(H, W); Cl(2:H-1, 2:W-1) = image(2:H-1, 1:W-2);
Cr = zeros(H, W); Cr(2:H-1, 2:W-1) = image(2:H-1, 3:W);
Cu = zeros(H, W); Cu(2:H-1, 2:W-1) = image(1:H-2, 2:W-1);
Cd = zeros(H, W); Cd(2:H-1, 2:W-1) = image(3:H, 2:W-1);

% Create edge vectors
L1 = image(2:H-1, 1);
L2 = image(2:H-1, 2);
L3 = image(2:H-1, 3);
L4 = image(2:H-1, 4);
R1 = image(2:H-1, W);
R2 = image(2:H-1, W-1);
R3 = image(2:H-1, W-2);
R4 = image(2:H-1, W-3);
U1 = image(1, 2:W-1);
U2 = image(2, 2:W-1);
U3 = image(3, 2:W-1);
U4 = image(4, 2:W-1);
D1 = image(H, 2:W-1);
D2 = image(H-1, 2:W-1);
D3 = image(H-2, 2:W-1);
D4 = image(H-3, 2:W-1);

% Ix, Iy, Ixx, & Iyy internals
Ix = (Cr - Cl) / 2;
Iy = (Cd - Cu) / 2;
Ixx = Cr - 2 * C + Cl;
Iyy = Cd - 2 * C + Cu;

% Ix, Iy, Ixx, & Iyy edges
Ix(2:H-1, 1)  = (-3 * L1 + 4 * L2 - L3) / 2;
Ix(2:H-1, W)  = (-3 * R1 + 4 * R2 - R3) / 2;
Iy(1, 2:W-1)  = (-3 * U1 + 4 * U2 - U3) / 2;
Iy(H, 2:W-1)  = (-3 * D1 + 4 * D2 - D3) / 2;
Ixx(2:H-1, 1) = 2 * L1 - 5 * L2 + 4 * L3 - L4;
Ixx(2:H-1, W) = 2 * R1 - 5 * R2 + 4 * R3 - R4;
Iyy(1, 2:W-1) = 2 * U1 - 5 * U2 + 4 * U3 - U4;
Iyy(H, 2:W-1) = 2 * D1 - 5 * D2 + 4 * D3 - D4;

% new Cu, Cd, Us, and Ds for Ixy
Cu(2:H-1, 2:W-1) = Ix(1:H-2, 2:W-1); 
Cd(2:H-1, 2:W-1) = Ix(3:H, 2:W-1);
U1 = Ix(1, 2:W-1);
U2 = Ix(2, 2:W-1);
U3 = Ix(3, 2:W-1);
D1 = Ix(H, 2:W-1);
D2 = Ix(H-1, 2:W-1);
D3 = Ix(H-2, 2:W-1);

% Ixy internals and edges
Ixy = (Cd - Cu) / 2;
Ixy(1, 2:W-1) = (-3 * U1 + 4 * U2 - U3) / 2;
Ixy(H, 2:W-1) = (-3 * D1 + 4 * D2 - D3) / 2;

diffx = [
    0 0 0;
   -1 0 1;
    0 0 0;
    ];

diffy = [
    0 -1 0;
    0  0 0;
    0  1 0;
    ];

ix = conv2(image, diffx, 'same');
iy = conv2(image, diffy, 'same');
ixx = conv2(ix, diffx, 'same');
ixy = conv2(ix, diffy, 'same');
iyy = conv2(iy, diffy, 'same');

% Elliptic feature
Hxy1 = boolean(Ixx .* Iyy - Ixy.^2 > 0);
Hxy2 = boolean(Ixx + Iyy > 0);

hxy1 = boolean(ixx .* iyy - ixy.^2 > 0);
hxy2 = boolean(ixx + iyy > 0);

ellip1 = Hxy1 .* double(image);
Ellip1 = hxy1 .* double(image);
ellip2 = Hxy2 .* double(image);
Ellip2 = hxy2 .* double(image);
% ellip = Hxy1 .* Hxy2 .* double(image);
toc

% figure(), imshow(Hxy1), title('Hxy1')
figure(), imshow(Hxy2), title('Where Ixx + Iyy > 0')
figure(), imshow(hxy2), title('Where ixx + iyy > 0')
% figure(), imshow(ellip, [0 255]), title('Both Filtered')
figure(), imshow(ellip1, [0 255]), title('Hxy1 Filtered')
figure(), imshow(Ellip1, [0 255]), title('hxy1 Filtered')
figure(), imshow(ellip2, [0 255]), title('Hxy2 Filtered')
figure(), imshow(Ellip2, [0 255]), title('hxy2 Filtered')

% RCT

tic
% Precalculate constants
A   = 1;                          % alpha (weight = 1)
D   = 1;                          % delta (change = 1 pixel)
a   = Ixx.^2 + Ixy.^2 + 4 * A;
b   = Ixx .* Ixy + Ixy .* Iyy;
c   = b;
d   = Ixy.^2 + Iyy.^2 + 4 * A;
e0  = -Ix .* Ixx - Iy .* Ixy;
f0  = -Ix .* Ixy - Iy .* Iyy;
u0  = (d .* e0 - b .* f0) / D;
v0  = (-c .* e0 + a .* f0) / D;
k11 = 4 * A * d / D;
k12 = -4 * A * b / D;
k21 = -4 * A * c / D;
k22 = 4 * A * a / D;
toc

tic
% initialize RCT vectors
iter = 75;
u    = zeros(H, W, iter);
v    = u;
w    = zeros(H, W);
uL   = zeros(H, W);
uR   = uL;
uU   = uL;
uD   = uL;
vL   = uL;
vR   = uL;
vU   = uL;
vD   = uL;
toc
%%
tic
for i = 1:iter-1
    % New shifted matrices
    uL(2:H-1, 2:W-1) = u(2:H-1, 1:W-2, i);
    uR(2:H-1, 2:W-1) = u(2:H-1, 3:W, i);
    uU(2:H-1, 2:W-1) = u(1:H-2, 2:W-1, i);
    uD(2:H-1, 2:W-1) = u(3:H, 2:W-1, i);
    vL(2:H-1, 2:W-1) = v(2:H-1, 1:W-2, i);
    vR(2:H-1, 2:W-1) = v(2:H-1, 3:W, i);
    vU(2:H-1, 2:W-1) = v(1:H-2, 2:W-1, i);
    vD(2:H-1, 2:W-1) = v(3:H, 2:W-1, i);
    
    % Average around point
    u_bar = (uL + uR + uU + uD) / 4;
    v_bar = (vL + vR + vU + vD) / 4;
    
    % Interpolate next point
    u(:, :, i+1) = u0 + k11 .* u_bar + k12 .* v_bar;
    v(:, :, i+1) = v0 + k21 .* u_bar + k22 .* v_bar;
    
    % Test for equilibrium
    w = u(:, :, i+1) + v(:, :, i+1);
    if i > 25
        for j = 1:H
            for k = 1:W
                if w(j, k) == 0
                    fprintf('Eq point @ (%u, %u, %u)\n', j, k, i)
                end
            end
        end
    end
end
toc

figure(), imshow(u(:, :, end))
figure(), imshow(v(:, :, end))


%% Detect and measure circular objects in an image
close all
clear all
clc

% From MATLAB examples online:
% http://www.mathworks.com/help/images/examples/detect-and-measure-circular-objects-in-an-image.html?prodcode=IP&language=en

% This seems to have identified virtually all of the healthy cells in the image, including
% overalapping ones. It did not identify apoptotic cells. It may be overdoing it on identifying
% two cells when really it is just an elongated cell

% Import image
image = imread('C:\Users\Ross\Google Drive\Capstone Images\87 - 11.21.13 - DDA C1 (EB)\Images\Image-8.tif');
image = rgb2gray(image);
figure(), imshow(image), title('Original Image')

tic
% Find circles
radUp  = 20;
radLow = 10;
[centers, radii] = imfindcircles(image, [radLow radUp], 'ObjectPolarity', 'bright', 'Sensitivity', 0.92);
h = viscircles(centers, radii);
toc

tic
delta = 1;
while delta > .001
    radAvg = mean(radii);
    radStd = std(radii);
    radUp  = round(radAvg + pi * radStd);
    radLow = round(radAvg - pi * radStd);
    [centers, radNew] = imfindcircles(image, [radLow radUp], 'ObjectPolarity', 'bright', 'Sensitivity', 0.92);
    radNewAvg = mean(radNew);
    delta = abs(radAvg - radNewAvg);
    radii = radNew;
end
toc

fprintf('Mean: %f\nStdev: %f\nLower bound: %i\nUpper bound: %i\n', mean(radii), std(radii), radLow, radUp)
figure(), imshow(image), title('Optimized Circles')
h2 = viscircles(centers, radii);


%% Find round objects
close all
clear all
clc

% From MATLAB examples online: 
% http://www.mathworks.com/help/images/examples/identifying-round-objects.html?prodcode=IP&language=en

% After evaluation, this method will not work, as it meshes together nearby and touching cells,
% causing them to appear to be one shape, and thus non-circular.

% Import image
image = imread('C:\Users\Ross\Google Drive\Capstone Images\87 - 11.21.13 - DDA C1 (EB)\Images\Image-9.tif');
image = rgb2gray(image);
figure(), imshow(image), title('Original Image')

% Create binary image
threshold = graythresh(image);
binary = im2bw(image, threshold);
figure(), imshow(binary), title('Binary Image')

% Remove the noise
area = 30;                              % Size of areas to be removed
binary = bwareaopen(binary, area);
strucElem = strel('disk', 2);           % Make structuring element to fill gaps
% binary = imclose(binary, strucElem);
binary = imfill(binary, 'holes');       % Fill in holes
figure(), imshow(binary), title('Noise Removed')

% Find the boundaries
[boundaries, labels] = bwboundaries(binary, 'noholes');
figure(), imshow(binary), title('All Boundaries'), hold on
for i = 1:length(boundaries)
    boundary = boundaries{i};
    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2)
end

% Determine which objects are round
stats = regionprops(labels, 'Area', 'Centroid');
threshold = 0.80;
for i = 1:length(boundaries)
    % obtain (X,Y) boundary coordinates corresponding to label 'k'
    boundary = boundaries{i};

    % compute a simple estimate of the object's perimeter
    delta_sq = diff(boundary).^2;
    perimeter = sum(sqrt(sum(delta_sq, 2)));

    % obtain the area calculation corresponding to label 'k'
    area = stats(i).Area;

    % compute the roundness metric
    metric = 4 * pi * area / perimeter^2;

    % display the results
    metricString = sprintf('%2.2f', metric);

    % mark objects above the threshold with a black circle
    if metric > threshold
        centroid = stats(i).Centroid;
        plot(centroid(1), centroid(2), 'ko');
    end

    text(boundary(1,2)-35, boundary(1,1)+13, metricString, 'Color', 'y', 'FontWeight', 'bold');

end




%% K-means Clustering
close all
clear all
clc

% From MATLAB examples online
% http://www.mathworks.com/help/images/examples/color-based-segmentation-using-k-means-clustering.html?prodcode=IP&language=en

% Works almost perfectly for FHA 2 (YOYO)
% Works terribly for DDA 1 (EB)
% Works well for DDA C1 (EB)

% Import image
image = imread('C:\Users\Ross\Google Drive\Capstone Images\87 - 10.7.13 - FHA 2 (YOYO)\Images\Image-1.tif');
% image = rgb2gray(image);
figure(), imshow(image), title('Original Image')

% Convert to a*b* space
cform = makecform('srgb2xyz');
Lab = applycform(image, cform);
ab = double(Lab(:, :, 1:3));
nrows = size(ab, 1);
ncols = size(ab, 2);
image = reshape(ab, nrows*ncols, 3);

% Perform clustering (3 times to avoid local minima)
[IDX, centroids, sumDists, dists] = kmeans(image, 3, 'distance', 'sqEuclidean', 'Replicates', 3);

% Label pixels with k-means results
labels = reshape(IDX, nrows, ncols);
figure(), imshow(labels, []), title('Image Labeled by Cluster Index')


%% Kernal testing
close all
clear all
clc

% Import image
image = imread('C:\Users\Ross\Google Drive\Capstone Images\87 - 10.7.13 - FHA 2 (YOYO)\Images\Image-1.tif');
image = rgb2gray(image);
figure(), imshow(image), title('Original Image')

laplacian = [
    0 1 0;
    1 -4 1;
    0 1 0;
    ];

laplacian2 = [
    1 1 1;
    1 -8 1;
    1 1 1;
    ];

sobelx = [
    1 0 -1;
    2 0 -2;
    1 0 -1;
    ];

sobely = [
    1 2 1;
    0 0 0;
    -1 -2 -1;
    ];

prewittx = [
    1 0 -1;
    1 0 -1;
    1 0 -1;
    ];

prewitty = [
    1 1 1;
    0 0 0;
    -1 -1 -1;
    ];

average = [
    1 1 1;
    1 1 1;
    1 1 1;
    ] / 9;


