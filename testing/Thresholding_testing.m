%% Thresholding Comparisons
close all
clear all
clc

warning('off', 'images:initSize:adjustingMag');
% image = imread('C:\Users\Ross\Google Drive\Capstone Images\Binh\Array Image-1.tif');
% image = imread('C:\Users\Ross\Google Drive\Capstone Images\Binh\X-ray for A172.tif');
% image = imread('C:\Users\Ross\Google Drive\Capstone Images\Binh\nanoparticle for HeLa.png');
% image = imread('C:\Users\Ross\Google Drive\Capstone Images\Binh\metal ion for HeLa.png');
image = imread('C:\Users\Ross\Google Drive\Capstone Images\Binh\figure pic.png');
% image = imread('C:\Users\Ross\Google Drive\Capstone Images\87 - 12.26.13 - FHA T5 (EB)\Images\Image-51.tif');
% image = imread('C:\Users\Ross\Google Drive\Capstone Images\97 - 05.11.14 0 D1 (EB)\Images\Image-35.tif');
% image = imread('C:\Users\Ross\Google Drive\Capstone Images\97 - 05.11.14 200 D1 (EB)\Images\Image-7.tif');

image = imfilter(rgb2gray(image), fspecial('average'));
% image = image(1:454, 1:708); % Crop our microscope images


%% Multithresh
close all
clc

tic
% Threshold into 3 levels w/ multithresh (Otsu's method)
otsu3 = imquantize(image, multithresh(image, 2));
totsu3 = toc;

tic
% Threshold into 4 levels w/ multithresh
otsu4 = imquantize(image, multithresh(image, 3));
totsu4 = toc

extract = uint8(imclearborder(boolean(otsu4 > 1))) .* image;
figure(), imshow(label2rgb(otsu4))
figure(), imshow(extract, [])


%% K-means
close all
clc

tic
% Convert to a*b* space
cform = makecform('srgb2xyz');
Lab = applycform(image, cform);
ab = double(Lab(:, :, 1:3));
nrows = size(ab, 1);
ncols = size(ab, 2);
reimage = reshape(ab, nrows*ncols, 3);

% Threshold into 3 levels w/ k-means
IDX = kmeans(reimage, 3, 'distance', 'sqEuclidean', 'Replicates', 3);
kmeans3 = reshape(IDX, nrows, ncols);
ktime = toc;

% Threshold into 4 levels w/ k-means
IDX = kmeans(reimage, 4, 'distance', 'sqEuclidean', 'Replicates', 3);
kmeans4 = reshape(IDX, nrows, ncols);

% Plot results
titles = {'Otsu''s Method -- 3 Levels', 'Otsu''s Method -- 4 Levels', 'K-means -- 3 levels', 'K-means -- 4 levels'};
images = {otsu3, otsu4, kmeans3, kmeans4};
for i = 1:4
    figure()
    imshow(image, []), title(titles{i}), hold on
    himage = imshow(label2rgb(images{i}));
    set(himage, 'AlphaData', 0.3);
    hold off
end

% Display time results
fprintf('Otsu time = %.4f seconds\nK-means time = %.4f seconds\n', otime, ktime)


%% Gradient magnitude
close all
clc

tic
% Get gradient magnitude
hy = fspecial('sobel');
Iy = imfilter(double(image), hy, 'replicate');
Ix = imfilter(double(image), hy', 'replicate');
gradient = sqrt(Ix.^2 + Iy.^2);
tgradmag = toc

figure(), imshow(gradient, [min(gradient(:)), max(gradient(:)) / 4])

tic
% Threshold gradient magnitude
imthresh = imquantize(gradient, multithresh(gradient, 2));
gradient = boolean(imthresh > 2);

% Dilate the image -- closing and opening smoothes edges
strel1   = strel('disk', 1);
gradient = imopen(gradient, strel1);
gradient = imdilate(gradient, strel('disk', 3));
gradient = imdilate(gradient, strel('diamond', 3));
gradient = imopen(gradient, strel1);
gradient = imclose(gradient, strel1);
gradient = imopen(gradient, strel1);
gradient = imclose(gradient, strel1);
tnucgrad = toc + tgradmag

% --- %

tic
gradient2 = boolean(imthresh > 1);

% Dilate the image -- closing and opening smoothes edges
strel1   = strel('disk', 1);
gradient2 = imopen(gradient2, strel1);
gradient2 = imdilate(gradient2, strel('disk', 3));
gradient2 = imdilate(gradient2, strel('diamond', 3));
gradient2 = imopen(gradient2, strel1);
gradient2 = imclose(gradient2, strel1);
gradient2 = imopen(gradient2, strel1);
gradient2 = imclose(gradient2, strel1);
tallgrad = toc + tgradmag

figure(), imshow(imclearborder(gradient))
figure(), imshow(image, [])
figure(), imshow(imclearborder(gradient2))
extract = uint8(imfill(imclearborder(gradient2), 'holes')) .* image;
figure(), imshow(extract, [])

