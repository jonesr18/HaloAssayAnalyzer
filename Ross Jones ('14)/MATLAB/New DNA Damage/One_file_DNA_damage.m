%% Image processing
close all
clear all
clc

% Ross Jones
% Single file DNA damage analysis
% Updated 1.11.14

% Import image
[file, filePath] = uigetfile('*.jpeg; *.png; *.tif', 'Pick image file(s)', 'MultiSelect', 'off');
% filePath = 'C:\Users\Ross\Google Drive\Capstone Images\Binh\';
% file = 'Array Image-1.tif';
image = rgb2gray(imread(strcat(filePath, file)));
image = imfilter(image, fspecial('average'));
figure(), imshow(image), title('Original Image')

% Threshold with multithresh, then combine two thresholds to make better binary image and make a
% binary image that only uses the nucleus threshold as well for later on. 
thresh = multithresh(image, 2);
imageThresh = uint8(imquantize(image, thresh));
binary = boolean(imageThresh > 1);
figure(), imshow(label2rgb(imageThresh)), title('Thresholded Image')
% figure(), imshow(binary), title('Adjusted Binary Image')

% Smoothen the image objects
strucElemDiam = strel('diamond', 1);
binary = imerode(binary, strucElemDiam);
binary = imerode(binary, strucElemDiam);
% figure(), imshow(binary), title('Smoothed Image')

% Find the boundaries
[boundaries, labels] = bwboundaries(binary, 'noholes');
figure(), imshow(binary), title('All Boundaries'), hold on
for i = 1:length(boundaries)
    boundary = boundaries{i};
    plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2)
end

% Determine which objects are round (just one cell)
stats = regionprops(labels, 'Area', 'Centroid', 'Eccentricity');
threshold = 0.510;
for i = 1:length(boundaries)
    boundary = boundaries{i};

    % compute the eccentricity metric
    eccentricity = stats(i).Eccentricity;

    % display the results
    eccenString = sprintf('%.3f', eccentricity);

    % mark objects below the threshold with a black circle
    if eccentricity < threshold && stats(i).Area > 50
        centroid = stats(i).Centroid;
        plot(centroid(1), centroid(2), 'ko');
    else
        binary = binary - boolean(labels == i);
    end

    text(boundary(1, 2) - 35, boundary(1, 1) + 50, eccenString, 'Color', 'g', 'FontWeight', 'bold');
end

% figure(), imshow(binary), title('Regions Removed')

% Remove regions from thresholded image --- note that removing from greyscale and rethresholding can
% cause issues, as sometimes there become dicontinginous regions. 
imageThresh = imageThresh .* uint8(binary);

% Obtain binary images containing all halos or all nuclei (note: halos include nuclei)
halos  = boolean(imageThresh > 1);
nuclei = boolean(imageThresh > 2);
% figure(), imshow(halos), title('Halos')
% figure(), imshow(nuclei), title('Nuclei')


%% Analysis
% close all
% clc

% Find areas of each nucleus and halo
haloStats = regionprops(halos, 'Area', 'Centroid', 'PixelList');
nucStats  = regionprops(nuclei, 'Area', 'Centroid', 'EquivDiameter', 'PixelList');
numHalos  = length(haloStats);
numNuclei = length(nucStats);
pair  = zeros(numNuclei, 1);
taken = zeros(numHalos, 1);
distance = pair;

% Find the closest centroid in the opposing list in order to ensure match
% # of nuclei is used as baseline because # halo >= # nuclei, so we only want halos w/ nuclei
%   - This should be changed to find which one is smaller and go with that
for i = 1:numNuclei
    min = 9999;
    for j = 1:numHalos
        if ~taken(j)
            dist = sum(abs(nucStats(i).Centroid - haloStats(j).Centroid));
            if dist < min
                min = dist;
                match = j;
            end
        end
    end
    pair(i) = match;        % Nuclear centroids: index --- Halo centroids: value
    taken(match) = true;
    distance(i) = min;
end

% If centroids are not within 30 pixels in x or y, print coordinates
nearby = boolean(distance < 30);
for i = 1:numNuclei
    if ~nearby(i)
        fprintf('Xn: %.2f, Yn: %.2f\n', nucStats(i).Centroid)
        fprintf('Xh: %.2f, Yh: %.2f\n', haloStats(pair(i)).Centroid)
    end
    fprintf('\n')
end

% Plot centroids with number next to it to compare (yellow and blue numbers should match)
figure(), imshow(image), title('Overlaid Centroids')
hold on
for i = 1:numNuclei
    cn = nucStats(i).Centroid;
    ch = haloStats(i).Centroid;
    cp = haloStats(pair(i)).Centroid;
    text(cn(1), cn(2), num2str(i), 'Color', 'r')
    text(ch(1), ch(2), num2str(i), 'Color', 'b')
    text(cp(1), cp(2), num2str(i), 'Color', 'y')
end
hold off

% Calculate damage per halo / nucleus pair
damage = zeros(numNuclei, 1);
for i = 1:numNuclei
    damage(i) = haloStats(pair(i)).Area / nucStats(i).Area;
end
meanNDF = mean(damage);                     % Nuclear Diffusion Factor as measured by Ma et al.
sterrNDF = std(damage) / sqrt(numNuclei);   % Standard error of the mean (SEM)

% Collect pixel coordinates that are only in halo region
haloPixels = struct;
for i = 1:numNuclei
    allPixels = haloStats(pair(i)).PixelList;
    nucPixels = nucStats(i).PixelList;
    remainingPixels = zeros(length(allPixels) - length(nucPixels), 2);
    k = 1;
    s = 1;
    for j = 1:length(allPixels)
        if allPixels(j, :) == nucPixels(k, :);
            if k < length(nucPixels)
                k = k + 1;
            end
            s = s + 1;
        else
            remainingPixels(j - s + 1, :) = allPixels(j, :);
        end               
    end
    haloPixels(i).PixelList = remainingPixels;
end

% Calculate halo moment (intensity * radius)
%   This essentially tells gives us a weighted metric of the DNA spread outside the nucleus.
% Calculate adjusted moment (intensity * radius / (# of pixels * nuclear radius))
%   This adjusts the metric to account for different resolution images - comparing different
%   resolution images should be avoided in general though.
% Calculate halo distance (radius / (# of pixels * nuclear radius))
%   This tells us the distance that a given pixel is from the nucleus, also adjusted for resolution.
haloMoments    = zeros(numNuclei, 1);
adjHaloMoments = haloMoments;
haloDistances  = haloMoments;
for i = 1:numNuclei
    cn         = nucStats(i).Centroid;
    nucRadius  = nucStats(i).EquivDiameter / 2;
    coords     = haloPixels(i).PixelList;
    numPixels  = length(coords);
    moments    = zeros(numPixels, 1);
    adjMoments = moments;
    distances  = moments;
    for j = 1:numPixels
        x = cn(1) - coords(j, 1);   % x distance
        y = cn(2) - coords(j, 2);   % y distance
        distance      = sqrt(x.^2 + y.^2); 
        intensity     = double(image(coords(j, 2), coords(j, 1)));
        moments(j)    = intensity * distance;
        adjMoments(j) = intensity * distance / nucRadius;
        distances(j)  = distance / nucRadius;
    end
    haloMoments(i)    = sum(moments);
    adjHaloMoments(i) = mean(adjMoments);
    haloDistances(i)  = mean(distances);
end
meanHaloMoment     = mean(haloMoments);
meanAdjHaloMoment  = mean(adjHaloMoments);
meanHaloDistance   = mean(haloDistances);
sterrHaloMoment    = std(haloMoments) / sqrt(numNuclei);
sterrAdjHaloMoment = std(adjHaloMoments) / sqrt(numNuclei);
sterrHaloDistance  = std(haloDistances) / sqrt(numNuclei);

% Display damage and halo moment next to each nucleus
figure(), imshow(image), title('NDF and Halo Moment for Each Measured Cell')
hold on
for i = 1:numNuclei
    cn = nucStats(i).Centroid;
    text(cn(1), cn(2), num2str(damage(i)), 'Color', 'g')
    text(cn(1), cn(2) + 25, sprintf('%.5g', haloMoments(i)), 'Color', 'y')
end
hold off

% Compute statistics for entire iamge simultaneous for comparison
threshStats   = regionprops(imageThresh, 'Area');      % Find areas of thresholded regions
nucSurface    = threshStats(3).Area;
haloSurface   = threshStats(2).Area;
overallDamage = (haloSurface + nucSurface) / nucSurface;


fprintf('NDF:\nTotal: %.4f\nMean: %.4f\nSEM: %.4f\n\n', overallDamage, meanNDF, sterrNDF)
fprintf('Halo Moment:\nMean: %.4g\nSEM: %.4g\n\n',meanHaloMoment, sterrHaloMoment)
fprintf('Adjusted Halo Moment:\nMean: %.4g\nSEM: %.4g\n\n', meanAdjHaloMoment, sterrAdjHaloMoment);
fprintf('Halo Distance:\nMean: %.4g\nSEM: %.4g\n\n', meanHaloDistance, sterrHaloDistance);


% -------------------------------------------------------------------
%
% oldDamage      = 3.6448; % Computed without removing large cell clusters
% oldDamage2     = 4.0641; % Computed removing cell cluseters using roundness metric
% oldMeanNDF     = 5.1986; % Same as above
% oldSterrNDF    = 0.2764; % " "
% oldMeanMoment  = 3.848 e+7; % " "
% oldSterrMoment = 8.67 e+6; % " "





