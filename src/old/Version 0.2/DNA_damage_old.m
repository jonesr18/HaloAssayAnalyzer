%% Yong old version
close all
clear all
clc

warning('off', 'images:initSize:adjustingMag');

% Import image(s) and plot
imArray = importImage();
figure(), imshow(imArray{1}), title('Original Image')

% Perform image processing
imRegions = processImage(imArray);
pairs = findMatches(imRegions);

% Calculate results -- you can change it to only calculate the NDF by removing the other names,
% likewise you will also have to remove them from the -for- loop below. 
% This also shows a plot of DNA damage metrics overlaid on the image when -true- is a parameter.
results = calcDamage(imRegions, pairs, 'NDF', 'moment', 'adjMoment', 'distance', true);

% Print out statistics
for i = 1:length(results)
    fprintf('\n------------------------------\n')
    fprintf('Image %d:\n\n', i)
    fprintf('NDF:\nTotal: %.4f\nMean: %.4f\nSEM: %.4f\n\n', results(i).TotalNDF, results(i).MeanNDF, results(i).SterrNDF)
    fprintf('Halo Moment:\nMean: %.4g\nSEM: %.4g\n\n', results(i).MeanHaloMoment, results(i).SterrHaloMoment)
    fprintf('Adjusted Halo Moment:\nMean: %.4g\nSEM: %.4g\n\n', results(i).MeanAdjHaloMom, results(i).SterrAdjHaloMom);
    fprintf('Halo Distance:\nMean: %.4g\nSEM: %.4g\n', results(i).MeanHaloDist, results(i).SterrHaloDist);
end