close all
clear all
clc

imArray = importImage();
% figure(), imshow(imArray{1}), title('Original Image')

tic

imRegions = processImage(imArray, false);
% figure(), imshow(label2rgb(imRegions(1).Thresh)), title('Thresholded Image')
% figure(), imshow(imRegions(1).Binary), title('Binary Image')
% figure(), imshow(imRegions(1).Nuclei), title('Nuclei Passing Eccentricity Criteria')
% figure(), imshow(imRegions(1).Halos),  title('Halos Passing Eccentricity Criteria')

pairs = findMatches(imRegions, false);

results = calcDamage(imRegions, pairs, 'NDF', 'moment', 'adjMoment', 'distance', false);

toc

% for i = 1:length(results)
%     fprintf('\n------------------------------\n')
%     fprintf('Image %d:\n\n', i)
%     fprintf('NDF:\nTotal: %.4f\nMean: %.4f\nSEM: %.4f\n\n', results(i).TotalNDF, results(i).MeanNDF, results(i).SterrNDF)
%     fprintf('Halo Moment:\nMean: %.4g\nSEM: %.4g\n\n', results(i).MeanHaloMoment, results(i).SterrHaloMoment)
%     fprintf('Adjusted Halo Moment:\nMean: %.4g\nSEM: %.4g\n\n', results(i).MeanAdjHaloMom, results(i).SterrAdjHaloMom);
%     fprintf('Halo Distance:\nMean: %.4g\nSEM: %.4g\n', results(i).MeanHaloDist, results(i).SterrHaloDist);
% end





