% Ross Jones
% Singh Lab
% University of Washington
% optimized circle detection using imfindcircles

% From MATLAB examples online:
% http://www.mathworks.com/help/images/examples/detect-and-measure-circular-objects-in-an-image.html?prodcode=IP&language=en

% This seems to have identified virtually all of the healthy cells in the image, including
% overalapping ones. It did not identify apoptotic cells. It may be overdoing it on identifying
% two cells when really it is just an elongated cell

function [] = Optimized_circle_finder(sequence, numFiles)

    tic
    for i = 1:numFiles
        image_i = sequence(:, :, i);
        strDisk = strel('disk', 5);
        image_i = imopen(image_i, strDisk);
        image_i = imclose(image_i, strDisk);
        
        % Initial guesses
        radLow = 30;
        radUp  = 50;
        delta  = 1;
        [centers, radii] = imfindcircles(image_i, [radLow radUp], 'ObjectPolarity', 'bright', 'Sensitivity', 0.92);
        
        % Find cirlces in each image
%         while delta > .01
%             
%             % New guesses
%             radAvg = mean(radii);
%             radUp  = round(radAvg + 5);
%             radLow = round(radAvg - 5);
%             
%             % Optimize range of radii to be searched for cirlces
%             [centers, radNew] = imfindcircles(image_i, [radLow radUp], 'ObjectPolarity', 'bright', 'Sensitivity', 0.92);
%             radNewAvg = mean(radNew);
%             
%             % Calculate new values
%             delta  = abs(radAvg - radNewAvg);
%             radii  = radNew;
%         end
        
        fprintf('Mean: %f\nStdev: %f\nLower bound: %i\nUpper bound: %i\n\n', mean(radii), std(radii), radLow, radUp)
        figure(), imshow(image_i), title('Original Image')
        figure(), imshow(image_i), title('Optimized Circles'), viscircles(centers, radii);
    end
    toc
end





