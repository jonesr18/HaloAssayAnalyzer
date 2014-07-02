% FuncDef
%   Details
function [] = gradientSegment(imRegions)
    
    image = imRegions.Image;
    
    % Edge detection of nuclei for segmentations
    close all

    % Get gradient magnitude
    hy = fspecial('sobel');
    Iy = imfilter(double(image), hy, 'replicate');
    Ix = imfilter(double(image), hy', 'replicate');
    gradmag = sqrt(Ix.^2 + Iy.^2);
    GRADMAG = gradmag;
    figure(), imshow(gradmag, [min(gradmag(:)), max(gradmag(:)/4)]), title('Gradient Magnitude')

    % Threshold gradient magnitude
    imthresh = imquantize(gradmag, multithresh(gradmag, 2));
    gradmag  = boolean(imthresh > 2);
    figure(), imshow(gradmag), title('Thresholded Gradient Magnitude')

    % Dilate the image
    strelDisk = strel('disk', 1);
    strelDiam = strel('diamond', 3);
    gradmag = imopen(gradmag, strelDisk);
    gradmag = imdilate(gradmag, strelDiam);
    gradmag = imdilate(gradmag, strelDiam);
    gradmag = imopen(gradmag, strelDisk);
    gradmag = imclose(gradmag, strelDisk);
    gradmag = imopen(gradmag, strelDisk);
    gradmag = imclose(gradmag, strelDisk); % closing and opening smoothes edges
    figure(), imshow(gradmag), title('Dilated Gradient Mask')

    % Fil interior gaps
    gradmag = imfill(gradmag, 'holes');
    figure(), imshow(gradmag), title('Binary Image with Filled Holes')

    % Remove connected objects on border
    gradmag = imclearborder(gradmag, 4);
    figure(), imshow(gradmag), title('Segmented Image')

    % Remove small objects and outline remaining ones
    [boundaries, labels] = bwboundaries(gradmag);
    figure(), imshow(gradmag), title('Outlined Original Image'), hold on
    stats = regionprops(labels, 'Area');
    for i = 1:length(stats)
        if stats(i).Area < 300
            gradmag = gradmag - boolean(labels == i);
        else
            boundary = boundaries{i};
            plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2)
        end
    end
    hold off

    binGradmag = boolean(gradmag);
    statsNew = regionprops(binGradmag, 'EquivDiameter');
    radList = [statsNew.EquivDiameter] / 2;
    rad = round(mean(radList));
    radstd = round(std(radList));

    % Fit circles to the image to find individual nuclei
    [centroids, radii] = imfindcircles(binGradmag, [rad - radstd, rad + radstd], 'ObjectPolarity', 'bright', 'Sensitivity', 0.92);

    % Match pixels to closest centroid
    blank = zeros(size(binary));
    for c = centroids'
        x = round(c(1));
        y = round(c(2));
        blank(y, x) = 1;
    end
    labels = watershed(bwdist(blank)) .* uint8(imRegions.Binary);

    % Plot estimated circles on image
    low  = min(GRADMAG(:));
    high = max(GRADMAG(:)/4);
    figure(), imshow(GRADMAG, [low, high]), title('Gradmag + Circles'), viscircles(centroids, radii);
    figure(), imshow(image, []), title('Original Image + Circles'), viscircles(centroids, radii);
    figure(), imshow(label2rgb(labels, 'jet', [.5 .5 .5])), title('Labelled Cells')
end
