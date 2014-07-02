% Ross Jones
% Singh Lab
% University of Washington
% Object roundness measurement

% From MATLAB examples online: 
% http://www.mathworks.com/help/images/examples/identifying-round-objects.html?prodcode=IP&language=en

% After evaluation, this method will not work, as it meshes together nearby and touching cells,
% causing them to appear to be one shape, and thus non-circular.

function [] = Measure_roundness(sequence, numFiles)

    for i = 1:numFiles

        % Create binary image
        image_i = sequence(:, :, i);
        threshold = graythresh(image_i);
        binary = im2bw(image_i, threshold);

        % Remove the noise
        area = 30;                              % Size of areas to be removed
        binary = bwareaopen(binary, area);
        strucElem = strel('disk', 2);           % Make structuring element to fill gaps
        binary = imfill(binary, 'holes');       % Fill in holes

        % Find the boundaries
        [boundaries, labels] = bwboundaries(binary, 'noholes');
        figure(), imshow(image_i), title('Original Image')
        figure(), imshow(binary), title('All Boundaries'), hold on
        for j = 1:length(boundaries)
            boundary = boundaries{j};
            plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2)
        end

        % Determine which objects are round
        stats = regionprops(labels, 'Area', 'Centroid');
        threshold = 0.80;
        for j = 1:length(boundaries)
            % obtain (X,Y) boundary coordinates corresponding to label 'k'
            boundary = boundaries{j};

            % compute a simple estimate of the object's perimeter
            delta_sq  = diff(boundary).^2;
            perimeter = sum(sqrt(sum(delta_sq, 2)));

            % obtain the area calculation corresponding to label 'k'
            area = stats(j).Area;

            % compute the roundness metric
            metric = 4 * pi * area / perimeter^2;

            % display the results
            metricString = sprintf('%2.2f', metric);

            % mark objects above the threshold with a black circle
            if metric > threshold
                centroid = stats(j).Centroid;
                plot(centroid(1), centroid(2), 'ko');
            end

            text(boundary(1, 2)-35, boundary(1, 1)+13, metricString, 'Color', 'y', 'FontWeight', 'bold');
        end
    end
end




