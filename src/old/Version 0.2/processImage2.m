% processImage2 identifies unclustered cells in the given image(s)
%   
%   imRegions = processImage2(imArray)
%
%       The function accepts a cell array of one or more grayscale images as input. 
%
%   imRegions = processImage2(imArray, showMetrics)
%   
%       Plotting behavior can be set by passing a boolean as the second argument:
%
%           true:  Plot each image in the given struct with boundaries and metrics.
%           false: Do not plot images (default given no input).
%
%   Output:
%
%       The function returns a Java ArrayList the size of the number of given images with each 
%       element containing a Java HashMap with certain keys. In order to improve the accuracy of DNA 
%       damage calculations, overlapping cells are removed by measuring the elliptical properties of
%       the cell objects. Highly elliptical binary regions are considered to be clustered cells. 
%
%   Keys created:
%       
%       'Image'        Original grayscale image passed into function
%       'Thresh'       3-level threshold of original image
%       'Binary'       Binary image constructed from the two thresholds comprising cell object
%       'Boundaries'   A cell array of boundary coordinates for each object in the binary image
%       'Eccentricity' An array of the eccentricity of each object in the binary image
%       'Halos'        The halos of non-overlapping cells (effectively a non-overlapping binary)
%       'Nuclei'       The nuclei of non-overlapping cells
%
%   Example:
%
%       % Display the thresholded image
%       image     = rgb2gray(imread('random.jpg'));
%       imArray   = {image};
%       imRegions = processImage(imArray);
%       imshow(label2rgb(imRegions.get(0).get('Thresh'))), title('Thresholded Image')
%
%
% Ross Jones
% Singh Laboratory
% University of Washington
% Updated 2.12.14

function imRegions = processImage2(imArray, varargin)
    
    % Check inputs and assign value to showMetrics (default = false)
    showMetrics = parseInputs(imArray, varargin);
    
    % Initialize return ArrayList and other variables used in the loop
    imRegions  = java.util.ArrayList(numel(imArray));
    strelDiam  = strel('diamond', 1);
    thresholds = 2;
    
    for n = 1:numel(imArray)
        imRegions.add(java.util.HashMap());
        
        % Threshold with multithresh, then combine two thresholds to make better binary image
        image    = imArray{n};
        thresh   = multithresh(image, thresholds);
        imThresh = uint8(imquantize(image, thresh));
        binary   = boolean(imThresh > 1);
        imRegions.get(n - 1).put('Image', image);
        imRegions.get(n - 1).put('Thresh', imThresh);
        imRegions.get(n - 1).put('Binary', binary);
        
        % Smoothen the image objects
        binary = imerode(binary, strelDiam);
        binary = imerode(binary, strelDiam);
        
        % Find the boundaries (returned in a cell array for each object)
        [boundaries, labels] = bwboundaries(binary, 'noholes');
        imRegions.get(n - 1).put('Boundaries', boundaries);
        
        % If user wants metrics to be show, an image will appear with boundaries drawn around each
        % cell, the eccentricity measure show next to each cell, and circles at the center of cells
        % which pass the criteria. 
        if showMetrics
            figure(), imshow(image), title('Cells Passing Eccentricity Criteria Marked ''o''')
            hold on
        end
        
        % Determine which objects are round, which ideally correlates with just one cell.
        % Eccentricity ranges between 0 (a line) and 1 (a circle). 
        stats = regionprops(labels, 'Area', 'Centroid', 'Eccentricity');
        eccentricity = zeros(length(stats), 1);
        threshold = 0.500;
        for i = 1:length(stats)
            eccentricity(i) = stats(i).Eccentricity;
            isObject = eccentricity(i) < threshold && stats(i).Area > 100;
            
            % Plot boundaries, metrics, and circles as noted above
            if showMetrics
                boundary = boundaries{i};
                eccenString = sprintf('%.3f', eccentricity(i));
                plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2)
                text(boundary(1, 2) - 35, boundary(1, 1) + 50, eccenString, 'Color', 'g');
                if isObject
                    centroid = stats(i).Centroid;
                    plot(centroid(1), centroid(2), 'ko');
                end
            end
            
            % Remove objects from binary image that do not pass criteria
            if ~isObject
                binary = binary - boolean(labels == i);
            end
        end
        imRegions.get(n - 1).put('Eccentricity', eccentricity);
        
        % Remove regions from thresholded image --- note that removing from greyscale and 
        % rethresholding can cause issues, as sometimes there become dicontinginous regions. 
        imThresh = imThresh .* uint8(binary);
        
        % Obtain binary images containing all halos or all nuclei (note: halos include nuclei)
        imRegions.get(n - 1).put('Halos', boolean(imThresh > 1));
        imRegions.get(n - 1).put('Nuclei', boolean(imThresh > thresholds));
    end
end


% Check if the given inputs are valid, returns boolean showMetrics for use in the main function.
function showMetrics = parseInputs(imArray, args)
    
    % Default values
    showMetrics = false; 
    
    % Check inputs
    if nargin > 0
        if ~isempty(args)
            if nargin > 2
                error('Capstone:ArgChk', 'Function only accepts two input arguments.')
            end
            if islogical(args{1})
                showMetrics = args{1};
            else
                error('Capstone:ArgChk', 'Second input argument must be of type boolean')
            end
        end
        if ~iscell(imArray)
            error('Capstone:ArgChk', 'First input argument must be of type cell')
        end
    else
        error('Capstone:ArgChk', 'Not enough input arguments.')
    end
    
end
