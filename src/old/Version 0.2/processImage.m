% processImage identifies unclustered cells in the given image(s)
%   
%   imRegions = processImage(imArray)
%
%       The function accepts either a cell array of one or more images, or just a single image as 
%       the first input. 
%
%   imRegions = processImage(imArray, OPTIONS)
%   
%       The default behavior of the function is to not plot each image with object boundaries and
%       eccentricity metrics overlaid, apply 2 thresholds to the image to achieve 3 levels, and use 
%       a sensitivity of 0.500 for the processing of overlapped cells. This behavior can be changed
%       by adding optional arguments after the first input.The sensitivity must be a number between
%       0 and 1. A larger number will incorporate more overlapping cells and a smaller number will
%       reject more cells based off of an eccentricity metric. 
%       
%           true/false      Set plotting behavior                                (default = false)
%           number >= 1     Set the number of thresholds used                    (default = 2)
%           0 < number < 1  Set sensitivity of overlapping cell removal          (default = 0.500)
%
%   Output:
%
%       The function returns a struct the size of the number of given images with several fields and
%       a cell array of CellRegion objects are returned which constitute each image. In order to
%       improve the accuracy of DNA damage calculations, overlapping cells are removed by measuring
%       the elliptical properties of the cell objects. Highly elliptical binary regions are
%       considered to be clustered cells. 
%
%   Fields created:
%       
%       'Image'        Original grayscale image passed into function
%       'Thresh'       3-level threshold of original image
%       'Binary'       Binary image constructed from the two thresholds comprising cell object
%       'Halos'        The halos of non-overlapping cells (effectively a non-overlapping binary)
%       'Nuclei'       The nuclei of non-overlapping cells
%
%   Example:
%
%       % Display the thresholded image
%       image     = rgb2gray(imread('random.jpg'));
%       imArray   = {image};
%       imRegions = processImage(imArray);
%       imshow(label2rgb(imRegions(1).Thresh)), title('Thresholded Image')
%
%
% Ross Jones
% Singh Laboratory
% University of Washington
% Updated 2.20.14

function imRegions = processImage(imArray, varargin)
    
    % Check inputs and assign value to showMetrics (default = false)
    [showMetrics, thresholds, sensitivity] = parseInputs(imArray, varargin{:});
    if iscell(imArray)
        numImages = numel(imArray);
    else
        numImages = 1;
    end
    
    % Initialize return arrays and other variables used in the loop
    imRegions = struct();
    strelDisk = strel('disk', 20);
    
    % If the input is a single RGB image, convert it to greyscale
    if ~iscell(imArray)
        if size(imArray, 3) > 1 % multiple dimensions (r, g, b)
            imArray = rgb2gray(imArray);
        end
    end
    
    for n = 1:numImages
        % Exctract image and ensure it is grayscale
        if iscell(imArray)
            image = imArray{n};
            if size(image, 3) > 1 % multiple dimensions (r, g, b)
                image = rgb2gray(image);
            end
        else
            image = imArray;
        end
        
        % Threshold with multithresh, then combine three thresholds to make better binary image
        thresh   = multithresh(image, thresholds);
        imThresh = uint8(imquantize(image, thresh));
        binary   = boolean(imThresh > 1);
        imRegions(n).Image  = image;
        imRegions(n).Thresh = imThresh;
        imRegions(n).Binary = binary;
        
        % Smoothen the image objects
        binary = imclose(binary, strelDisk);
        binary = imopen(binary, strelDisk);
        
        % Find the boundaries (returned in a cell array for each object)
        [boundaries, labels] = bwboundaries(binary, 'noholes');
        
        % If user wants metrics to be shown, an image will appear with boundaries drawn around each
        % cell, the eccentricity measure show next to each cell, and circles at the center of cells
        % which pass the criteria. 
        if showMetrics
            figure(), imshow(image, []), title('Cells Passing Eccentricity Criteria Marked ''o''')
            hold on
        end
        
        % Determine which objects pass the given eccentricity criteria.
        % Eccentricity ranges between 0 (a line) and 1 (a circle). 
        stats      = regionprops(labels, 'Area', 'Centroid', 'Eccentricity');
        numRegions = length(stats);
        for i = 1:numRegions
            eccentricity = stats(i).Eccentricity;
            isObject = eccentricity < sensitivity && stats(i).Area > 100;
            
            % Plot boundaries, metrics, and circles as noted above
            if showMetrics
                boundary = boundaries{i};
                eccenString = sprintf('%.3f', eccentricity);
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
        
        % Remove regions from thresholded image --- note that removing from greyscale and 
        % rethresholding can cause issues, as sometimes there become dicontinginous regions. 
        imThresh = imThresh .* uint8(binary);
        
        % Obtain binary images containing all halos or all nuclei
        nuclei = boolean(imThresh > thresholds);
        imRegions(n).Nuclei = nuclei;
        imRegions(n).Halos = boolean(imThresh > 1) - nuclei;
    end
end


% Check if the given inputs are valid, returns boolean showMetrics for use in the main function.
function [showMetrics, thresholds, sensitivity]  = parseInputs(varargin)
    
    % Default values
    showMetrics = false; 
    thresholds  = 2;
    sensitivity = 0.500;
    
    % First input must be an RGB image, a greyscale image, or a cell array of either
    validateattributes(varargin{1}, {'uint8', 'cell'}, {}, mfilename, 'imArray', 1);
    if iscell(varargin{1})
        for im = varargin{1}
            validateattributes(im{:}, {'uint8'}, {}, mfilename, 'cell array', 1)
        end
    end
    
    % Second, third...etc inputs must be an number > 0 && < 1, an integer >= 1, or a boolean
    for i = 2:nargin
        arg = varargin{i};
        validateattributes(arg, {'logical', 'numeric'}, {'nonnegative'}, mfilename, 'OPTIONS', i);
        if islogical(arg)
            showMetrics = arg;
        elseif arg < 1
            sensitivity = arg;
        else
            validateattributes(arg, {'numeric'}, {'integer'}, mfilename, 'OPTIONS', i)
        end
    end
end
