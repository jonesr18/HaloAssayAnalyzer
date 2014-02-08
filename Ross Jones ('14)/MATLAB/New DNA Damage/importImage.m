% importImage runs the GUI to import one or multiple images
%
%   imstruct = importImage()
%
%       No arguments are accepted as input, the user can select one or more images from the same
%       directory using the GUI that pops up and both cases are handled by the function.
%
%   Output:
%
%       The function returns a cell array the size of the number of files selcted with each element
%       of the array being an image that is pre-processed with an averaging filter. 
%
%   Example:
% 
%       % Call function to select images, plot the first one
%       image = importImage();
%       imshow(image{1}), title('First image')
%
% 
% Ross Jones
% Singh Laboratory
% University of Washington
% Updated 1.17.14

function imArray = importImage()
    
    % Run GUI to select file(s)
    [files, filePath] = uigetfile('*.jpeg; *.png; *.tif', 'Pick image file(s)', 'MultiSelect', 'on');
    
    if iscell(files)
        numFiles = numel(files);
        imArray = cell(numFiles, 1);
        for i = 1:numFiles
            image_i = rgb2gray(imread(fullfile(filePath, files{i})));
            imArray{i} = imfilter(image_i, fspecial('average'));
        end
    else
        image = rgb2gray(imread(strcat(filePath, files)));
        imArray = {imfilter(image, fspecial('average'))};
    end
end

