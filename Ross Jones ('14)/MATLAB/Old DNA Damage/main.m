% Ross Jones
% Singh Lab
% University of Washington
% Main file: image import and algorithm calling

% PROBLEM: random upside down image of a child appears...?

function [] = main()
    close all
    clear all
    clc
    
    % Select and import files
    % PROBLEM: doesn't work for just one image
    [files, filePath] = uigetfile('*.jpeg; *.png; *.tif', 'Pick image file(s)', 'MultiSelect', 'on');
    numFiles = numel(files);

    % Preallocate the arrays
    image    = rgb2gray(imread(fullfile(filePath, files{1})));
    sequence = zeros([size(image) numFiles], class(image));

    % Create image sequence array
    % PROBLEM: doesn't work with different sized images
    for i = 1:numFiles
        image_i = rgb2gray(imread(fullfile(filePath, files{i}))); 
        sequence(:, :, i) = image_i;
    end
    
    % Calculate and display DNA damage using multithresh
    % NOTE: This will return two measurements for all images combined and is a fast algorthim
    tic
    DNA_damage_multithresh(sequence, numFiles)
    toc
    
    % Calculate and display DNA damage using kmeans clustering
    % NOTE: This will return two measurements for all images combined and is a slow algorithm
    tic
%     DNA_damage_kmeans(files, filePath, numFiles, sequence)
    toc
    
    % Find circles in images
    % NOTE: This will show 2x the number of images selected
    tic
    Optimized_circle_finder(sequence, numFiles)
    toc
    
    % Measure roundness of each object in image
    % NOTE: This will show 2x the number of images selected
    tic
    Measure_roundness(sequence, numFiles)
    toc
end

















