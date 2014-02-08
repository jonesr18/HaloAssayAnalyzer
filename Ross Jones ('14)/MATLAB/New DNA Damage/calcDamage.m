% calcDamage calculates the DNA damage using various metrics
%
%   results = calcDamage(imRegions, pairs)
%
%       The function accepts a struct array and a cell array as input. For all tests to be performed
%       correctly, the first input (the struct) must at least contain the fields 'Image', 'Nuclei',
%       and 'Halos' (case-insensitive). The struct array can be any size, each element representing 
%       a single image. The second input (the cell array) must have a size of Nx3, with the first 
%       element being a boolean and the second and third elements being numeric arrays. The second
%       element must be an array of indexes for corresponding halo or nuclei pairs. 
%
%   results = calcDamage(imRegions, pairs, OPTIONS)
%
%       The default behavior of the function is to calculate just the nuclear diffusion factor
%       (NDF), and not to plot any of the DNA damage metrics on the original greyscale image. This
%       behavior can be changed by adding arguments after the second input. Available options:
%
%           'all'       All DNA damage metrics are turned on
%           'NDF'       Calculate the nuclear diffusion factor metric            (defualt = on)
%           'moment'    Calculate the halo moment metric                         (default = off)
%           'adjMoment' Calculate the adjusted halo moment metric                (default = off)
%           'distance'  Calculate the halo distance metric                       (default = off)
%           true/false  Set plotting behavior                                    (default = false)
%
%       These options can be input in any order and are case-insensitive. The default of NDF being
%       turned on is reversed if any other tests (beside 'all') are passed as arguments, and thus 
%       must be manually turned back on by passing 'NDF' as another argument. 
%
%   Output:
%
%       The funtion returns a struct array where each element in the array corresponds with an image
%       that is given in the first input argument struct array. The fields in this array give the
%       data for the metrics selected when the function is called.
%
%   Fields created:
%
%       'TotalNDF'        The total NDF not calculated on a cell-by-cell basis   (defualt = on)
%       'MeanNDF'         The mean NDF calculated for each cell                  (defualt = on)
%       'SterrNDF'        The standard error of the mean (SEM) for the mean NDF  (defualt = on)
%       'MeanHaloMoment'  The mean halo moment calculated for each cell          (defualt = off)
%       'SterrHaloMoment' The SEM for the mean halo moment                       (defualt = off)
%       'MeanAdjHaloMom'  The mean adjusted halo moment calculated for each cell (defualt = off)
%       'SterrAdjHaloMom' The SEM for the mean adjusted halo moment              (defualt = off)
%       'MeanHaloDist'    The mean halo distance calculated for each cell        (defualt = off)
%       'SterrHaloDist'   The SEM for the mean halo distance                     (defualt = off)
%
%   Damage metrics:
%
%       NDF:
%
%           The NDF is the ratio of the combined area of a cell's halo and nucleus to the nucleus
%           alone, and is used in most publications assaying with alkaline halo assays. 
%
%       Halo moment:
%
%           The halo moment is calculated by finding the sum of the intensity of each pixel in the
%           halo multiplied by its distance from the centroid of the nucleus.
%
%       Adjusted halo moment:
%
%           This metric adjusts the halo moment by dividing by the total number of cells and the
%           radius in pixels of the nucleus. This allows for comparisons between images of differing
%           resolutions and for normalization to the spread of the nuclear DNA, should some cells
%           with similar halo profiles have disparate nuclei.
%
%       Halo distance:
%
%           This is essentially the adjusted halo moment, but without the intensity calculation.
%           This allows one to understand how far the average pixel in the halo is from the nuclear
%           centroid in terms of nuclear radii. Thus, it is a normalized metric purely for
%           understanding the level of DNA spread for a cell.
%
%       Example:
%
%           % Given an image and its halos and nuclei that are paired, calculate halo moment
%           imRegions = struct('Image', image, 'Halos', halos, 'Nuclei', nuclei);
%           pairs = {moreHalos, indexes, distances}; % cell array
%           results = calcDamage(imRegions, pairs, 'moment');
%           average = results(1).MeanHaloMoment;
%           sterr   = results(1).SterrHaloMoment;
%           fprintf('Halo Moment:\nMean: %.4g\nSEM: %.4g\n\n', average, sterr)
%
%
% Ross Jones
% Singh Laboratory
% University of Washington
% Updated 1.28.14

function results = calcDamage(imRegions, pairs, varargin)
    
    % Check inputs and assign value to getNDF, getMoment, getAdjMoment, and getDistance
    % (defaults = [true, false, false, false])
    [showStats, getNDF, getMoment, getAdjMom, getDist] = parseInputs(imRegions, pairs, varargin);
    
    % Initialize results struct array
    results = struct;
    
    for n = 1:length(imRegions)
        
        nuclei = imRegions(n).Nuclei;
        halos  = imRegions(n).Halos;
        if getMoment || getAdjMom
            image = imRegions(n).Image;
        end
        
        % Extract relevant information from nuclei and halos
        nucStats  = regionprops(nuclei, 'Area', 'PixelList', 'Centroid', 'EquivDiameter');
        haloStats = regionprops(halos, 'Area', 'PixelList', 'Centroid', 'EquivDiameter');
        
        % Extract relevant information from pairs cell array
        moreHalos = pairs{n, 1};
        pair = pairs{n, 2};
        if moreHalos
            numPairs = length(nucStats);
        else
            numPairs = length(haloStats);
        end
        
        % Calculate NDF for each per halo-nucleus pair -- O(n)
        if getNDF
            damage = zeros(numPairs, 1);
            if moreHalos
                for i = 1:numPairs
                    damage(i) = haloStats(pair(i)).Area / nucStats(i).Area;
                end
            else
                for i = 1:numPairs
                    damage(i) = haloStats(i).Area / nucStats(pair(i)).Area;
                end
            end
            results(n).MeanNDF  = mean(damage);
            results(n).SterrNDF = std(damage) / sqrt(numPairs);
            
            % Compute statistics for entire iamge simultaneously for comparison
            threshStats = regionprops(imRegions(n).Thresh, 'Area');
            nucSurface  = threshStats(3).Area;
            haloSurface = threshStats(2).Area;
            results(n).TotalNDF = (haloSurface + nucSurface) / nucSurface;  
        end
        
        % These tests all rely on many of the same image manipulations
        if getMoment || getAdjMom || getDist
                        
            % Collect pixel coordinates that are only in halo region -- O(n^2)
            % Conisder doing this by syubtracting images instead -- should be much faster
            haloPixels = struct;
            for i = 1:numPairs
                if moreHalos
                    allPixels = haloStats(pair(i)).PixelList;
                    nucPixels = nucStats(i).PixelList;
                else
                    allPixels = haloStats(i).PixelList;
                    nucPixels = nucStats(pair(i)).PixelList;
                end
                remainingPixels = zeros(length(allPixels) - length(nucPixels), 2);
                overlaps = ismember(allPixels, nucPixels, 'rows');
                k = 1;
                for j = 1:length(allPixels)
                    if ~overlaps(j)
                        remainingPixels(k, :) = allPixels(j, :);
                        k = k + 1;
                    end
                end
                haloPixels(i).PixelList = remainingPixels;
                haloPixels(i).Overlaps = overlaps;
            end
            
            % Calculate halo moment, adjusted halo moment, and/or halo distance -- O(n^2)
            if getMoment
                haloMoments = zeros(numPairs, 1);
            end
            if getAdjMom
                adjHaloMom = zeros(numPairs, 1);
            end
            if getDist
                haloDists = zeros(numPairs, 1);
            end
            for i = 1:numPairs
                cn        = nucStats(i).Centroid;
                coords    = haloPixels(i).PixelList;
                numPixels = size(coords, 2);
                if getMoment || getAdjMom
                    moments = zeros(numPixels, 1);
                end
                if getDist
                    distances = zeros(numPixels, 1);
                end
                % Need this try block in case there are no halo coordinates outside of nucleus
                try
                    % Calculate distances from centroid and intensities
                    for j = 1:numPixels
                        x = cn(1) - coords(j, 1);   % x distance
                        y = cn(2) - coords(j, 2);   % y distance
                        distance   = sqrt(x.^2 + y.^2); 
                        intensity  = double(image(coords(j, 2), coords(j, 1)));
                        if getMoment || getAdjMom
                            moments(j) = intensity * distance;
                        end
                        if getDist
                            distances(j) = distance;
                        end
                    end
                    % Finish calculations by summing and performing other operations
                    if getMoment
                        haloMoments(i) = sum(moments);
                    end
                    if getAdjMom || getDist
                        nucRadius = nucStats(i).EquivDiameter / 2;
                        if getAdjMom
                            adjHaloMom(i) = mean(moments) / nucRadius;
                        end
                        if getDist
                            haloDists(i) = mean(distances) / nucRadius;
                        end
                    end
                catch
                    % Execute these commands if an error would occur
                    if getMoment
                        haloMoments(i) = 0;
                    end
                    if getAdjMom
                        adjHaloMom(i) = 0;
                    end
                    if getDist
                        haloDists(i) = 1;
                    end
                    if moreHalos
                        num = pair(i);
                    else
                        num = i;
                    end
                    warning('Halo %d in image %d has no pixels outside of the nucleus', num, n)
                end
            end
            % Store results for each image
            if getMoment
                results(n).MeanHaloMoment  = mean(haloMoments);
                results(n).SterrHaloMoment = std(haloMoments);
            end
            if getAdjMom
                results(n).MeanAdjHaloMom  = mean(adjHaloMom);
                results(n).SterrAdjHaloMom = std(adjHaloMom) / sqrt(numPairs);
            end
            if getDist
                results(n).MeanHaloDist  = mean(haloDists);
                results(n).SterrHaloDist = std(haloDists) / sqrt(numPairs);
            end
        end
        
        % Display damage and halo moment next to each nucleus -- O(n)
        if showStats
            figure(), imshow(image), title('Data Metrics for Each Measured Cell')
            hold on
            for i = 1:numPairs
                cn = nucStats(i).Centroid;
                shift = 0;
                if getNDF
                    text(cn(1), cn(2), num2str(damage(i)), 'Color', 'g')
                    text(0, 11, 'NDF', 'Color', 'g')
                    shift = shift + 22;
                end
                if getMoment
                    text(cn(1), cn(2) + shift, sprintf('%.5g', haloMoments(i)), 'Color', 'y')
                    text(0, 11 + shift, 'Moment', 'Color', 'y')
                    shift = shift + 22;
                end
                if getAdjMom
                    text(cn(1), cn(2) + shift, sprintf('%.5g', adjHaloMom(i)), 'Color', 'c')
                    text(0, 11 + shift, 'Adj Moment', 'Color', 'c')
                    shift = shift + 22;
                end
                if getDist
                    text(cn(1), cn(2) + shift, sprintf('%.5g', haloDists(i)), 'Color', 'm')
                    text(0, 11 + shift, 'Distance', 'Color', 'm')
                end
            end
            hold off
        end      
    end
end


% Check if the given inputs are valid, returns booleans showStats, getNDF, getMoment, getAdjMom,
% and getDist for use in the main function.
function [showStats, getNDF, getMoment, getAdjMom, getDist] = parseInputs(imRegions, pairs, args)
    
    % Default Values
    showStats = false;
    getNDF    = true;
    getMoment = false;
    getAdjMom = false;
    getDist   = false;
    
    first = true;
    if nargin > 1
        if nargin > 2
            for i = 1:numel(args)
                if islogical(args{i})
                    showStats = args{i};
                elseif ischar(args{i})
                    if first
                        getNDF = false;
                        first  = false;
                    end
                    if strcmpi(args{i}, 'all')
                        getNDF    = true;
                        getMoment = true;
                        getAdjMom = true;
                        getDist   = true;
                        break
                    elseif strcmpi(args{i}, 'NDF')
                        getNDF = true;
                    elseif strcmpi(args{i}, 'moment')
                        getMoment = true;
                    elseif strcmpi(args{i}, 'adjMoment')
                        getAdjMom = true;
                    elseif strcmpi(args{i}, 'distance')
                        getDist = true;
                    else
                        error('Capstone:ArgChk', 'Test name not recognized: %s', args{i})
                    end
                else
                    error('Capston:ArgChk', 'Test names must be of type string.')
                end
            end
            if ~any([getNDF, getMoment, getAdjMom, getDist])
                error('Capstone:ArgChk', 'At least one test must be selected')
            end
        end
        if isstruct(imRegions)
            if ~any(strcmpi(fieldnames(imRegions), 'Halos'))
                error('Capstone:ArgChk', 'Struct input argument must contain field: ''Halos''')
            end
            if ~any(strcmpi(fieldnames(imRegions), 'Nuclei'))
                error('Capstone:ArgChk', 'Struct input argument must contain field: ''Nuclei''')
            end
            if getNDF
                if ~any(strcmpi(fieldnames(imRegions), 'Thresh'))
                    error('Capstone:ArgChk', 'Struct input argumetn must contain field: ''Thresh''')
                end 
            end
            if getMoment || getAdjMom
                if ~any(strcmpi(fieldnames(imRegions), 'Image'))
                    error('Capstone:ArgChk', 'Struct input argumetn must contain field: ''Image''')
                end
            end
        else
            error('Capstone:ArgChk', 'First input argument must be of type struct')
        end
        if iscell(pairs)
            if size(pairs, 2) ~= 3
                error('Capstone:ArgChk', 'Cell array ''pairs'' must contain 3 elements')
            end
            for i = 1:size(pairs, 1)
                if ~islogical(pairs{i, 1})
                    error('Capstone:ArgChk', 'First element in each cell row must be of type boolean')
                end
                if ~isnumeric(pairs{i, 2})
                    error('Capstone:ArgChk', 'Second element in each cell row must be of type numeric')
                end
                if ~isnumeric(pairs{i, 3})
                    error('Capstone:ArgChk', 'Third element in each cell row must be of type numeric')
                end
            end
        else
            error('Capstone:ArgChk', 'Second input argument must be of type cell')
        end
    else
        error('Capstone:ArgChk', 'Not enough input arguments.')
    end
    
end