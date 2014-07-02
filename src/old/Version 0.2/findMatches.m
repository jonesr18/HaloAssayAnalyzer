% findMatches matches each halo and nucleus object in one or multiple images
%
%   pairs = findMatches(imRegions)
%
%       The fucntion accepts a struct contianing binary images of cellular halos and nuclei. The
%       given struct must contain the fields 'Halos', 'Nuclei', and 'Image' (case-insensitive) in 
%       order for this to work properly. The default behavoir is to not plot the matches, but a 
%       warning will be shown to the user mismatches are suspected.
%
%   pairs = findMatches(imRegions, debug)
%
%       Plotting behavior can be set by passing a boolean as the second argument:
%
%           true:  Plot each image in the given struct with nuclear and halo centorids labeled
%           false: Do not plot images
%
%   Output:
%
%       The function returns a cell array with three columns (one row per image). For a given row in
%       the cell array, the first column is a boolean that is true if there were more halos than 
%       nuclei, the second column contains the indexes of the paired items in the corresponding 
%       image, and the third column contains the distances between the halo and nuclear centroids
%       forming these pairs (distance is given as city block distance). The indexes of the paired
%       items are used for the halos if the boolean is true and for the nuclei otherwise.
%
%   Example:
%
%       % Given an image and its halos and nuclei:
%       imRegions = struct('Image', image, 'Halos', halos, 'Nuclei', nuclei);
%       pairs = findMatches(imRegions, true); 
%       moreHalos = pairs{1, 1};
%       indexes   = pairs{1, 2};
%       distances = pairs{1, 3};
%       if moreHalos
%           haloIndexes = indexes;
%       else
%           nucIndexes = indexes;
%       end
%
% Ross Jones
% Singh Laboratory
% University of Washington
% Updated 2.12.14

function pairs = findMatches(imRegions, varargin)
    
    % Check inputs and assign value to debug (default = false)
    debug = parseInputs(imRegions, varargin);
    
    % Initialize return cell array
    pairs = cell(length(imRegions), 3);
    
    for n = 1:length(imRegions)
        
        halos  = imRegions(n).Halos;
        nuclei = imRegions(n).Nuclei;
        
        % Find centroid of each nucleus and halo
        haloStats = regionprops(halos, 'Centroid', 'Area', 'PixelList');
        nucStats  = regionprops(nuclei, 'Centroid', 'Area', 'PixelList', 'EquivDiameter');
        
        % Compare number of halos to nuclei
        moreHalos   = length(haloStats) > length(nucStats); 
        pairs{n, 1} = moreHalos;
        if moreHalos
            highNum = length(haloStats);
            lowNum  = length(nucStats);
        else % Includes if they are equal because the names don't matter in that case
            highNum = length(nucStats);
            lowNum  = length(haloStats);
        end
        
        pair  = zeros(lowNum, 1);
        dists = zeros(lowNum, 1);
        taken = zeros(highNum, 1);
        max   = 30; % max distance that matched nuclei are assumed to be away from one another

        % Find the closest centroid in the opposing list in order to ensure match
        for i = 1:lowNum
            min = max;
            for j = 1:highNum
                if ~taken(j)
                    if moreHalos
                        dist = sum(abs(nucStats(i).Centroid - haloStats(j).Centroid));
                    else
                        dist = sum(abs(haloStats(i).Centroid - nucStats(j).Centroid));
                    end
                    if dist < min
                        min = dist;
                        match = j;
                    end
                end
            end
            pair(i, 1)   = match;
            dists(i, 1)  = min;
            taken(match) = true;
        end
        pairs{n, 2} = pair;
        pairs{n, 3} = dists;
                    
        % If centroids are not within 30 pixels in x or y, print coordinates and warn user
        nearby = boolean(dists(:, 1) < max);
        warn = true;
        for i = 1:lowNum
            if ~nearby(i)
                if warn
                    warning('At least one nucleus-halo pair is misaligned in image %d', n)
                    fprintf('Length of haloStats = %d/n', length(haloStats))
                    fprintf('Length of nucStats  = %d/n', length(nucStats))
                end
                if moreHalos
                    fprintf('Xn: %.2f, Yn: %.2f\n', nucStats(i).Centroid)
                    fprintf('Xh: %.2f, Yh: %.2f\n', haloStats(pair(i, 1)).Centroid)
                else
                    fprintf('Xn: %.2f, Yn: %.2f\n', nucStats(pair(i, 1)).Centroid)
                    fprintf('Xh: %.2f, Yh: %.2f\n', haloStats(i).Centroid)
                end
            end
        end
        
         % Plot centroids with number next to it to compare. Colors are unspecific to regions, and
         % yellow corresponds with the adjusted value. 
        if debug
            figure(), imshow(imRegions(n).Image, []), title('Overlaid Centroids')
            hold on
            for i = 1:lowNum
                cn = nucStats(i).Centroid;
                ch = haloStats(i).Centroid;
                if moreHalos
                    cp = haloStats(pair(i, 1)).Centroid;
                    text(cn(1), cn(2), num2str(i), 'Color', 'r')
                else
                    cp = nucStats(pair(i, 1)).Centroid;
                    text(ch(1), ch(2), num2str(i), 'Color', 'r')
                end
                text(cp(1), cp(2), num2str(i), 'Color', 'y')
            end
            hold off
        end
    end
end


% Check if the given inputs are valid, returns boolean debug for use in the main function.
function debug = parseInputs(imRegions, args)
    
    % Default values
    debug = false;
    
    if nargin > 0
        if ~isempty(args)
            if nargin > 2
                error('Capstone:ArgChk', 'Function only accepts two input arguments.')
            end
            if islogical(args{1})
                debug = args{1};
            else
                error('Capstone:ArgChk', 'Second input argument must be of type boolean')
            end
        end
        if isstruct(imRegions)
            if ~any(strcmpi(fieldnames(imRegions), 'Halos'))
                error('Capstone:ArgChk', 'Struct input argument must contain field: ''Halos''')
            end
            if ~any(strcmpi(fieldnames(imRegions), 'Nuclei'))
                error('Capstone:ArgChk', 'Struct input argument must contain field: ''Nuclei''')
            end
            if debug
                if ~any(strcmpi(fieldnames(imRegions), 'Image'))
                    error('Capstone:ArgChk', 'Struct input argumetn must contain field: ''Image''')
                end
            end
        else
            error('Capstone:ArgChk', 'First input argument must be of type struct')
        end
    else
        error('Capstone:ArgChk', 'Not enough input arguments.')
    end
end
