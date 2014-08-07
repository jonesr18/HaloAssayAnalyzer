classdef HaloImage < handle
    % Manages and processes images of cells captured for assessment with the Fast Halo Assay. 
    % 
    %   Methods:
    %
    %       HaloImage   - Constructs and returns a handle to a HaloImage object.
    %       crop        - Crops the image.
    %       invert      - Inverts the image.
    %       setThresh   - Sets the values to be used for thresholding the image into 3 levels.
    %       setMinSize  - Sets the minimum object size to be evaluated in the image in terms of 
    %                     pixel area.
    %       process     - Processes the image to prepare it for segmentation. 
    %       setCriteria - Sets the criteria for exclusion of cells by eccentricity metric.
    %       makeCells   - Segments cells in the image & creates a cell array of appropriate objects.
    %       setCells    - Sets the fields 'aCells' and 'dCells' to appropriately passed arrays.
    %       calcDamage  - Calcuates the desired DNA damage parameters and statistics for each cell 
    %                     in the image.
    %       classify    - Classifies cells as apoptotic, healthy, and necrotic.
    %       plot        - Plots the given image properties in new figures.
    %       plotDamage  - Plots the desired DNA damage parameters on top of the grayscale image.
    %
    %   A handle to this object can be created as such:
    %
    %       image = rgb2gray(imread('filename'));
    %       dCell = DmgCell(image);
    %
    %
    % Copyright (C) 2014 Ross Jones, Singh Laboratory, University of Washington
    % 
    % Contact: jonesr18@gmail.com
    % 
    % This program is free software: you can redistribute it and/or modify it under the terms of the
    % GNU General Public License as published by the Free Software Foundation, either version 3 of 
    % the License, or (at your option) any later version.
    % 
    % This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
    % without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    % See the GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License along with this program. If
    % not see http://www.gnu.org/licenses/gpl.html
    % 
    % Updated 8.6.14
    
    properties (SetAccess = private)
        % Image properties
        image;      % Greyscale image
        thresh;     % Labels resulting from thresholding
        levels;     % Threshold levels
        binary;     % Binary image
        gradmag;    % Gradient magnitude of image
        
        % Cell properties
        labels;     % Labels resulting from segmentation
        boundaries; % Boundaries of cell objects
        eccen;      % Eccentricities of cell objects
        centroids;  % Centroids of each cell object
        radii;      % Radii of each cell object
        
        % User settings
        criteria;   % Eccentricity criteria
        minSize;    % Minimum number of pixels constituting a cell
        
        % Data
        aCells;     % Cell array of ApopCells constituting image
        dCells;     % Cell array of DmgCells constituting image
        damage;     % Struct array containing DNA damage metrics data
        cellTypes;  % Struct array containing cell classification data
    end
    
    properties (Access = private)
        processed;  % Becomes true once process() is called
        reset;      % Becomes true when setThresh() or setMinSize() are called
    end
    
    methods (Access = public)
    %% Constructor
        function self = HaloImage(image)
            % Constructs and returns a handle to a HaloImage object.
            % 
            %   hImage = HaloImage(image)
            %
            %       image   - any 8-bit image
            %       hImage  - A handle to the HaloImage object constructed (optional)
            
            % Check input
            validateattributes(image, {'uint8'}, {'real', 'nonsparse'}, mfilename, 'image');
            
            % Assign fields
            self.image = image;
            self.criteria = .500;       % default value
            self.minSize = 300;         % default value
            self.processed = false;
            self.reset = false;
            
            % Check if image requires inversion
            bwTest = im2bw(image, graythresh(image));
            if sum(bwTest(:)) > sum(~bwTest(:))
                self.invert();
            end
        end
        
        
    %% public void crop(int, int, int, int)
        function [] = crop(self, y1, y2, x1, x2)
            % Crops the image.
            %
            % The image being managed by the handle 'hImage' can be easily cropped as such:
            %
            %   hImage.crop(y1, y2, x1, x2)
            %       
            %       The arguments are the bounds of the cropped image, x is horozontal and y
            %       verticle. The 1 indicates the lower and 2 upper bound values. 
            %
            %   Note that this method will throw an error if called after processing.
            
            % Ensure image has not been processed
            if self.processed
                error('Method cannot be invoked after processing');
            end
            
            % Check inputs
            validateattributes(y1, {'numeric'}, {'positive'}, mfilename, 'y1', 1);
            validateattributes(y2, {'numeric'}, {'positive'}, mfilename, 'y2', 2);
            validateattributes(x1, {'numeric'}, {'positive'}, mfilename, 'x1', 3);
            validateattributes(x1, {'numeric'}, {'positive'}, mfilename, 'x2', 4);
            
            % Force inputs to be an integer
            y1 = round(y1);
            y2 = round(y2);
            x1 = round(x1);
            x2 = round(x2);
            
            % Crop image
            [H, W] = size(self.image);
            if y2 > H
                y2 = H;
            end
            if x2 > W
                x2 = W;
            end
            self.image = self.image(y1:y2, x1:x2);
        end
        
        
    %% public void invert()
        function [] = invert(self)
            % Inverts the image.
            %
            %   An image being managed by the handle 'hImage' can be easily inverted as such:
            %
            %   hImage.invert()
            %
            %   Note that this method will throw an error if called after processing. 
            
            % Ensure image has not been processed
            if self.processed
                error('Method cannot be invoked after processing');
            end
            
            % Invert image
            self.image = 255 - self.image;
        end
        
        
    %% public void setThresh(uint8, uint8)
        function [] = setThresh(self, lower, upper)
            % Sets the values to be used for thresholding the image into 3 levels.
            %
            %   Thresholding by default will be automatically completed using Otzu's method. To
            %   manually set the threshold levels, this method can be used.
            %
            %   For a HaloImage object with handle 'hImage':
            %
            %   hImage.setThresh(lower, upper)
            %
            %       Sets the thresholding values to the lower and upper values given. The values 
            %       will be automatically converted to an unsigned integer between 0 and 255.
            
            % Ensure image has not been processed
            if self.processed
                warning('Image has already been processed, reprocess to apply this change');
            end
            
            % Check input
            validateattributes(lower, {'numeric'}, {}, mfilename, 'lower', 1);
            validateattributes(upper, {'numeric'}, {}, mfilename, 'upper', 2);
            
            % Force levels to be an integer between 0 and 255
            lower = uint8(lower);
            upper = uint8(upper);
            
            % Assign levels
            self.levels = [lower, upper];
            self.reset = true;
        end
        
        
    %% public void setMinSize(int)
        function [] = setMinSize(self, pixels)
            % Sets the minimum object size to be evaluated in the image in terms of pixel area.
            %
            %   For a HaloImage object with handle 'hImage':
            %
            %   hImage.setMinSize(pixels)
            %
            %       Sets the minimum size to the value given by pixels.
            
            % Ensure image has not been processed
            if self.processed
                warning('Image has already been processed, reprocess to apply this change');
            end
            
            % Check input
            validateattributes(pixels, {'numeric'}, {'nonnegative'}, mfilename, 'pixels', 1);
            
            % Round to the nearest integer and assign field value
            self.minSize = round(pixels);
            self.reset = true;
        end
        
        
    %% public void process()
        function [] = process(self)
            % Processes the image to prepare it for segmentation. 
            %
            %   This method performs the actual image thesholding, addresss an unexplained error,
            %   and calculates the gradient magnitude of the image. 
            %
            %   For a HaloImage object with handle 'hImage':
            %
            %   hImage.process()
            
            % Warn if image has already been processed
            if self.processed && ~self.reset
                warning('Image has already been processed');
            end
            
            % Threshold
            if isempty(self.levels)
                self.levels = multithresh(self.image, 3);
            end
            TH = imquantize(self.image, self.levels);
                        
            % Reduce image borders should strange error occur
            [H, W] = size(TH);
            if sum(TH(1, :)) >= 2 * W
                TH(1, :) = TH(1, :) - 1;
            end
            if sum(TH(end, :)) >= 2 * W
                TH(end, :) = TH(end, :) - 1;
            end
            if sum(TH(:, 1)) >= 2 * H
                TH(:, 1) = TH(:, 1) - 1;
            end
            if sum(TH(:, end)) >= 2 * H
                TH(end, :) = TH(end, :) - 1;
            end
            self.thresh = TH;
            
            % Get gradient magnitude
            hy = fspecial('sobel');
            Iy = imfilter(double(self.image), hy, 'replicate');
            Ix = imfilter(double(self.image), hy', 'replicate');
            self.gradmag = sqrt(Ix.^2 + Iy.^2);
                                    
            % Update fields
            self.processed = true;
            self.reset = false;
        end      
        
        
    %% public void setCriteria(double)
        function [] = setCriteria(self, criteria)
            % Sets the criteria for exclusion of cells by eccentricity metric.
            %
            %   Cells with an eccentricity calculated greater than the metric will not be used for
            %   analysis if the parameter 'filterEccen' is passed to the method makeCells. This
            %   method can be used to set the criteria for this filtering. By default, the value is
            %   set to 0.500. 
            %
            %   For a HaloImage object with handle 'hImage':
            %
            %   hImage.setCriteria(criteria)
            %
            %       Sets the eccentricity passing criteria to the value given by criteria.
            
            % Check input
            validateattributes(criteria, {'numeric'}, {'>=', 0, '<=', 1}, mfilename, 'criteria', 1);
            
            % Uopdate field
            self.criteria = criteria;
        end
        
        
    %% public void makeCells(string, string, string, boolean)
        function [] = makeCells(self, type, mode, scheme, filterEccen)
            % Segments cells in the image & creates a cell array of appropriate objects.
            %
            %   The image must be processed before this method can be called, otherwise an error
            %   will be thrown. To do so, call the method process().
            %
            %   For a HaloImage object with handle 'hImage':
            %
            %   hImage.makeCells(type, mode, scheme, filterEccen)
            %   
            %       type            - The type of cells to make (damage/dmg, apoptosis/apop)
            %       mode            - The mode for extracting cells (binary/bw, fit/circles)
            %       scheme          - The scheme to make a binary (threshold/thresh, gradient/grad)
            %       filterEccen     - Boolean: Mark true if filtering cells by eccentricity
            %
            %   When making apoptosis cells, a hanning window is used to extract cells from the
            %   image to reduce FFT artifacts. When making damage cells, images are directly 
            %   extracted with no change in values to preseve them for calculations. 
            %
            %   Apoptosis cells are stored in the field 'aCells'.
            %   Damage cells are stored in the field 'dCells'.
            
            % Ensure image has been processed
            if ~self.processed
                error(['Image must be processed before cells can be made.\n',...
                       'Call the method process() before continuing.']);
            end
            
            % Check inputs
            type = lower(type);
            validatestring(type, {'damage', 'dmg', 'apoptosis', 'apop'}, mfilename, 'type', 1);
            mode = lower(mode);
            validatestring(mode, {'binary', 'bw', 'cirlces', 'fit'}, mfilename, 'mode', 2);
            scheme = lower(scheme);
            validatestring(scheme, {'threshold', 'thresh', 'gradient', 'grad'}, mfilename, 'scheme', 3);
            validateattributes(filterEccen, {'logical', 'numeric'}, {}, mfilename, 'filterEccen', 4);
            filterEccen = logical(filterEccen);
            
            % Get labels from segmentation scheme
            switch mode
                case {'binary', 'bw'}
                    [L, BW] = self.binarySegment(scheme, 1);
                case {'circles', 'fit'}
                    [L, BW] = self.circleSegment(scheme);
            end
            
            % Set binary image
            self.binary = BW;
            
            % If specified, only pass cell objects with eccentricity less than the criteria
            if filterEccen
                BW = (L > 0);
                [self.boundaries, Lbound] = bwboundaries(BW);
                stats = regionprops(Lbound, 'Eccentricity');
                for i = 1:length(stats)
                    self.eccen(i, 1) = stats(i).Eccentricity;
                    if stats(i).Eccentricity > self.criteria
                        BW = BW - (Lbound == i);
                    end
                end
                L = bwlabel(logical(BW), 4);
            end
            
            % Set labels
            self.labels = L;
            
            % Extract relevant regions from image and create desired cells
            switch type
                case {'damage', 'dmg'}
                    % Damage cells are created by simply extracting the region on the image that was
                    % labelled with the segmentation scheme. 
                    stats = regionprops(L, 'BoundingBox');
                    self.dCells = cell(1, length(stats));
                    for i = 1:length(stats);
                        % Find edges of cell on image
                        bb = stats(i).BoundingBox;
                        x1 = round(bb(1));
                        x2 = x1 + round(bb(3)) - 1;
                        y1 = round(bb(2));
                        y2 = y1 + round(bb(4)) - 1;

                        % Extract cell image and create damage cells
                        cellImage = self.image .* uint8(L == i);
                        self.dCells{i} = DmgCell(cellImage(y1:y2, x1:x2));
                    end
                
                case {'apoptosis', 'apop'}
                    % Apoptosis cells are created by passing a tukey window over the region on the 
                    % image that was labelled with the segmentation scheme. The reason for this is
                    % to smoothen the edges to reduce artifacting for the Fourier transform.
                    stats = regionprops(L, 'Area', 'Centroid');
                    self.aCells = cell(1, length(stats));
                    for i = 1:length(stats) 
                        % Estimate window length and force it to be even
                        len = round(4 * sqrt(stats(i).Area / pi));
                        if rem(len, 2) == 1
                            len = len + 1;
                        end
                        
                        % Create tukey hanning window
                        window = tukeywin(len, 0.5);
                        lin = linspace(-2, 2, len); % Goes from -2 : 2 in (len) steps
                        [x,y] = meshgrid(lin);
                        r = sqrt(x.^2 + y.^2);
                        
                        % Create 2D window
                        win2D = zeros(len);
                        win2D(:) = interp1(lin, window, r(:));
                        win2D(isnan(win2D)) = 0;
                        
                        % Apply window to image
                        cn = round(stats(i).Centroid);
                        x1 = cn(1) - len/2;
                        x2 = cn(1) + len/2 - 1;
                        y1 = cn(2) - len/2;
                        y2 = cn(2) + len/2 - 1;
                        try
                            cellImage = double(self.image(y1:y2, x1:x2)) .* win2D;
                        catch
                            [H, W] = size(self.image);
                            warning('Window out of bounds, shifting to fix cell %d', i)
                            
                            % Check x-bounds
                            if x1 < 1
                                x2 = x2 - x1 + 1;
                                x1 = 1;
                            elseif x2 > W
                                x1 = x1 - x2 + W;
                                x2 = W;
                            end
                            
                            % Check y-bounds
                            if y1 < 1
                                y2 = y2 - y1 + 1;
                                y1 = 1;
                            elseif y2 > H
                                y1 = y1 - y2 + H;
                                y2 = H;
                            end
                            
                            % Retry applying window
                            cellImage = double(self.image(y1:y2, x1:x2)) .* win2D;
                        end
                        
                        % Create apoptosis cells
                        self.aCells{i} = ApopCell(uint8(cellImage));
                    end
            end
        end
        
    %% public void setCells(cell)
        function [] = setCells(self, cellArray)
            % Sets the fields 'aCells' and 'dCells' to appropriately passed arrays.
            %
            %   The image must be processed before this method can be called, otherwise an error
            %   will be thrown. To do so, call the method process().
            %
            %   For a HaloImage object with handle 'hImage':
            %
            %   hImage.setCells(cellArray)
            %   
            %       cellArray can be a cell array of ApopCell or DmgCell objects. Depending on which
            %       type of object is contianed in the array, the method will replace the field
            %       'aCells' or 'dCells' with the array given. An error will be thrown if all
            %       objects in the cell array are not of the same class.
            %
            %   Apoptosis cells are stored in the field 'aCells'.
            %   Damage cells are stored in the field 'dCells'.
            %
            %   If aCells are given classified, the field 'cellTypes' will be updated accordingly.
            
            % Ensure image has been processed
            if ~self.processed
                error(['Image must be processed before cells can be made.\n',...
                       'Call the method process() before continuing.']);
            end
            
            % Check input
            validateattributes(cellArray, {'cell'}, {}, mfilename, 'cellArray', 1);
            if isa(cellArray{1}, 'DmgCell')
                arrayType = 'dCells';
            else
                arrayType = 'aCells';
            end
            
            % Check every cell in array is of consistent type
            if strcmp(arrayType, 'dCells')
                for i = 1:numel(cellArray)
                    if ~isa(cellArray{i}, 'DmgCell')
                        error('cellArray does not contain consistent cell object types');
                    end
                end
            else
                classifications = cell(1, numel(cellArray));
                for i = 1:numel(cellArray)
                    if ~isa(cellArray{i}, 'ApopCell')
                        error('cellArray does not contain consistent cell object types');
                    end
                    classifications{i} = cellArray{i}.cellType;
                end
                
                % Save number of each cell type if applicable
                for type = {'apoptotic', 'necrotic', 'healthy', 'ignored'}
                    self.cellTypes.(type{:}) = sum(strcmp(classifications, type{:}));
                end
            end
            
            % Set appropriate field
            self.(arrayType) = cellArray;
        end
        
    %% public void calcDamage(<string>)
        function [] = calcDamage(self, varargin)
            % Calcuates the desired DNA damage parameters and statistics for each cell in the image.
            %
            %   DmgCell objects must be created prior to using this method, otherwise an error will
            %   be thrown. To make them, call makeCells(...) with 'dmg' as the first parameter.
            %
            %   For a HaloImage object with handle 'hImage':
            %
            %   hImage.calcDamage(...)
            %
            %       Valid input arguments (case insensitive, separate with commas):
            %
            %           'ndf'   - Nuclear Diffusion Factor
            %           'hd'    - Halo Distance
            %           'hm'    - Halo Moment
            %           'ahm'   - Adjusted Halo Moment
            %           'ihi'   - Integrated Halo Intensity
            %           'ahi'   - Adjusted Halo Intensity
            %           'rhi'   - Relative Halo Intensity (DNA % Halo)
            %           'all'   - Calculate all metrics 
            %
            %   All calculated values for each cell are stored in the each respective cell object in
            %   the field 'dCells'. Means, standard deviations, and standard errors for all cells 
            %   are stored in the struct field 'damage' under the names 'Mean', 'Stdev', and 'SEM'.

            % Make sure damage cell objects have been created
            if isempty(self.dCells)
                error('DmgCell objects need to be initialized')
            end
            
            % Check Inputs
            tests = {'ndf', 'hm', 'ahm', 'hd', 'ihi', 'ahi', 'rhi'};
            if any(strcmpi(varargin, 'all'))
                varargin = tests;
            else
                for i = 1:numel(varargin)
                    varargin{i} = lower(varargin{i});
                    validatestring(varargin{i}, tests, mfilename, 'metric', i);
                end
            end

            % Calculate damage
            for c = self.dCells
                c{:}.calcDamage(varargin{:});
            end

            % Calculate statistics
            for test = varargin
                data = zeros(size(self.dCells));
                for i = 1:numel(self.dCells)
                    data(i) = self.dCells{i}.(test{:});
                end
                self.damage.(test{:}) = struct('Mean', mean(data),...
                                               'Stdev', std(data),...
                                               'SEM', std(data) / sqrt(numel(self.dCells)));
            end
        end
        
    %% public void classify(uint8, string, <string>)
        function [] = classify(self, data, answers, varargin)
            % Classifies cells as apoptotic, healthy, and necrotic.
            %
            %   ApopCell objects must be created prior to using this method, otherwise an error will
            %   be thrown. To make them, call makeCells(...) with 'apop' as the first parameter.
            %
            %   For a HaloImage object with handle 'hImage':
            %
            %   hImage.classify(data, answers)
            %
            %       Builds a classifier (k-nearest neighbors by default) using the given data and
            %       answers. 'data' is a N x M matrix contianing M degrees of magnitude spectrum
            %       data from N scored cells. 'answers' is an N x 1 cell array of answers which can
            %       only contain the strings 'apoptotic', 'necrotic', 'healthy', and 'ignored'.
            %
            %       The classifier is used to label each ApopCell in the field 'aCells' as either
            %       'apoptotic', 'necrotic', 'healthy' or 'ignored'.
            %
            %   hImage.classify(data, answers, OPTIONS)
            %
            %       The following strings can be passed as options:
            %
            %           Select Classifier:
            %           'knn'       - k-nearest neighbors algorithm                     (default)
            %           'kmeans'    - k-means algoritm*
            %           'ctree'     - classification (decision) tree
            %
            %           Select Data Types:
            %           'magspec'   - use select frequencies of magnitude spectrum
            %                         * Note that if this option is selected, then a Mx1 matrix of
            %                           frequency coordinates must be passed as well. Coordinates
            %                           are indexes of the 1D reshaped frequency spectrum.
            %           'allbins'   - use binned magnitude spectrum
            %           'innerbins' - use binned inner magnitude spectrum               (defualt)
            %
            %           * will not give 'ignored' as an answer
            %
            %   Cell types in ApopCell objects can be found in their field, 'cellType'.
            %   
            %   The number of each cell type in the image is stored in the field 'cellTypes', which
            %   is a struct array with fieldnames 'apoptotic', 'necrotic', 'healthy', and 'ignored'. 
            
            % Make sure apoptosis cell objects have been created
            if isempty(self.aCells)
                error('ApopCell objects need to be initialized')
            end
            
            % Check inputs
            narginchk(3, 6);
            [classifier, dataType, freqCoords] = checkClassifyInputs(varargin{:});
            %   Note that freqCoords is returned empty if not using 'magspec' dataType
            
            % Find indexes of each cell type in data
            a = strcmp(answers, 'apoptotic');
            n = strcmp(answers, 'necrotic');
            h = strcmp(answers, 'healthy');
            i = strcmp(answers, 'ignored');
            
            % Make desired classifier and classify cells accordingly
            predictions = cell(1, numel(self.aCells));
            switch classifier
                case 'knn'
                    % Classifies cells based on known data that is most similar by distance.
                    
                    % Build classifier
                    knn = ClassificationKNN.fit(data, answers, 'DistanceWeight', 'squaredinverse');
                    
                    % Classify each cell
                    for j = 1:length(self.aCells)
                        [fourier, allbins, innerbins] = self.aCells{j}.fft();
                        cellData = makeCellData(dataType, freqCoords);
                        prediction = knn.predict(cellData);
                        self.aCells{j}.setType(prediction{:});
                        predictions(j) = prediction;
                    end
                case 'kmeans'
                    % Makes a classifier using supervised K-means clustering. Cells are classified
                    % based on which Voronoi fragment they belong to. 
                    
                    % Find centroids of each cluster
                    [~, cn_a] = kmeans(data(a, :), 3, 'replicates', 3);
                    [~, cn_n] = kmeans(data(n, :), 3, 'replicates', 3);
                    [~, cn_h] = kmeans(data(h, :), 3, 'replicates', 3);
                    [~, cn_i] = kmeans(data(i, :), 3, 'replicated', 3);
                    
                    % Create lookup array for quick prediction reference
                    lookup = cell(1, 12);
                    lookup(1:3) = {'apoptotic'};    lookup(4:6) = {'necrotic'}; 
                    lookup(7:9) = {'healthy'};      lookup(10:end) = {'ignored'};
                    
                    % Classify each cell
                    for j = 1:numel(self.aCells)
                        [fourier, allbins, innerbins] = self.aCells{j}.fft();
                        cellData = makeCellData(dataType, freqCoords);
                        dists = dist([cellData, cn_a', cn_n', cn_h', cn_i']);
                        dists = dists(2:end, 1:end - 1); % TODO fix
                        [~, idx] = min(dists);
                        prediction = lookup{idx};
                        self.aCells{j}.setType(prediction);
                        predictions{j} = prediction;
                    end                    
                case 'ctree'
                    % Makes a decision tree to classify cell based on known data.
                    
                    % Build classifier
                    tree = ClassificationTree.fit(data, answers');
                    
                    % Classify each cell
                    for j = 1:numel(self.aCells)
                        [fourier, allbins, innerbins] = self.aCells{j}.fft();
                        cellData = makeCellData(dataType, freqCoords);
                        prediction = tree.predict(cellData);
                        self.aCells{j}.setType(prediction{:});
                        predictions(j) = prediction;                    
                    end
            end
            
            % Save number of each cell type
            for type = {'apoptotic', 'necrotic', 'healthy', 'ignored'}
                self.cellTypes.(type{:}) = sum(strcmp(predictions, type{:}));
            end
            
            
            % Checks the inputs for this function
            function [classifier, dataType, freqCoords] = checkClassifyInputs(varargin)
                
                % Validate data and answers
                validateattributes(data, {'numeric'}, {'positive'}, mfilename, 'allMags', 1);
                validateattributes(answers, {'cell'}, {}, mfilename, 'answers', 2);
                [nObs, nDeg] = size(data);
                if nObs ~= numel(answers)
                    error(['Mismatched number of observations in data and answers',...
                        'Data: %d\nAnswers: %d\n'], nObs, numel(answers));
                end
                
                % Check for optional inputs
                if isempty(varargin)
                    classifier = 'knn';         % Default is KNN classifier
                    dataType = 'innerbins';     % Default is to use inner bins
                else
                    for k = 1:numel(varargin)
                        if ischar(varargin{k})
                            arg = lower(varargin{k});
                            if any(strcmp(arg, {'knn', 'kmeans', 'ctree'}));
                                % Classifier option
                                classifier = arg;
                            elseif any(strcmp(arg, {'magspec', 'allbins', 'innerbins'}))
                                % Data type option
                                dataType = arg;
                            else
                                % Incorrect option
                                error('Argument %d did not match any of the following strings: \n%s',...
                                    k + 3, 'knn, kmeans, ctree, magspec, allbins, innerbins');
                            end
                        else
                            % The only non-string inputs should be frequency coordinates, which must
                            % be passed if the user selects the 'magspec' data type option
                            freqCoords = uint64(varargin{k});
                            validateattributes(freqCoords, {'numeric'}, {'ncols', 2, 'nrows', nDeg},...
                                mfilename, 'freqCoords', k);
                        end
                    end
                end
                
                % Check if frequency coordinates were passed if user is using 'magspec' data type
                if strcmp(dataType, 'magspec')
                    if ~exist('freqCoords', 'var')
                        error(['If using ''magspec'' data type, you must pass freqency coordinates ',...
                            'to the function']);
                    end
                else
                    freqCoords = [];
                end
            end
            
            
            % Makes the cell data for classification purposes
            function cellData = makeCellData(dataType, freqCoords)
                switch dataType
                    case 'magspec'
                        % Find which frequencies to use
                        cellData = zeros(1, size(data, 2));
                        for k = 1:length(cellData)
                            f = reshape(fourier, 1, []);
                            cellData(k) = f(freqCoords(k));
                        end
                    case 'allbins'
                        cellData = reshape(allbins, 1, []);
                    case 'innerbins'
                        cellData = reshape(innerbins, 1, []);
                end
            end
        end
        
    %% public void plot(<string>)
        function [] = plot(self, varargin)
            % Plots the given image properties in new figures.
            %
            %   For a HaloImage object with handle 'hImage':
            %
            %   hImage.plot(...)
            %
            %       Valid input arguments (case insensitive, separate with commas):
            %
            %           'image'     - Plots the grayscale cell image
            %           'contour'   - Plots the contour image
            %           'thresh'    - Plots the thresholded image
            %           'binary'    - Plots the binary image resulting from thresholding
            %           'eccen'     - Plots lines outlining each object and prints the eccentricity 
            %                         metric calculated next to objects
            %           'gradmag'   - Plots the gradient magnitude of the image
            %           'circles'   - Plots circles fitted to the gradient magnitude or nuclei on 
            %                         top of both the gradient magnitude and grayscale image
            %           'labels'    - Plots labels for cells passing criteria and used for analysis
            %           'all'       - Plots all options

            
            % Check inputs
            modes = {'image', 'contour', 'thresh', 'binary', 'eccen', 'gradmag', 'circles', 'labels'};
            if any(strcmpi(varargin, 'all'))
                varargin = modes;
            else
                for i = 1:numel(varargin)
                    varargin{i} = lower(varargin{i});
                    validatestring(varargin{i}, modes, mfilename, 'mode', i);
                end
            end
            
            % Plot requested images
            for mode = varargin
                switch mode{:}
                    case 'image'
                        figure(), imshow(self.image, []), title('Grayscale Image')
                    case 'contour'
                        figure(), contour(self.image), title('Contour Map of Image')
                    case 'thresh'
                        if ~self.processed
                            warning('Image thresholds not yet calculated; call process()')
                        else
                            figure(), imshow(label2rgb(self.thresh)), title('Thresholded Image')
                        end
                    case 'binary'
                        if isempty(self.binary)
                            warning('Binary image not yet calculated; call makeCells()')
                        else
                            figure(), imshow(self.binary), title('Binary Image')
                        end
                    case 'eccen'
                        if isempty(self.eccen)
                            warning('Cell object eccentricities not yet calculated')
                        else
                            [H, W] = size(self.image);
                            figure(), imshow(self.image, []), hold on
                            for i = 1:numel(self.boundaries)
                                boundary = self.boundaries{i};
                                eccenStr = sprintf('%.3f', self.eccen(i));
                                plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2)
                                text(boundary(1, 2) - W / 30, boundary(1, 1) + H / 30, eccenStr, 'Color', 'g');
                                if self.eccen(i) < self.criteria
                                    plot(mean(boundary(:, 2)), mean(boundary(:, 1)), 'ko');
                                end
                            end
                            title('Cell Objects Passing Eccentricity Criteria Marked ''o'''), hold off
                        end
                    case 'gradmag'
                        if isempty(self.gradmag)
                            warning('Gradient magnitude not yet calculated; call process()')
                        else
                            gm = self.gradmag;
                            figure(), imshow(gm, [min(gm(:)), max(gm(:)/4)]), title('Gradient Magnitude')
                        end
                    case 'circles'
                        if isempty(self.centroids)
                            warning('Circles not yet fitted')
                        else
                            gm = self.gradmag;
                            figure(), imshow(gm, [min(gm(:)), max(gm(:)/4)])
                            viscircles(self.centroids, self.radii);
                            title('Fitted Circles on Gradient Magnitude')
                            figure(), imshow(self.image, []), viscircles(self.centroids, self.radii);
                            title('Fitted Circles on Original Image')
                        end
                    case 'labels'
                        if isempty(self.labels)
                            warning('Labels not yet created; call makeCells()')
                        else
                            figure(), imshow(label2rgb(self.labels)), title('Cell Labels')
                        end
                end
            end
        end
        
    %% public void plotDamage(<string>)
        function [] = plotDamage(self, varargin)
            % Plots DNA damage metric scores on top of the grayscale image.
            %
            %   DNA damage metrics must be calculated prior to using this method, otherwise an error
            %   will be thrown. To do so, call calcDamage(...). 
            %
            %   For a HaloImage object with handle 'hImage':
            %
            %   hImage.plotDamage(...)
            %
            %       Valid input arguments (case insensitive, separate with commas):
            %
            %           'ndf'   - Nuclear Diffusion Factor
            %           'hd'    - Halo Distance
            %           'hm'    - Halo Moment
            %           'ahm'   - Adjusted Halo Moment
            %           'ihi'   - Integrated Halo Intensity
            %           'ahi'   - Adjusted Halo Intensity
            %           'rhi'   - Relative Halo Intensity (DNA % Halo)
            %           'all'   - Plot all metrics 
            %
            %       Metrics not computed but passed as arguments will be ignored.
            %
            %   Metrics are printed color-coded with their name, which is in the top left corner.
            
            % Plot desired metrics
            if isempty(self.damage)
                warning('Damage not calculated')
                return
            end

            % Check inputs
            modes = {'ndf', 'hd', 'hm', 'ahm', 'ihi', 'ahi', 'rhi'};
            if any(strcmpi(varargin, 'all'))
                varargin = modes;
            else
                for i = 1:numel(varargin)
                    varargin{i} = lower(varargin{i});
                    validatestring(varargin{i}, modes, mfilename, 'mode', i);
                end
            end

            % Adjust for image size
            [H, ~] = size(self.image);
            spacer = round(H / 50);

            % Plot data on each cell
            stats = regionprops(self.labels, 'Centroid');
            colors = 'rygcbmw';
            figure(), imshow(self.image, []), hold on
            for i = 1:length(stats)
                c = self.dCells{i};
                cn = stats(i).Centroid;
                shift = 0;
                for j = 1:numel(varargin)
                    if ~isempty(c.(varargin{j}))
                        text(cn(1), cn(2) + shift, sprintf('%.4g',c.(varargin{j})), 'Color', colors(j))
                        shift = shift + spacer;
                    end
                end
            end

            % Plot names of metrics in top left corner, color-coded
            shift = 0;
            for i = 1:numel(varargin)
                switch varargin{i}
                    case 'ndf'
                        text(1, 10 + shift, 'Nuclear Diffusion Factor', 'Color', colors(i))
                    case 'hd'
                        text(1, 10 + shift, 'Halo Distance', 'Color', colors(i))
                    case 'hm'
                        text(1, 10 + shift, 'Halo Moment', 'Color', colors(i))
                    case 'ahm'
                        text(1, 10 + shift, 'Adjusted Halo Moment', 'Color', colors(i))
                    case 'ihi'
                        text(1, 10 + shift, 'Integrated Halo Intensity', 'Color', colors(i))
                    case 'ahi'
                        text(1, 10 + shift, 'Adjusted Halo Intensity', 'Color', colors(i))
                    case 'rhi'
                        text(1, 10 + shift, 'Relative Halo Intensity', 'Color', colors(i))
                end
                shift = shift + spacer;
            end
            title('Damage Metrics for Each Measured Cell'), hold off
        end
        
    %% public void plotTypes()
        function [] = plotTypes(self)
            % Plots the type of scored cell on top of the grayscale image.
            % 
            %   Cell classification must be completed prior to using this method, otherwise an error
            %   will be thrown. To do so, call classify(...). 
            %
            %   For a HaloImage object with handle 'hImage':
            %
            %   hImage.plotTypes(...)
            
            if isempty(self.cellTypes)
                warning('Cells not classified')
                return
            end
            
            % Plot label on each cell
            stats = regionprops(self.labels, 'Centroid');
            figure(), imshow(self.image, []), hold on
            for i = 1:length(stats)
                c = self.aCells{i};
                cn = stats(i).Centroid;
                text(cn(1), cn(2), c.cellType, 'Color', 'r')
            end
        end
    end
    
    methods (Access = private)
    %% private <uint8, boolean> binarySegment(string, int)
        function [varargout] = binarySegment(self, scheme, level)
            % Segments the image using a binary image extraction technique.
            %
            % Optional returns are the resulting labels and binary image.
            
            switch scheme
                case {'threshold', 'thresh'}
                    % Segment image by extracting regions with greatest pixel intensity
                    BW = (self.thresh > level);
                case {'gradient', 'grad'}
                    % Segment image by edge detection of cell regions
                    gradThresh = imquantize(self.gradmag, multithresh(self.gradmag, 2));
                    BW = (gradThresh > level);
                    
                    % Dilate the gradient -- closing and opening smoothes edges
                    strel1   = strel('disk', 1);
                    BW = imopen(BW, strel1);
                    BW = imdilate(BW, strel('disk', 3));
                    BW = imdilate(BW, strel('diamond', 3));
                    BW = imopen(BW, strel1);
                    BW = imclose(BW, strel1);
                    BW = imopen(BW, strel1);
                    BW = imclose(BW, strel1); 
            end
            
            % Remove regions with area < minSize
            BW = imclearborder(BW, 4);
            BW = imfill(BW, 'holes');
            BW = ~imdilate(~BW, strel('disk', 2)); % Remove very tiny specs to improve speed
            L = bwlabel(BW, 4);
            for i = 1:max(L(:))
                if sum(L(L == i)) / i < self.minSize
                    BW = BW - (L == i);
                end
            end
            
            % Return items
            BW = logical(BW);
            varargout{1} = bwlabel(BW, 4);
            varargout{2} = BW;
        end
        
    %% private <uint8, boolean> circleSegment(string)
        function [varargout] = circleSegment(self, scheme)
            % Segment image by fitting circles to an aspect of the image that is determined by the
            % given scheme.
            %
            % Optional returns are the resulting labels and binary image.
            
            % Get binary images from threshSegment
            [~, BWfit] = self.binarySegment(scheme, 2);
            [~, BW] = self.binarySegment(scheme, 1);
            
            % Collect stats on newly processed gradient
            statsNew = regionprops(BWfit, 'EquivDiameter');
            radList = [statsNew.EquivDiameter] / 2;
            rad = round(median(radList));
            radstd = round(std(radList));

            % Fit circles to the image to find individual nuclei
            % NOTE: Circle finding is difficult with low-resolution pictures. It also had
            %       issues when CMYK was converted to RGB due to poor contrast.
            range = [rad - radstd, rad + radstd];
            [cn, r] = imfindcircles(BWfit, range, 'ObjectPolarity', 'bright', 'Sensitivity', 0.92);
            if isempty(cn)
                % If no circles found, just return labelled BW image
                warning('No circles were able to be fit')
                varargout{1} = bwlabel(BW, 4);
                varargout{2} = BW;
            else % If circles are found, use them
                % Combine overfit circles
                L = self.tidyFits(cn, r, BW);

                % Return items
                varargout{1} = L;
                varargout{2} = BW;
            end
        end
        
    %% private uint8 tidyFits(double, double)    
        function [labels] = tidyFits(self, cn, r, BW)
            % Makes optimizations for circles fitted. 
            
            % Combine centroids within 20 pixels (overfitted circles)
            allDistances = dist(cn');
            [x, y] = find(allDistances < 20);
            indexes = find(x ~= y);
            for i = indexes'
                average = mean([cn(x(i), :); cn(y(i), :)]);
                cn(x(i), :) = average;
                cn(y(i), :) = average;
                radAvg = mean([r(x(i)), r(y(i))]);
                r(x(i)) = radAvg;
                r(y(i)) = radAvg;
            end
            self.centroids = cn;
            self.radii = r;

            % Create blank binary image with 1s at centroid locations
            blank = zeros(size(BW));
            for c = cn'
                x = round(c(1));
                y = round(c(2));
                blank(y, x) = 1;
            end
            
            % Only pass cell objects with one centroid located within them.
            L = bwlabel(BW, 4);
            for i = 1:max(L(:))
                idx = L == i;
                initArea = sum(L(idx)) / i;
                L = L + blank;
                newArea = sum(L(idx)) / i;
                diff = abs(newArea - initArea - 1 / i);
                if diff > .0001
                    BW = BW - logical(idx);
                end
                L = L - blank;
            end
            labels = bwlabel(logical(BW), 4);
        end
    end
end