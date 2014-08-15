classdef DmgCell < handle
    % Stores and manages data for cells being analyzed for DNA damage by the Fast Halo Assay. 
    % 
    %   Methods:
    %
    %       DmgCell     - Constructs and returns a handle to a DmgCell object.
    %       calcDamage  - Computes and returns the magnitude spectrum of the fourier transform.
    %       plot        - Plots the given image properties of the DmgCell in new figures.
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
    % Updated 7.1.14
    
    properties (SetAccess = private)
        % General
        image;      % Grayscale image of the cell region
        thresh;     % Thresholded image
        binary;     % Binary image
        area;       % Area of binary image
        radius;     % Equivalent radius of cell, given a circle with equal area
        pixels;     % List of pixel coordinates
        
        % Nucleus
        nucleus;    % Binary image of nucleus only
        nucArea;    % Area of nucleus
        nucRadius;  % Equivalent radius of nucleus, given a circle with equal area
        centroid;   % Center of the nucleus
        
        % Halo
        halo;       % Binary image of halo only
        haloPixels; % List of pixel coordinates in the halo
        
        % Damage Metrics
        ndf;        % Nuclear diffusion factor measurement
        hd;         % Halo distance measurement
        hm;         % Halo moment measurement
        ahm;        % Adjusted halo moment measurement
        ihi;        % Integrated halo intensity measurement
        ahi;        % Adjusted halo intensity measurement
        rhi;        % Relative halo intensity measurement
    end
    
    methods (Access = public)
    %% Constructor
        function self = DmgCell(image)
            % Constructs and returns a handle to a DmgCell object.
            % 
            %   dCell = ApopCell(image)
            %
            %       image   - any 8-bit image
            %       dCell   - A handle to the DmgCell object constructed (optional)

            
            % Check input
            validateattributes(image, {'uint8'}, {}, mfilename, 'image', 1);
            
            % Extract basic image info
            numThresh = 3;
            self.image = image;
            self.thresh = imquantize(self.image, multithresh(self.image, numThresh));
            self.binary = (self.thresh > 1);
            self.nucleus = (self.thresh > numThresh);
            self.halo = (self.binary - self.nucleus); 
            
            % Compute region statistics
            self.area = sum(self.binary(:));
            self.radius = sqrt(self.area / pi);
            self.nucArea = sum(self.nucleus(:));
            self.nucRadius = sqrt(self.nucArea / pi);
            [row, col] = find(self.nucleus);
            self.centroid = [mean(row), mean(col)];
            [row, col] = find(self.halo);
            self.haloPixels = [row, col];
        end
        
    %% public void calcDamage
        function varargout = calcDamage(self, varargin)
            % Calculates and returns the given DNA damage metrics for the cell.
            % 
            %   For a DmgCell object with handle 'dCell':
            %
            %   aCell.calcDamage(...)
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
            %   Scores are returned in the same order requested. Returns can be selectively ignored
            %   by replacing the output with a ~. all returns are optional to include. For example:
            %
            %       [ndf, ~, hm] = dCells.calcDamage('ndf', 'ahi', 'hm');
            %
            %   Will return the NDF and HM while calculating all three tests. 
            %
            %   All calculated values are stored in a field of the same name.   
            
            % Defaults
            getNDF = false;
            getHD = false;
            getHM = false;
            getAHM = false;
            getIHI = false;
            getAHI = false;
            getRHI = false;
            
            % Check inputs 
            if any(strcmpi(varargin, 'all'))
                varargin = {'ndf', 'hd', 'hm', 'ihi', 'ahi', 'rhi'};
            end
            for i = 1:numel(varargin)
                varargin{i} = lower(varargin{i});
                switch varargin{i}
                    case 'ndf'
                        getNDF = true;
                    case 'hd'
                        getHD = true;
                    case 'hm'
                        getHM = true;
                    case 'ahm'
                        getAHM = true;
                    case 'ihi'
                        getIHI = true;
                    case 'ahi'
                        getAHI = true;
                    case 'rhi'
                        getRHI = true;
                end
            end
            
            % Calculate nuclear diffusion factor
            if getNDF
                self.ndf = self.area / self.nucArea;
            end
           
            % Calculate halo distance
            if getHD
                self.hd = self.radius / self.nucRadius;
            end
            
            % Calculate halo moment, adjusted halo moment, and/or integrated halo intensity. These
            % are grouped together because they all rely on many of the same methods.
            if getHM || getAHM || getIHI || getAHI
                
                % All calculations rely on pixel intensity in halo
                intensities = self.image .* uint8(self.halo);
                
                % Moment calcs rely on the distance of each pixel to the centroid
                if getHM || getAHM
                    blank = zeros(size(self.image));
                    blank(round(self.centroid(1)), round(self.centroid(2))) = 1;
                    dist = bwdist(blank);
                    moments = double(intensities) .* dist;
                    if getHM
                        self.hm = sum(moments(:));
                    end
                    if getAHM
                        self.ahm = sum(moments(:)) / (self.nucRadius * sum(self.halo(:)));
                    end
                end
                
                % Intensity calcs
                if getIHI
                    self.ihi = sum(intensities(:));
                end
                if getAHI
                    self.ahi = sum(intensities(:)) / sum(self.halo(:));
                end
                if getRHI
                    self.rhi = sum(intensities(:)) / sum(self.image(:));
                end
            end
            
            % Return requested results
            varargout = cell(1, numel(varargin));
            for i = 1:numel(varargin)
                varargout{i} = self.(varargin{i});
            end
        end
        
    %% public void plot
        function [] = plot(self, varargin)
            % Plots the given image properties in new figures.
            %
            %   For a DmgCell object with handle 'dCell':
            %
            %   dCell.plot(...)
            %
            %       Valid input arguments (case insensitive, separate with commas):
            %
            %           'image'     - Plots the grayscale cell image
            %           'thresh'    - Plots the thresholded image
            %           'binary'    - Plots the binary image resulting from thresholding
            %           'halo'      - Plots the halo region as a binary image
            %           'nucleus'   - Plots the nucleus region as a binary image
            %           'centroid'  - Plots the centroid with lines extending to the projected edges
            %                         of the halo and the nucleus
            %           'damage'    - Plots the DNA damage metrics overlaid on the grayscale image
            %           'all'       - Plots all options
            
            % Check inputs
            modes = {'image', 'thresh', 'binary', 'halo', 'nucleus', 'centroid', 'damage'};
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
                    case 'thresh'
                        figure(), imshow(label2rgb(self.thresh)), title('Thresholded Image')
                    case 'binary'
                        figure(), imshow(self.binary), title('Binary Image')
                    case 'halo'
                        figure(), imshow(self.halo), title('Binary Halo')
                    case 'nuclei'
                        figure(), imshow(self.nucleus), title('Binary Nucleus')
                    case 'centroid'
                        figure(), imshow(self.image, []), hold on
                        c = self.centroid;
                        plot(c(2), c(1), 'ko')
                        plot([c(2), c(2)], [c(1), c(1) + self.radius], 'b')
                        plot([c(2), c(2) + self.nucRadius], [c(1), c(1)], 'r')
                        title('Centroid and halo/nuclear radii'), hold off
                    case 'damage'
                        figure(), imshow(self.image, []), hold on
                        shift = 0;
                        if ~isempty(self.ndf)
                            text(1, 8, sprintf('NDF: %d', self.ndf), 'Color', 'g')
                            shift = shift + 8;
                        end
                        if ~isempty(self.hd)
                            text(1, 8 + shift, sprintf('Halo Distance: %.5g', self.hd), 'Color', 'm')
                            shift = shift + 8;
                        end
                        if ~isempty(self.hm)
                            text(1, 8 + shift, sprintf('Halo Moment: %.5g', self.hm), 'Color', 'y')
                            shift = shift + 8;
                        end
                        if ~isempty(self.ahm)
                            text(1, 8 + shift, sprintf('Adj Halo Moment: %.5g', self.ahm), 'Color', 'c')
                            shift = shift + 8;
                        end
                        if ~isempty(self.ihi)
                            text(1, 8 + shift, sprintf('Halo Intensity : %d', self.ihi), 'Color', 'r')
                            shift = shift + 8;
                        end
                        if ~isempty(self.ahi)
                            text(1, 8 + shift, sprintf('Adj Halo Intensity : %d', self.ahi), 'Color', 'b')
                            shift = shift + 9;
                        end
                        if ~isempty(self.rhi)
                            text(1, 8 + shift, sprintf('Rel Halo Intensity : %d', self.rhi), 'Color', 'w')
                        end
                        title('Damage Metrics'), hold off
                end 
            end
        end
    end
end