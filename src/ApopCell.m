classdef ApopCell < handle
    % Stores and manages data for cells being analyzed for apoptosis / cell classification by the
    % Fast Halo Assay. 
    % 
    %   Methods:
    %
    %       ApopCell    - Constructs and returns a handle to an ApopCell object.
    %       fft         - Computes and returns the magnitude spectrum of the fourier transform.
    %       setType     - Sets the cell type.
    %       plot        - Plots the given image properties of the ApopCell in new figures.
    %
    %   A handle to this object can be created as such:
    %
    %       image = rgb2gray(imread('filename'));
    %       aCell = ApopCell(image);
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
        image;      % Grayscale image
        fourier;    % 2D FFT magnitude spectrum
        magbins;    % Magnitude spectrum averaged into bins of 10 each
        cellType;   % The type of cell represented by the region (control, apoptotic, necrotic)
    end
    
    methods (Access = public)
    %% Constructor
        function self = ApopCell(image)
            % Constructs and returns a handle to an ApopCell object.
            % 
            %   aCell = ApopCell(image)
            %
            %       image   - any 8-bit image
            %       aCell   - A handle to the ApopCell object constructed (optional)
            
            % Check input
            validateattributes(image, {'uint8'}, {}, mfilename, 'image', 1);
            
            % Assign image field
            self.image  = image;
        end
        
    %% public <double, double> fft()
        function varargout = fft(self)
            % Computes and returns the magnitude spectrum of the fourier transform.
            %
            %   The magnitude spectrum is normalized to a size of 200x200 pixels to facilitate 
            %   comparisons between images. 
            %
            %   For an ApopCell object with handle 'aCell':
            %
            %   [fourier, magbins] = aCell.fft()
            %
            %       Returns the reduced bins of the fourier transform magnitude spectrumm. The first
            %       returned variable can be selectively ignored by replacing it with a ~. Both 
            %       returns are optional to include.
            %
            %   The magnitude spectrum is stored in the field 'fourier'.
            %   The binned magnitudes are stored in the field 'magbins'.
            
            if isempty(self.fourier)
                self.fourier = abs(fftshift(fft2(double(self.image), 200, 200)));
                self.magbins = zeros(10, 10);
                for i = 1:10
                    for j = 1:10
                        rows = 20 * (i - 1) + 1 : 20 * i;
                        cols = 20 * (j - 1) + 1 : 20 * j;
                        self.magbins(i, j) = mean(mean(self.fourier(rows, cols)));
                    end
                end
            end
            varargout = {self.fourier, self.magbins};
        end
        
    %% public void setType(string)
        function [] = setType(self, type)
            % Sets the cell type.
            %
            %   For an ApopCell object with handle 'aCell':
            %
            %   aCell.setType(type)
            %
            %       Type can be any of the following (case insensitive):
            %
            %           'apoptotic'
            %           'necrotic'
            %           'healthy'
            %           'ignored'
            
            % Check input
            type = lower(type);
            validatestring(type, {'apoptotic', 'necrotic', 'healthy', 'ignored'}, mfilename, 'type', 1);
            
            % Set cell type
            self.cellType = type;
        end
        
    %% public void plot(string)
        function [] = plot(self, varargin)
            % Plots the given image properties in new figures.
            %
            %   For an ApopCell object with handle 'aCell':
            %
            %   aCell.plot(...)
            %
            %       Valid input arguments (case insensitive, separate with commas):
            %
            %           'image'     - Plots the grayscale cell image
            %           'fourier'   - Plots the Fourier transform magnitude spectrum
            %           'magbins'   - Plots the binnned magnitude spectrum
            %           'all'       - Plots all options
            
            % Check inputs
            modes = {'image', 'fourier', 'magbins'};
            if any(strcmpi(varargin, 'all'))
                varargin = modes;
            else
                for i = 1:numel(varargin)
                    varargin{i} = lower(varargin{i});
                    validatestring(varargin{i}, modes, mfilename, 'mode', i);
                end
            end
            
            warned = false;
            % Plot requested images
            for mode = varargin
                switch mode{:}
                    case 'image'
                        figure(), imshow(self.image, []), title('Grayscale Image')
                    case 'fourier'
                        if isempty(self.fourier) 
                            if ~warned
                                warning('Fourier transform not yet computed')
                                warned = true;
                            end
                        else
                            figure(), imshow(self.fourier, [0, 10000]), title('Magnitude Spectrum')
                        end
                    case 'magbins'
                        if isempty(self.magbins)
                            if ~warned
                                warning('Fourier transform not yet computed')
                                warned = true;
                            end
                        else
                            figure(), imshow(self.magbins, [0, 1000]), title('Magnitude Bins for Sorting')
                        end
                end
            end
        end
    end
end