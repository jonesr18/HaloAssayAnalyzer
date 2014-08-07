classdef Exports < handle
    % A functional class for saving data in standard CSV format from the FHA.
    %
    %   Methods (static):
    %
    %       haloData     - Exports data to a standard fast halo assay data file.
    %       learningData - Writes learning data to standard data files.
    %       saveConfig   - Saves the configuration data of the Initialization GUI.
    %
    %   A handle to this object can be created as such:
    %
    %       out = Exports();
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
    
    methods (Static)
    %% public static void haloData(string, cell)
        function [] = haloData(fullfile, data)
            % Writes data to a standard fast halo assay data file. This method is static and thus
            % can be invoked without first creating a handle by calling Exports.haloData(...).
            %
            %   For an Exports object with handle 'out':
            %
            %   out.haloData(fullfile, data)
            %
            %       fullfile    - The full filename to be written ('.csv' auto-added if omitted)
            %       data        - A cell array containing all data and labels in correct form
            %
            %   The correct format for the parameter data is a cell array as shown below. Incorrect
            %   formatting will cause the output have random numbers (character values) and numbers
            %   which are formatted as strings. The top level is primarily data labels, and the 
            %   cells containing doubles hold data and statistics. The left two columns designate 
            %   cell number and statistic names.
            %
            %       [string]    [string]    [string]    [string]    ....    [string]
            %       [string]    [string]    [double]    [double]    ....    [double]
            %         ....        ....        ....        ....      ....      ....
            %       [string]    [string]    [double]    [double]    ....    [double]
            
            % Check inputs
            validateattributes(fullfile, {'char'}, {}, mfilename, 'fullfile', 1);
            validateattributes(data, {'cell'}, {'2d'}, mfilename, 'data', 2);
            if ~strcmp(fullfile(end - 3 : end), '.csv')
                fullfile = strcat(fullfile, '.csv');
            end
            
            % Write file
            fileID = fopen(fullfile, 'wt');
            [H, W] = size(data);
            fprintf(fileID, '%s,', data{1, 1:W - 1});            % First row is always all strings
            fprintf(fileID, '%s\n', data{1, W});                 % Last value goes to next line
            for i = 2:H
                fprintf(fileID, '%s,', data{i, 1:2});            % First two columns are always strings
                if strcmp(data{i, 3}, 'N/A')
                    fprintf(fileID, '%s,', data{i, 3:W - 1});    % Rest of columns are all strings (N/A)
                    fprintf(fileID, '%s\n', data{i, W});         % Last value goes to next line
                else
                    fprintf(fileID, '%f,', data{i, 3:W - 1});    % Rest of the columns are all numbers
                    fprintf(fileID, '%f\n', data{i, W});         % Last value goes to next line
                end
            end
            fclose(fileID);
        end
        
    %% public static void learningData(int, double, double, double, cell) 
        function [] = learningData(magSpec, allbins, innerbins, answers)
            % Writes learning data to standard data files. This method is static and thus can be
            % invoked without first creating a handle by calling Exports.learningData(...).
            %
            %   For an Exports object with handle 'out':
            %
            %   out.learningData(freqCoords, magSpec, allbins, innerbins, answers)
            %
            %       magSpec     - The complete 2D FFT magnitude spectrum of each cell
            %       allbins     - The complete magnitude spectrum averaged into 100 20x20 bins
            %       innerbins   - Low-frequency magnitude spectrum averaged into 100 4x4 bins
            %       answers     - Manual cell classification
            %
            %   Each of magSpec, allbins, and innerbins is and N x M array, where N is the number of
            %   observations, corresponding with the length of answers. M is variable for magSpec
            %   and 100 for both allbins and innerbins, corresponding with the reshaped (and 
            %   potentially averaged) parts of the magnitude spectrum.

            
            % Check inputs
            validateattributes(magSpec, {'numeric'}, {'>=' 0}, mfilename, 'magSpec', 2);
            validateattributes(allbins, {'numeric'}, {'ncols', 100, '>=' 0}, mfilename, 'allbins', 3);
            validateattributes(innerbins, {'numeric'}, {'ncols', 100, '>=', 0}, mfilename, 'innerbins', 4);
            validateattributes(answers, {'cell'}, {}, mfilename, 'answers', 5);
            freqCoords = findFreqCoords(magSpec, answers);
            selMags = magSpec(:, freqCoords);
            
            % Validate size of input
            nAns = numel(answers);
            errmsg = 'The number of observations and answers must be the same.\nnObs (%s): %d\nnAns: %d';
            namepairs = {'Selected Freqs', 'All Freq Bins', 'Inner Freq Bins'
                          selMags,          allbins,         innerbins};
            for i = 1:size(namepairs, 2)
                if nAns ~= size(namepairs{2, i}, 1)
                    error(errmsg, namepairs{1, i}, size(namepairs{2, i}, 1), nAns);
                end
            end
            
            % Create directory
            f = datestr(now);
            f(f == ':') = ';';
            f(f == '-') = ' ';
            fParent = 'Learning Data';
            if ~exist(strcat(pwd, '\', fParent), 'dir')
                mkdir(fParent)
            end
            mkdir(fParent, f);
            filepath = strcat(pwd, '\', fParent, '\', f, '\');
            
            % Write magnitude specturm contingent files
            for filename = namepairs;
                fileID = fopen(strcat(filepath, filename{1}, '.csv'), 'wt');
                fprintf(fileID, '%s %s\n', f, filename{1});
                for i = 1:size(filename{2}, 1)
                    fprintf(fileID, '%d,', filename{2}(i, 1:end-1)); % All values are numbers
                    try
                        fprintf(fileID, '%d\n', filename{2}(i, end));    % Last value goes to next line
                    catch
                        if isempty(freqCoords)
                            warning('No frequencies selected, ''Frequency Coords.csv'' will be empty')
                            break
                        end
                    end
                end
                fclose(fileID);
            end
            
            % Write frequency coordinates file
            fileID = fopen(strcat(filepath, 'Frequency Coords.csv'), 'wt');
            fprintf(fileID, '%s %s\n', f, 'Frequency Coords');
            for i = 1:length(freqCoords)
                fprintf(fileID, '%d\n', freqCoords(i));
            end
            fclose(fileID);
            
            % Write answers file
            fileID = fopen(strcat(filepath, 'Answers.csv'), 'wt');
            fprintf(fileID, '%s %s\n', f, 'Answers');
            for i = 1:numel(answers)
                fprintf(fileID, '%s\n', answers{i});
            end
            fclose(fileID);
        end
        
    %% public static void saveConfig(struct)
        function [] = saveConfig(handles)
            % Saves the configuration data of the Initialization GUI.
            %
            %   For an Exports object with handle 'out':
            %
            %   out.saveConfig(handles)
            %
            %       handles     - The handles data struct from the GUI
            %
            %   The data is saved as a struct named configData in the file config.mat. To access
            %   this data later, simply do the following:
            %
            %       load config.mat
            %
            %   and the variable configData will be loaded into the workspace appropriately.
            
            % Handles to ignore
            ignored = {'figure1', 'licence', 'about', 'text10', 'segmentationPanel',...
                'uipanel11', 'uipanel10', 'uipanel9', 'uipanel3', 'run', 'uipanel6', 'uipanel1',...
                'text9', 'text8', 'text7', 'outputBrowse', 'text5', 'text4', 'text6', 'axes',...
                'images', 'currImage', 'output', 'imageMenu'};
            
            % Handles which require accessing their string
            strings = {'lowerThresh', 'upperThresh', 'eccenMetric', 'outputPath', 'sampleName'};
            
            % Handles which have direct data
            direct = {'outdir', 'customName'};
                        
            % Create data structure to save to
            configData = struct();
            for f = fieldnames(handles)'
                switch f{:}
                    case ignored    % Do nothing
                    case strings    % Save string
                        configData.(f{:}) = get(handles.(f{:}), 'string');
                    case direct     % Save data
                        configData.(f{:}) = handles.(f{:});
                    otherwise       % Save value
                        configData.(f{:}) = get(handles.(f{:}), 'value');
                end
            end
            
            % Save data to .mat file
            save config.mat configData
        end
    end
end

function freqCoords = findFreqCoords(fourier, answers)
    % Finds the frequency coordinats at which the cell types are most different
    % 
    %   fourier   reshaped magnitude spectrum of many cells
    %   answers   corresponding cell types

    % Match answers to data
    a = find(strcmp(answers, 'apoptotic'));
    n = find(strcmp(answers, 'necrotic'));
    h = find(strcmp(answers, 'healthy'));
    i = find(strcmp(answers, 'ignored'));
    
    % average all from each set
    apopData = mean(fourier(a, :));
    necrData = mean(fourier(n, :));
    healData = mean(fourier(h, :));
    ignrData = mean(fourier(i, :));
    
    % Compare data sets
    an = apopData ./ necrData;
    ah = apopData ./ healData;
    ai = apopData ./ ignrData;
    na = necrData ./ apopData;
    nh = necrData ./ healData;
    ni = necrData ./ ignrData;
    ha = healData ./ apopData;
    hn = healData ./ necrData;
    hi = healData ./ ignrData;
    ia = ignrData ./ apopData;
    in = ignrData ./ necrData;
    ih = ignrData ./ healData;
    
    % Select points >= 85% of highest RQ differences
    vars = {an, ah, ai, na nh, ni, ha, hn, hi, ia, in, ih};
    combined = zeros(1, 0);
    for var = vars
        combined = [combined, var{:}(var{:} > 0.85*max(var{:}))]; %#ok<AGROW>
    end
    freqCoords = unique(combined);
end