classdef Exports < handle
    % A functional class for saving data in standard CSV format from the FHA.
    %
    %   Methods (static):
    %
    %       writeCSV    - Exports data to a standard fast halo assay data file.
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
    % Updated 7.1.14
    
    methods (Static)
    %% public static void writeCSV(string, cell)
        function [] = writeCSV(fullfile, data)
            % Writes data to a standard fast halo assay data file. This method is static and thus
            % can be invoked without first creating a handle by calling Exports.writeCSV(...).
            %
            %   For an Exports object with handle 'out':
            %
            %   out.writeCSV(fullfile, data)
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
    end
end