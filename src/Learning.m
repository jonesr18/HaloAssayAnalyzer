function varargout = Learning(varargin)
% LEARNING MATLAB code for Learning.fig
%      LEARNING, by itself, creates a new LEARNING or raises the existing
%      singleton*.
%
%      H = LEARNING returns the handle to a new LEARNING or the handle to
%      the existing singleton*.
%
%      LEARNING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEARNING.M with the given input arguments.
%
%      LEARNING('Property','Value',...) creates a new LEARNING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Learning_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Learning_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
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

% Last Modified by GUIDE v2.5 29-Jul-2014 16:34:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Learning_OpeningFcn, ...
                   'gui_OutputFcn',  @Learning_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before Learning is made visible.
function Learning_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Learning (see VARARGIN)

% Load data and add to handles
data = load('dataForLearning.mat');
handles.cellNumber = 0;
handles.numCells = numel(data.aCells);
handles.aCells = data.aCells;
handles.labels = data.labels;
handles.stats = regionprops(data.labels, 'BoundingBox');
handles.image = data.image;
handles = updateImages(handles);

% Set image number
set(handles.imageText, 'string', sprintf('Image %d/%d', data.imNumber, data.numImages));

% Choose default command line output for Learning
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Learning wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Learning_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Updates the images on each "turn"
function handles = updateImages(handles)
% handles    structure with handles and user data (see GUIDATA)

% Update cell number
handles.cellNumber = handles.cellNumber + 1;
if handles.cellNumber > handles.numCells
    return
end
cellNum = handles.cellNumber;
set(handles.cellText, 'string', sprintf('Cell %d/%d', cellNum, handles.numCells));

% Draw box around the current cell
bb = handles.stats(cellNum).BoundingBox;
x1 = round(bb(1));
x2 = x1 + round(bb(3)) - 1;
y1 = round(bb(2));
y2 = y1 + round(bb(4)) - 1;
axes(handles.mainImage), imshow(handles.image, []), hold on
plot([x1 x2], [y1 y1], 'r', [x2 x2], [y1, y2], 'r', [x1 x2], [y2, y2], 'r', [x1 x1], [y1 y2], 'r')

% Show zoomed up current cell
axes(handles.zoomCell), imshow(handles.aCells{cellNum}.image, [])

% Calculate FFT of current cell image and show it
fourier = handles.aCells{cellNum}.fft();
axes(handles.magCell), imshow(fourier, [0, 5000])

% --- Closes the window, returning control to Initialization
function closeWindow(hObject, handles)
% hObject    handle to figure
% handles    structure with handles and user data (see GUIDATA)

aCells = handles.aCells;
magSpec = zeros(handles.numCells, 200^2);
allbins = zeros(handles.numCells, 100);
innerbins = zeros(handles.numCells, 100);
answers = cell(1, handles.numCells);
% 'reshape' back and forth between horizonal and normal data for ease
for i = 1:handles.numCells
    magSpec(i, :) = reshape(aCells{i}.fourier, 1, []);
    allbins(i, :) = reshape(aCells{i}.magbins, 1, []);
    innerbins(i, :) = reshape(aCells{i}.innerbins, 1, []);
    answers{i} = aCells{i}.cellType;
end
save learningResults.mat aCells magSpec allbins innerbins answers
close(hObject)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            %
%     CALLBACK FUNCTIONS     %
%                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --- Licence & Information --- &                                                               AL



% --- Executes on button press in about.
function about_Callback(~, ~, ~) %#ok<*DEFNU>
% hObject    handle to about (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
text = sprintf(['Copyright (C) 2014 Ross Jones, Singh Laboratory, University of Washington\n\n',...
                'Contact: jonesr18@gmail.com\n\n',...
                'This program is free software: you can redistribute it and/or ',...
                'modify it under the terms of the GNU General Public License as ',...
                'published by the Free Software Foundation, either version 3 of ',...
                'the License, or (at your option) any later version.\n\n',...
                'This program is distributed in the hope that it will be useful, ',...
                'but WITHOUT ANY WARRANTY; without even the implied warranty of ',...
                'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n\n',...
                'See the GNU General Public License for more details.\n\n',...
                'You should have received a copy of the GNU General Public ',...
                'License along with this program. If not see: ',...
                'http://www.gnu.org/licenses/gpl.html']);
msgbox(text, 'About', 'help')

% --- Executes on button press in licence.
function licence_Callback(~, ~, ~)
% hObject    handle to licence (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
open('GNU General Public Licence v3.0.txt')



% --- Buttons --- &                                                                             B



% --- Executes on button press in markHealthy.
function markHealthy_Callback(hObject, ~, handles)
% hObject    handle to markHealthy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.aCells{handles.cellNumber}.setType('healthy');
handles = updateImages(handles);
guidata(hObject, handles)
if handles.cellNumber > handles.numCells
    closeWindow(handles.figure1, handles)
end

% --- Executes on button press in markIgnore.
function markIgnore_Callback(hObject, ~, handles)
% hObject    handle to markIgnore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.aCells{handles.cellNumber}.setType('ignored');
handles = updateImages(handles);
guidata(hObject, handles)
if handles.cellNumber > handles.numCells
    closeWindow(handles.figure1, handles)
end

% --- Executes on button press in markNecr.
function markNecr_Callback(hObject, ~, handles)
% hObject    handle to markNecr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.aCells{handles.cellNumber}.setType('necrotic');
handles = updateImages(handles);
guidata(hObject, handles)
if handles.cellNumber > handles.numCells
    closeWindow(handles.figure1, handles)
end

% --- Executes on button press in markApop.
function markApop_Callback(hObject, ~, handles)
% hObject    handle to markApop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.aCells{handles.cellNumber}.setType('apoptotic');
handles = updateImages(handles);
guidata(hObject, handles)
if handles.cellNumber > handles.numCells
    closeWindow(handles.figure1, handles)
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          %
% --- Create Functions --- %
%                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
