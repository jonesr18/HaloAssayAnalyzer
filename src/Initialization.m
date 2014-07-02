function varargout = Initialization(varargin)
% INITIALIZATION MATLAB code for Initialization.fig
%      INITIALIZATION, by itself, creates a new INITIALIZATION or raises the existing
%      singleton*.
%
%      H = INITIALIZATION returns the handle to a new INITIALIZATION or the handle to
%      the existing singleton*.
%
%      INITIALIZATION('CALLBACK',hObject,~,handles,...) calls the local
%      function named CALLBACK in INITIALIZATION.M with the given input arguments.
%
%      INITIALIZATION('Property','Value',...) creates a new INITIALIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Initialization_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Initialization_OpeningFcn via varargin.
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
% Updated 7.1.14

% Last Modified by GUIDE v2.5 01-Jul-2014 11:15:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Initialization_OpeningFcn, ...
                   'gui_OutputFcn',  @Initialization_OutputFcn, ...
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


% --- Executes just before Initialization is made visible.
function Initialization_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Initialization (see VARARGIN)

warning('off', 'images:initSize:adjustingMag');

% Set handles fields
handles.images = {};
handles.currImage = [];
handles.outdir = strcat(pwd, '\');
handles.customName = false;
set(handles.lowerSlider, 'value', 10); % needed to make it go to the bottom for some reason
set(handles.lowerSlider, 'value', 0);
set(handles.upperSlider, 'value', 255);
set(handles.lowerThresh, 'string', num2str(0));
set(handles.upperThresh, 'string', num2str(255));
set(handles.outputPath, 'string', handles.outdir);
set(handles.filterEccen, 'value', 1);   % By default, filter cells by eccentricity

% Choose default command line output for Initialization
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Initialization wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Initialization_OutputFcn(~, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



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


% --- Prepare Images --- %                                                                      PI



% --- Executes on button press in imageBrowse.
function imageBrowse_Callback(hObject, ~, handles)
% hObject    handle to imageBrowse (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
imArray = Imports.images();
menuItems = cell(1, numel(imArray));
handles.images = cell(1, numel(imArray));
for i = 1:numel(imArray)
    handles.images{i} = HaloImage(imArray{i});
    menuItems{i} = sprintf('Image %d', i);
end
set(handles.imageMenu, 'String', char(menuItems));
handles.currImage = 1;
updateImage(handles)
guidata(hObject, handles)

% --- Executes on selection change in imageMenu.
function imageMenu_Callback(hObject, ~, handles)
% hObject    handle to imageMenu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns imageMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imageMenu
handles.currImage = get(hObject, 'Value');
updateImage(handles)
guidata(hObject, handles)

% --- Executes on button press in crop.
function crop_Callback(hObject, ~, handles)
% hObject    handle to crop (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[~, rect] = imcrop(handles.images{handles.currImage}.image, []);
if ~isempty(rect)
    x1 = round(rect(1));
    x2 = x1 + round(rect(3)) - 1;
    y1 = round(rect(2));
    y2 = y1 + round(rect(4)) - 1;
    for im = handles.images
        im{:}.crop(y1, y2, x1, x2);
    end
end
updateImage(handles)
guidata(hObject, handles)

% --- Executes on button press in invert.
function invert_Callback(hObject, ~, handles)
% hObject    handle to invert (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for im = handles.images
    im{:}.invert();
end
updateImage(handles)
guidata(hObject, handles)

% --- Executes on slider movement.
function lowerSlider_Callback(hObject, ~, handles)
% hObject    handle to lowerSlider (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = round(get(hObject, 'value'));
set(hObject, 'value', val);
set(handles.lowerThresh, 'string', num2str(val));
updateImage(handles)
guidata(hObject, handles)

% --- Executes on slider movement.
function upperSlider_Callback(hObject, ~, handles)
% hObject    handle to upperSlider (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = round(get(hObject, 'value'));
set(hObject, 'value', val);
set(handles.upperThresh, 'string', num2str(val));
updateImage(handles)
guidata(hObject, handles)

function lowerThresh_Callback(hObject, ~, handles)
% hObject    handle to lowerThresh (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowerThresh as text
%        str2double(get(hObject,'String')) returns contents of lowerThresh as a double
val = uint8(str2double(get(hObject, 'string')));
set(hObject, 'string', num2str(val));
set(handles.lowerSlider, 'value', val);
updateImage(handles)
guidata(hObject, handles)

function upperThresh_Callback(hObject, ~, handles)
% hObject    handle to upperThresh (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upperThresh as text
%        str2double(get(hObject,'String')) returns contents of upperThresh as a double
val = uint8(str2double(get(hObject, 'string')));
set(hObject, 'string', num2str(val));
set(handles.upperSlider, 'value', val);
updateImage(handles)
guidata(hObject, handles)

% --- Executes on button press in manualThresh.
function manualThresh_Callback(hObject, ~, handles)
% hObject    handle to manualThresh (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of manualThresh
updateImage(handles)
guidata(hObject, handles)

function updateImage(handles)
% Updates the current displayed image. 
% If manual thresh is selected, thresholds are overlaid on the image. 
image = handles.images{handles.currImage}.image;
imshow(image, [])
if get(handles.manualThresh, 'value')
    hold on
    levels = [get(handles.lowerSlider, 'value'), get(handles.upperSlider, 'value')];
    thresh = imquantize(image, levels);
    himage = imshow(label2rgb(thresh));
    set(himage, 'AlphaData', 0.3);
    hold off
end



% --- DNA Damage Metrics: Calculate --- %                                                       DDMc



% --- Executes on button press in calcAll.
function calcAll_Callback(hObject, ~, handles)
% hObject    handle to calcAll (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calcAll
val = get(hObject, 'value');
set(handles.calcNDF, 'value', val);
set(handles.calcHD, 'value', val);
set(handles.calcHM, 'value', val);
set(handles.calcAHM, 'value', val);
set(handles.calcIHI, 'value', val);
set(handles.calcAHI, 'value', val);
set(handles.calcRHI, 'value', val);
if ~get(hObject, 'value')
    set(handles.plotAllDmg, 'value', 0)
    plotAllDmg_Callback(handles.plotAllDmg, [], handles);
end
guidata(hObject, handles)

% --- Executes on button press in calcNDF.
function calcNDF_Callback(hObject, ~, handles)
% hObject    handle to calcNDF (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calcNDF
if ~get(hObject, 'value')
    if get(handles.calcAll, 'value')
        set(handles.calcAll, 'value', 0);
    end
    set(handles.plotNDF, 'value', 0)
    plotNDF_Callback(handles.plotNDF, [], handles);
end
guidata(hObject, handles)

% --- Executes on button press in calcHD.
function calcHD_Callback(hObject, ~, handles)
% hObject    handle to calcHD (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calcHD
if ~get(hObject, 'value')
    if get(handles.calcAll, 'value')
        set(handles.calcAll, 'value', 0);
    end
    set(handles.plotHD, 'value', 0)
    plotHD_Callback(handles.plotHD, [], handles);
end
guidata(hObject, handles)

% --- Executes on button press in calcHM.
function calcHM_Callback(hObject, ~, handles)
% hObject    handle to calcHM (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calcHM
if ~get(hObject, 'value')
    if get(handles.calcAll, 'value')
        set(handles.calcAll, 'value', 0);
    end
    set(handles.plotHM, 'value', 0)
    plotHM_Callback(handles.plotHM, [], handles);
end
guidata(hObject, handles)

% --- Executes on button press in calcAHM.
function calcAHM_Callback(hObject, ~, handles)
% hObject    handle to calcAHM (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calcAHM
if ~get(hObject, 'value')
    if get(handles.calcAll, 'value')
        set(handles.calcAll, 'value', 0);
    end
    set(handles.plotAHM, 'value', 0)
    plotAHM_Callback(handles.plotAHM, [], handles);
end
guidata(hObject, handles)

% --- Executes on button press in calcIHI.
function calcIHI_Callback(hObject, ~, handles)
% hObject    handle to calcIHI (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calcIHI
if ~get(hObject, 'value')
    if get(handles.calcAll, 'value')
        set(handles.calcAll, 'value', 0);
    end
    set(handles.plotIHI, 'value', 0)
    plotIHI_Callback(handles.plotIHI, [], handles);
end
guidata(hObject, handles)

% --- Executes on button press in calcAHI.
function calcAHI_Callback(hObject, ~, handles)
% hObject    handle to calcAHI (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calcAHI
if ~get(hObject, 'value')
    if get(handles.calcAll, 'value')
        set(handles.calcAll, 'value', 0);
    end
    set(handles.plotAHI, 'value', 0)
    plotAHI_Callback(handles.plotAHI, [], handles);
end
guidata(hObject, handles)

% --- Executes on button press in calcRHI.
function calcRHI_Callback(hObject, ~, handles)
% hObject    handle to calcRHI (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of calcRHI
if ~get(hObject, 'value')
    if get(handles.calcAll, 'value')
        set(handles.calcAll, 'value', 0);
    end
    set(handles.plotRHI, 'value', 0)
    plotRHI_Callback(handles.plotRHI, [], handles);
end
guidata(hObject, handles)



% --- DNA Damage Metrics: Plot --- %                                                            DDMp



% --- Executes on button press in plotAllDmg.
function plotAllDmg_Callback(hObject, ~, handles)
% hObject    handle to plotAllDmg (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotAllDmg
val = get(hObject, 'value');
set(handles.plotNDF, 'value', val);
plotNDF_Callback(handles.plotNDF, [], handles);

set(handles.plotHD, 'value', val);
plotHD_Callback(handles.plotHD, [], handles)

set(handles.plotHM, 'value', val);
plotHM_Callback(handles.plotHM, [], handles)

set(handles.plotAHM, 'value', val);
plotAHM_Callback(handles.plotAHM, [], handles)

set(handles.plotIHI, 'value', val);
plotIHI_Callback(handles.plotIHI, [], handles)

set(handles.plotAHI, 'value', val);
plotAHI_Callback(handles.plotAHI, [], handles)

set(handles.plotRHI, 'value', val);
plotRHI_Callback(handles.plotRHI, [], handles)

if ~get(handles.calcAll, 'value')
    set(hObject, 'value', 0);
end
guidata(hObject, handles)

% --- Executes on button press in plotNDF.
function plotNDF_Callback(hObject, ~, handles)
% hObject    handle to plotNDF (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotNDF
if ~get(handles.calcNDF, 'value')
    set(hObject, 'value', 0);
end
if ~get(hObject, 'value')
    if get(handles.plotAllDmg, 'value')
        set(handles.plotAllDmg, 'value', 0);
    end
end
guidata(hObject, handles)

% --- Executes on button press in plotHD.
function plotHD_Callback(hObject, ~, handles)
% hObject    handle to plotHD (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotHD
if ~get(handles.calcHD, 'value')
    set(hObject, 'value', 0);
end
if ~get(hObject, 'value')
    if get(handles.plotAllDmg, 'value')
        set(handles.plotAllDmg, 'value', 0);
    end
end
guidata(hObject, handles)

% --- Executes on button press in plotHM.
function plotHM_Callback(hObject, ~, handles)
% hObject    handle to plotHM (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotHM
if ~get(handles.calcHM, 'value')
    set(hObject, 'value', 0);
end
if ~get(hObject, 'value')
    if get(handles.plotAllDmg, 'value')
        set(handles.plotAllDmg, 'value', 0);
    end
end
guidata(hObject, handles)

% --- Executes on button press in plotAHM.
function plotAHM_Callback(hObject, ~, handles)
% hObject    handle to plotAHM (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotAHM
if ~get(handles.calcAHM, 'value')
    set(hObject, 'value', 0);
end
if ~get(hObject, 'value')
    if get(handles.plotAllDmg, 'value')
        set(handles.plotAllDmg, 'value', 0);
    end
end
guidata(hObject, handles)

% --- Executes on button press in plotIHI.
function plotIHI_Callback(hObject, ~, handles)
% hObject    handle to plotIHI (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotIHI
if ~get(handles.calcIHI, 'value')
    set(hObject, 'value', 0);
end
if ~get(hObject, 'value')
    if get(handles.plotAllDmg, 'value')
        set(handles.plotAllDmg, 'value', 0);
    end
end
guidata(hObject, handles)

% --- Executes on button press in plotAHI.
function plotAHI_Callback(hObject, ~, handles)
% hObject    handle to plotAHI (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotAHI
if ~get(handles.calcAHI, 'value')
    set(hObject, 'value', 0);
end
if ~get(hObject, 'value')
    if get(handles.plotAllDmg, 'value')
        set(handles.plotAllDmg, 'value', 0);
    end
end
guidata(hObject, handles)

% --- Executes on button press in plotRHI.
function plotRHI_Callback(hObject, ~, handles)
% hObject    handle to plotRHI (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotRHI
if ~get(handles.calcRHI, 'value')
    set(hObject, 'value', 0);
end
if ~get(hObject, 'value')
    if get(handles.plotAllDmg, 'value')
        set(handles.plotAllDmg, 'value', 0);
    end
end
guidata(hObject, handles)



% --- Image Plots --- %                                                                         IP



% --- Executes on button press in plotAllImage.
function plotAllImage_Callback(hObject, ~, handles)
% hObject    handle to plotAllImage (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotAllImage
val = get(hObject, 'value');
set(handles.plotImage, 'value', val);
set(handles.plotContour, 'value', val);
set(handles.plotThresh, 'value', val);
set(handles.plotBinary, 'value', val);
set(handles.plotEccen, 'value', val);
set(handles.plotGradmag, 'value', val);
set(handles.plotCircles, 'value', val);
set(handles.plotLabels, 'value', val);
plotCircles_Callback(handles.plotCircles, [], handles)
guidata(hObject, handles)

% --- Executes on button press in plotImage.
function plotImage_Callback(hObject, ~, handles)
% hObject    handle to plotImage (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotImage
if ~get(hObject, 'value')
    if get(handles.plotAllImage, 'value')
        set(handles.plotAllImage, 'value', 0);
    end
end
guidata(hObject, handles)

% --- Executes on button press in plotContour.
function plotContour_Callback(hObject, ~, handles)
% hObject    handle to plotContour (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotContour
if ~get(hObject, 'value')
    if get(handles.plotAllImage, 'value')
        set(handles.plotAllImage, 'value', 0);
    end
end
guidata(hObject, handles)

% --- Executes on button press in plotThresh.
function plotThresh_Callback(hObject, ~, handles)
% hObject    handle to plotThresh (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotThresh
if ~get(hObject, 'value')
    if get(handles.plotAllImage, 'value')
        set(handles.plotAllImage, 'value', 0);
    end
end
guidata(hObject, handles)

% --- Executes on button press in plotBinary.
function plotBinary_Callback(hObject, ~, handles)
% hObject    handle to plotBinary (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotBinary
if ~get(hObject, 'value')
    if get(handles.plotAllImage, 'value')
        set(handles.plotAllImage, 'value', 0);
    end
end
guidata(hObject, handles)

% --- Executes on button press in plotEccen.
function plotEccen_Callback(hObject, ~, handles)
% hObject    handle to plotEccen (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotEccen
if ~get(handles.filterEccen, 'value')
    set(hObject, 'value', 0);
end
if ~get(hObject, 'value')
    if get(handles.plotAllImage, 'value')
        set(handles.plotAllImage, 'value', 0);
    end    
end
guidata(hObject, handles)

% --- Executes on button press in plotGradmag.
function plotGradmag_Callback(hObject, ~, handles)
% hObject    handle to plotGradmag (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotGradmag
if ~get(hObject, 'value')
    if get(handles.plotAllImage, 'value')
        set(handles.plotAllImage, 'value', 0);
    end
end
guidata(hObject, handles)

% --- Executes on button press in plotCircles.
function plotCircles_Callback(hObject, ~, handles)
% hObject    handle to plotCircles (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotCircles
contents = cellstr(get(handles.mode, 'string'));
selected = contents{get(handles.mode, 'value')};
if ~strcmp(selected, 'Fit Circles')
    set(hObject, 'value', 0);
end
if ~get(hObject, 'value')
    if get(handles.plotAllImage, 'value')
        set(handles.plotAllImage, 'value', 0);
    end
end
guidata(hObject, handles)

% --- Executes on button press in plotLabels.
function plotLabels_Callback(hObject, ~, handles)
% hObject    handle to plotLabels (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotLabels
if ~get(hObject, 'value')
    if get(handles.plotAllImage, 'value')
        set(handles.plotAllImage, 'value', 0);
    end
end
guidata(hObject, handles)



% --- Segmentation --- %                                                                        S



% --- Executes on selection change in mode.
function mode_Callback(hObject, ~, handles)
% hObject    handle to mode (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns mode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mode
contents = cellstr(get(hObject, 'string'));
selected = contents{get(hObject, 'value')};
if ~strcmp(selected, 'Fit Circles')
    set(handles.plotCircles, 'value', 0);
    plotCircles_Callback(handles.plotCircles, [], handles);
end
guidata(hObject, handles)

% --- Executes on selection change in scheme.
function scheme_Callback(hObject, ~, handles)
% hObject    handle to scheme (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns scheme contents as cell array
%        contents{get(hObject,'Value')} returns selected item from scheme
guidata(hObject, handles)

% --- Executes on button press in filterEccen.
function filterEccen_Callback(hObject, ~, handles)
% hObject    handle to filterEccen (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filterEccen
if ~get(hObject, 'value')
    set(handles.plotEccen, 'value', 0);
    plotEccen_Callback(handles.plotEccen, [], handles);
end
guidata(hObject, handles)

function eccenMetric_Callback(hObject, ~, handles)
% hObject    handle to eccenMetric (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of eccenMetric as text
%        str2double(get(hObject,'String')) returns contents of eccenMetric as a double
val = str2double(get(hObject, 'string'));
if val > 1
    val = 1;
elseif val < 0
    val = 0;
end
set(hObject, 'string', sprintf('%.3f', val));
guidata(hObject, handles)



% --- Cell Classification --- %                                                                 CC



% --- Executes on selection change in classifier.
function classifier_Callback(hObject, ~, handles)
% hObject    handle to classifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns classifier contents as cell array
%        contents{get(hObject,'Value')} returns selected item from classifier
contents = cellstr(get(hObject, 'string'));
if strcmp(contents{get(hObject, 'value')}, 'Select Classifier')
    % Off state
    set(handles.plotClass, 'value', 0);
end
guidata(hObject, handles)

% --- Executes on button press in plotClass.
function plotClass_Callback(hObject, ~, handles)
% hObject    handle to plotClass (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotClass
contents = cellstr(get(handles.classifier, 'string'));
if strcmp(contents{get(handles.classifier, 'value')}, 'Select Classifier')
    % Off state
    set(hObject, 'value', 0);
end
guidata(hObject, handles)



% --- Output Directory --- %                                                                    OD



function outputPath_Callback(hObject, ~, handles)
% hObject    handle to outputPath (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outputPath as text
%        str2double(get(hObject,'String')) returns contents of outputPath as a double
outdir = get(hObject, 'string');
if exist(outdir, 'dir')
    if outdir(end) ~= '\'
        outdir = strcat(outdir, '\');
    end
    handles.outdir = outdir;
end
set(hObject, 'string', handles.outdir);
guidata(hObject, handles)

% --- Executes on button press in outputBrowse.
function outputBrowse_Callback(hObject, ~, handles)
% hObject    handle to outputBrowse (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
outdir = uigetdir(pwd, 'Output File Directory');
if ischar(outdir)
    set(handles.outputPath, 'string', outdir);
    outputPath_Callback(handles.outputPath, [], handles);
end
guidata(hObject, handles)




% --- Sample Info --- %                                                                         SI



function sampleName_Callback(hObject, ~, handles)
% hObject    handle to sampleName (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sampleName as text
%        str2double(get(hObject,'String')) returns contents of sampleName as a double
handles.customName = true;
guidata(hObject, handles)



% --- Run Program --- %                                                                         Run



% --- Executes on button press in run.
function run_Callback(~, ~, handles)
% hObject    handle to run (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

metrics = {'NDF', 'HD', 'HM', 'AHM', 'IHI', 'AHI', 'RHI'};
plotTypes = {'Image', 'Contour', 'Thresh', 'Binary', 'Eccen', 'Gradmag', 'Circles', 'Labels'};

% Check DNA damage calculation options
damageCalcs = {};
index = 1;
if get(handles.calcAll, 'value')
    damageCalcs = {'all'};
else
    for metric = metrics
        if get(handles.(sprintf('calc%s', metric{:})), 'value')
            damageCalcs{index} = metric{:}; %#ok<AGROW>
            index = index + 1;
        end
    end
end

% Check DNA damage plotting options
damagePlots = {};
index = 1;
if get(handles.plotAllDmg, 'value')
    damagePlots = {'all'};
else
    for metric = metrics
        if get(handles.(sprintf('plot%s', metric{:})), 'value')
            damagePlots{index} = metric{:}; %#ok<AGROW>
            index = index + 1;
        end
    end
end

% Check image plotting options
imagePlots = {};
index = 1;
if get(handles.plotAllImage, 'value')
    imagePlots = {'all'};
else
    for type = plotTypes
        if get(handles.(sprintf('plot%s', type{:})), 'value')
            imagePlots{index} = type{:}; %#ok<AGROW>
            index = index + 1;
        end
    end
end

% Check for manual threshold
if get(handles.manualThresh, 'value')
    for im = handles.images
        im{:}.setThresh(get(handles.lowerSlider, 'value'), get(handles.upperSlider, 'value'));
    end
end

% Check segmentation scheme
contents = cellstr(get(handles.mode, 'string'));
switch contents{get(handles.mode, 'value')};
    case {'Binary Image'}
        mode = 'bw';
    case {'Fit Circles'}
        mode = 'fit';
end
contents = cellstr(get(handles.scheme,'string'));
switch contents{get(handles.scheme, 'value')};
    case {'Threshold'}
        scheme = 'thresh';
    case {'Gradient Mag'}
        scheme = 'grad';
end
if get(handles.filterEccen, 'value')
    filterEccen = true;
    for im = handles.images
        im{:}.setCriteria(str2double(get(handles.eccenMetric, 'string')));
    end
else
    filterEccen = false;
end

% Prepare cell matrix for DNA damage Excel output
%   DD is for "DNA Damage" hereout
if ~isempty(damageCalcs) 
    if strcmp(damageCalcs{1}, 'all')
        metricsUsed = metrics;
    else
        metricsUsed = damageCalcs;
    end
    titlesDD = [{'Image #', 'Cell #'}, metricsUsed];
    combinedDD = cell(1, numel(handles.images));
    numCellsDD = 0; 
    allStatsDataDD = cell(0, numel(metricsUsed));
end

% Prepare for cell classification output
%   CC is for "Cell Classification" hereout
contents = cellstr(get(handles.classifier, 'string'));
classifier = contents{get(handles.classifier, 'value')};
if ~strcmp(classifier, 'Select Classifier')
    
    % Check which classifier to use
    classify = true;
    if strncmp(classifier, 'k-NN', 4)
        classifier = 'knn';
    elseif strcmp(classifier, 'k-Means')
        classifier = 'kmeans';
    elseif strncmp(classifer, 'Classif', 7)
        classifier = 'ctree';
    else
        warning('Classifier not correctly found: %s', classifier)
    end
    
    % Prepare cell matrix
    cellTypes = {'Apoptotic', 'Necrotic', 'Healthy'};
    titlesCC = [{'Image #', 'Cell #'}, cellTypes];
    combinedCC = cell(1, numel(handles.images)); 
    numCellsCC = 0; 
    allStatsDataCC = num2cell(zeros(1, numel(cellTypes)));
else
    classify = false;
end


% Cycle through each iamge and perform actions as requested
%   Process images, calculate metrics, plot images, classify cells
for i = 1:numel(handles.images)
    
    % Process image
    im = handles.images{i};
    im.process();
    
    % Make apoptosis cell objects and analyze if requested
    if classify
        im.makeCells('apop', mode, scheme, filterEccen);
        if ~isempty(im.aCells)
            am = load('non_adj_H2O2_allMags.mat');
            an = load('non_adj_H2O2_answers.mat');
            im.classify(am.allMags, an.answers, classifier);
        else
            warning('No cells found in image %d for scoring apoptosis', i);
        end
        
        % Collect image data
        imageDataCC = cell(0, numel(titlesCC));
        j = 0;
        for c = im.aCells
            if any(strcmpi(c{:}.cellType, cellTypes))
                j = j + 1;
                numCellsCC = numCellsCC + 1;
                
                % Get cell type
                typeData = num2cell(double(strcmpi(cellTypes, c{:}.cellType)));

                % Store data
                imageDataCC = [imageDataCC; [{'', sprintf('Cell %d', j)}, typeData]]; %#ok<AGROW>
            end
        end
        countedCells = j;

        % Collect image stats
        statsLabelsCC = {
            sprintf('Image %d', i), 'num';
            sprintf('n = %d', j), '%'};
        statsDataCC = cell(2, numel(cellTypes));
        for j = 1:numel(cellTypes)
            if isempty(im.aCells)
                statsDataCC{1, j} = 'N/A';
                statsDataCC{2, j} = 'N/A';
            else
                statsDataCC{1, j} = im.cellTypes.(lower(cellTypes{j}));
                allStatsDataCC{1, j} = allStatsDataCC{1, j} + statsDataCC{1, j};
                statsDataCC{2, j} = im.cellTypes.(lower(cellTypes{j})) / countedCells * 100;
            end
        end

        % Merge stats and data
        combinedCC{i} = [
            statsLabelsCC, statsDataCC; 
            imageDataCC];
    end
    
    % Make damage cell objects and analyze if requested
    if ~isempty(damageCalcs)
        im.makeCells('dmg', mode, scheme, filterEccen);
        if ~isempty(im.dCells)
            im.calcDamage(damageCalcs{:});
        else
            warning('No cells found in image %d for scoring DNA damage', i);
        end
        
        % Collect image data
        imageDataDD = cell(numel(im.dCells), numel(titlesDD));
        for j = 1:numel(im.dCells)
            c = im.dCells{j};
            numCellsDD = numCellsDD + 1;

            % Get damage values
            damageData = cell(1, numel(metricsUsed));
            for k = 1:numel(metricsUsed)
                damageData{k} = c.(lower(metricsUsed{k}));
            end

            % Store data
            allStatsDataDD = [allStatsDataDD; damageData]; %#ok<AGROW>
            imageDataDD(j, :) = [{'', sprintf('Cell %d', j)}, damageData];
        end

        % Collect image stats
        statsLabelsDD = {
            sprintf('Image %d', i), 'mean';
            sprintf('n = %d', j), 'stdev';
            '', 'SEM'};
        statsDataDD = cell(3, numel(metricsUsed));
        for j = 1:numel(metricsUsed)
            if isempty(im.dCells)
                statsDataDD{1, j} = 'N/A';
                statsDataDD{2, j} = 'N/A';
                statsDataDD{3, j} = 'N/A';
            else
                statsDataDD{1, j} = im.damage.(lower(metricsUsed{j})).Mean;
                statsDataDD{2, j} = im.damage.(lower(metricsUsed{j})).Stdev;
                statsDataDD{3, j} = im.damage.(lower(metricsUsed{j})).SEM;
            end
        end

        % Merge stats and data
        combinedDD{i} = [
            statsLabelsDD, statsDataDD; 
            imageDataDD];
    end
    
    % Plot requested figures
    if ~isempty(damagePlots)
        im.plotDamage(damagePlots{:})
    end
    if ~isempty(imagePlots)
        im.plot(imagePlots{:})
    end
    if get(handles.plotClass, 'value')
        im.plotTypes();
    end
end

% Combine overall statistics for cell classification
if classify
    
    % Prepare statistics labels
    allStatsLabelsCC = {
        'All Images', 'num';
        sprintf('n = %d', numCellsCC), '%'};
    
    % Finalize all cell statistics data
    percentCellType = num2cell(cell2mat(allStatsDataCC) / numCellsCC * 100);
    allStatsDataCC = [allStatsDataCC; percentCellType];
    
    % Combine all data
    allDataCC = [
        titlesCC;
        allStatsLabelsCC, allStatsDataCC];
    for extract = combinedCC
        allDataCC = [allDataCC; extract{:}]; %#ok<AGROW>
    end
end

% Compute overall statistics for DNA damage
if ~isempty(damageCalcs)
    
    % Prepare statistics labels
    allStatsLabelsDD = {
        'All Images', 'mean';
        sprintf('n = %d', numCellsDD), 'stdev'
        '', 'SEM'};
    
    % Extract raw data from each cell and compute statistics
    statsDataDD = cell(3, numel(metricsUsed));
    for i = 1:numel(metricsUsed)
        col = cell2mat(allStatsDataDD(:, i));
        statsDataDD{1, i} = mean(col(:));
        statsDataDD{2, i} = std(col(:));
        statsDataDD{3, i} = std(col(:)) / sqrt(numCellsDD);
    end

    % Combine all data
    allDataDD = [
        titlesDD;
        allStatsLabelsDD, statsDataDD];
    for extract = combinedDD
        allDataDD = [allDataDD; extract{:}]; %#ok<AGROW>
    end
end


% Make filenames for CSV output
if handles.customName
    filename = get(handles.sampleName, 'string');
    if exist(strcat(handles.outdir, filename, ' damage.csv'), 'file') ||...
       exist(strcat(handles.outdir, filename, ' classification.csv'), 'file')
        filename = strcat(filename, ' X');
    end
else
    f = datestr(now);
    f(f == ':') = ';';
    f(f == '-') = ' ';
    filename = f;
end
fullfileDD = strcat(handles.outdir, filename, ' damage', '.csv');
fullfileCC = strcat(handles.outdir, filename, ' classification', '.csv');

% Write files to CSV in standard output
out = Exports();
out.writeCSV(fullfileDD, allDataDD);
out.writeCSV(fullfileCC, allDataCC);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                          %
% --- Create Functions --- %
%                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% --- Executes during object creation, after setting all properties.
function lowerSlider_CreateFcn(hObject, ~, handles)
% hObject    handle to lowerSlider (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function upperSlider_CreateFcn(hObject, ~, handles)
% hObject    handle to upperSlider (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function sampleName_CreateFcn(hObject, ~, handles)
% hObject    handle to sampleName (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function sampleNumber_CreateFcn(hObject, ~, handles)
% hObject    handle to sampleNumber (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function outputPath_CreateFcn(hObject, ~, handles)
% hObject    handle to outputPath (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function imageMenu_CreateFcn(hObject, ~, handles)
% hObject    handle to imageMenu (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function upperThresh_CreateFcn(hObject, ~, handles)
% hObject    handle to upperThresh (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function lowerThresh_CreateFcn(hObject, ~, handles)
% hObject    handle to lowerThresh (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function eccenMetric_CreateFcn(hObject, ~, handles)
% hObject    handle to eccenMetric (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function segmentationPanel_CreateFcn(hObject, ~, handles)
% hObject    handle to segmentationPanel (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function scheme_CreateFcn(hObject, ~, handles)
% hObject    handle to scheme (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function mode_CreateFcn(hObject, ~, handles)
% hObject    handle to mode (see GCBO)
% ~  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles)

% --- Executes during object creation, after setting all properties.
function classifier_CreateFcn(hObject, ~, handles)
% hObject    handle to classifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
guidata(hObject, handles)
