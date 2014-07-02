%% DNA_damage
% A script for calculating DNA damage metrics in one or more image. I will have a GUI out soon,
% which should make this an easier process.
%
% !!! Avoid using CMYK color scheme images, as they will not work under the configuration here !!!
%
% With the method calcDamage(), you can change the argument 'all' to be one or more of the following
% (see writeup for info on what each one means):
%
%   'ndf'       Nuclear diffusion factor
%   'hd'        Halo distance
%   'hm'        Halo moment
%   'ahm'       Adjusted halo moment
%   'ihi'       Integrated halo intensity
%   'ahi'       Adjusted halo intensity
%   'rhi'       Relative halo intensity
%
% You can use the plot() method to see several things for each image. You can change the argument
% 'all' to one or more of the following:
%
%   'image'     Grayscale image     
%   'thresh'    Thresholded image
%   'binary'    Binary image
%   'halos'     Binary image of halos only
%   'nuclei'    Binary image of nuclei only
%   'eccen'     Outlined cells with eccentricity metric shown for each
%   'labels'    Cells labelled with different colors to show segmentation
%
% And you can use the plotDamage() method to plot the damage metrics on top of each cell. The
% argument 'all' can be changed to one or more of the same arguments passed to calcDamage() (see 
% above). If the numbers and words are too close together (likely to happen), maximize the image.
%
% As always, let me know if you have any issues.
%
% Ross Jones
% 4.12.14

close all
clear all
clc

% Import image(s) in a cell array
imArray = importImage();

% Perform image processing
for i = 1:numel(imArray)
    image = imArray{i};
    
    % Convert to RGB and pass over an averaging filter to smoothen image
    if size(image, 3) == 3
        image = rgb2gray(image);
    end
    image = imfilter(image, fspecial('average'));
    
    % Perform image processing and calculations
    cells = HaloImage(image);
    cells.process();
    cells.makeCells('damage');
    cells.calcDamage('all');            % --- Change 'all' if you want to specify metrics, such as ndf
    
    % Print out data
    fprintf('\n------------------------------------------------------------\nImage %d\n\n', i);
    names = fieldnames(cells.damage);
    for f = names'
        disp(f{:})
        disp(cells.damage.(f{:}))
    end
    
    % Plot images
    cells.plot('all')                   % --- Change 'all' if you want to specify images plotted
    cells.plotDamage('all')             % --- Change 'all' if you want to specify metrics displayed
end




