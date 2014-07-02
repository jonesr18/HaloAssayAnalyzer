%% Segmentation testing
close all
clear all
clc

warning('off', 'images:initSize:adjustingMag');
image = imread('C:\Users\Ross\Google Drive\Capstone Images\Binh\Array Image-1.tif');
image = imfilter(rgb2gray(image), fspecial('gaussian'));
cells = HaloImage(image);
cells.process();
cells.makeCells('dmg');


%% Hanning window (Tukey)
close all
clc

stats = regionprops(cells.dLabels, 'Area', 'Centroid');
for i = 1:max(cells.dLabels(:)); 
    % Estimate window length and force it to be even
    len = round(sqrt(stats(i).Area) * 1.5);
    if rem(len, 2) == 1
        len = len + 1;
    end
    
    % Create tukey hanning window
    window = tukeywin(len, 0.4);
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
    cellRegion = cells.image(y1:y2, x1:x2) .* uint8(win2D);
    
    % Plot results
    figure()
    subplot(1, 3, 1), imshow(win2D)
    subplot(1, 3, 2), imshow(cells.image(y1:y2, x1:x2), [])
    subplot(1, 3, 3), imshow(cellRegion, [])
end


%% Edge detection of nuclei for segmentations
close all
clc

% Get gradient magnitude
hy = fspecial('sobel');
Iy = imfilter(double(cells.image), hy, 'replicate');
Ix = imfilter(double(cells.image), hy', 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
% self.gradmag = gradmag;
figure(), imshow(gradmag, [min(gradmag(:)), max(gradmag(:)/4)]), title('Gradient Magnitude')

% Threshold gradient magnitude
imthresh = imquantize(gradmag, multithresh(gradmag, 2));
gradient  = boolean(imthresh > 2);
figure(), imshow(gradient), title('Thresholded Gradient Magnitude')

% Dilate the image -- closing and opening smoothes edges
strel1   = strel('disk', 1);
gradient = imopen(gradient, strel1);
gradient = imdilate(gradient, strel('disk', 3));
gradient = imdilate(gradient, strel('diamond', 3));
gradient = imopen(gradient, strel1);
gradient = imclose(gradient, strel1);
gradient = imopen(gradient, strel1);
gradient = imclose(gradient, strel1); 
% figure(), imshow(gradient), title('Dilated Gradient Mask')

% Fil interior gaps
gradient = imfill(gradient, 'holes');
% figure(), imshow(gradient), title('Binary Image with Filled Holes')

% Remove connected objects on border
gradient = imclearborder(gradient, 4);
figure(), imshow(gradient), title('Segmented Image')

% Remove small objects and outline remaining ones
[boundaries, labels] = bwboundaries(gradient);
% figure(), imshow(gradmag), title('Outlined Original Image'), hold on
stats = regionprops(labels, 'Area');
for i = 1:length(stats)
    if stats(i).Area < 300
        gradient = gradient - boolean(labels == i);
    else
%         boundary = boundaries{i};
%         plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2)
    end
end
% hold off

% Collect stats on newly processed gradient
binGradmag = boolean(gradmag);
statsNew = regionprops(binGradmag, 'EquivDiameter');
radList = [statsNew.EquivDiameter] / 2;
rad = round(mean(radList));
radstd = round(std(radList));

% Fit circles to the image to find individual nuclei
% NOTE: Circle finding is difficult with low-resolution pictures. It also had
%       issues when CMYK was converted to RGB due to poor contrast.
try 
    range = [rad - radstd, rad + radstd];
    [cn, r] = imfindcircles(gradient, range, 'ObjectPolarity', 'bright',...
                            'Sensitivity', 0.92);
    if isempty(cn)
        warning('No circles found using gradient, attempting with nuclei...')
        error('not disp')
    end
catch
    % In case the gradient magnitude wasn't strong enough or just didn't work,
    % use the binary nuclei image instead.
    statsCatch = regionprops(cells.nuclei, 'EquivDiameter');
    radList = [statsCatch.EquivDiameter] / 2;
    rad = round(median(radList));
    radstd = round(std(radList) / 2);
    range = [rad - radstd, rad + radstd];
    [cn, r] = imfindcircles(self.nuclei, range, 'ObjectPolarity', 'bright',...
                            'Sensitivity', 0.92);
    if isempty(cn)
        error('Could not fit circles to image, try a higher resolution')
    end
end
% self.centroids = cn;
% self.radii = r;

% Combine centroids within 20 pixels (overfitted circles)
allDistances = dist(cn');
[x, y] = find(allDistances < 20);
indexes = find(x ~= y);
for i = indexes'
    average = mean([cn(x(i), :); cn(y(i), :)]);
    cn(x(i), :) = average;
    cn(y(i), :) = average;
end
                    
% Match pixels to closest centroid
blank = zeros(size(cells.binary));
for c = cn'
    x = round(c(1));
    y = round(c(2));
    blank(y, x) = 1;
end

% Perform segmentation
labels = watershed(bwdist(blank));
labels = labels .* imclearborder(uint8(self.binary), 4);
% self.aLabels = labels;

% Plot estimated circles on image
low  = min(gradmag(:));
high = max(gradmag(:)/4);
figure(), imshow(gradmag, [low, high]), title('Gradmag + Circles'), viscircles(centroids, radii);
figure(), imshow(cells.image, []), title('Original Image + Circles'), viscircles(centroids, radii);
figure(), imshow(label2rgb(labels)), title('Labelled Cells')


%% Watershed transform -- Run the following before running the subsequent sections
close all
clc

% Pre-process image
imRegions = processImage(image, .70, 2);
binary = imRegions(1).Nuclei;
figure(), imshow(image, []), title('Pre-processed Image')

% May be a good idea to individually threshold each cell CLUSTER instead of whole image


%% Segment image with circleSegment
close all
clc

labels = circleSegment(imRegions, true);


%% Normal method
% Using nuclei allows for ID of more cells here, though extrapolation will be difficult
close all
clc

% Compute the distance transform of the complement of the binary image. 
distance = bwdist(~binary);
figure(), imshow(distance, []), title('Distance transform of ~bw')
 
% Complement the distance transform, and force pixels that don't belong to the objects to be at -Inf.
distance = -distance;
distance(~binary) = -Inf;
 
% Compute the watershed transform, and display the resulting label matrix as an RGB image.
L = watershed(imhmin(distance, 0.33)); 
rgb = label2rgb(L, 'jet', [.5 .5 .5]);
figure(), imshow(rgb), title('Watershed transform of distance')


%% With image itself
close all
clc

figure(), imshow(image, [])
image = imfilter(image, fspecial('average', 3));
image = imfilter(image, fspecial('disk', 5));
image = imfilter(image, fspecial('gaussian'));
figure(), imshow(image, [])

L = watershed(imhmin(image, .50));
rgb = label2rgb(L, 'jet', [.5 .5 .5]);
figure(), imshow(rgb), title('Watershed transform of image')


%% Method on MATLAB site (http://www.mathworks.com/matlabcentral/answers/2151)
% http://www.mathworks.com/products/image/examples.html?file=/products/demos/shipping/images/ipexwatershed.html
% Cannot use just nuclei with this method
close all
clc

% Use the gradient magnitude as the segmentation function
hy = fspecial('sobel');
Iy = imfilter(double(image), hy, 'replicate');
Ix = imfilter(double(image), hy', 'replicate');
gradmag = sqrt(Ix.^2 + Iy.^2);
L = watershed(gradmag);
figure(), imshow(gradmag, [min(gradmag(:)), max(gradmag(:)/4)]), title('Gradient magnitude (gradmag)')
figure(), imshow(label2rgb(L)), title('Watershed transform of gradient magnitude (Lrgb)')

% Mark the foreground objects
se = strel('disk', 20);
Io = imopen(image, se);
Ie = imerode(image, se);
Iobr = imreconstruct(Ie, image);
Ioc = imclose(Io, se);
Iobrd = imdilate(Iobr, se);
Iobrcbr = imreconstruct(imcomplement(Iobrd), imcomplement(Iobr));
Iobrcbr = imcomplement(Iobrcbr);
fgm = imregionalmax(Iobrcbr);
image2 = image;
image2(fgm) = 255;
se2 = strel(ones(5,5));
fgm2 = imclose(fgm, se2);
fgm3 = imerode(fgm2, se2);
fgm4 = bwareaopen(fgm3, 20);
image3 = image;
image3(fgm4) = 255;
figure(), imshow(Io), title('Opening (Io)')
figure(), imshow(Ioc), title('Opening-closing (Ioc)')
figure(), imshow(Iobr), title('Opening-by-reconstruction (Iobr)')
figure(), imshow(Iobrcbr), title('Opening-closing by reconstruction (Iobrcbr)')
figure(), imshow(fgm), title('Regional maxima of opening-closing by reconstruction (fgm)')
figure(), imshow(image2), title('Regional maxima superimposed on original image (I2)')
figure(), imshow(image3), title('Modified regional maxima superimposed on original image (fgm4)')

% Compute background markers
bw = im2bw(Iobrcbr, graythresh(Iobrcbr));
distance = bwdist(bw);
DL = watershed(distance);
bgm = DL == 0;
figure(), imshow(bw), title('Thresholded opening-closing by reconstruction (bw)')
figure(), imshow(bgm), title('Watershed ridge lines (bgm)')

% Compute the watershed transform of the segmentation function.
gradmag2 = imimposemin(gradmag, bgm | fgm4);
L = watershed(gradmag2);

% Visualize the result
image4 = image;
image4(imdilate(L == 0, ones(3, 3)) | bgm | fgm4) = 255;
Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
figure(), imshow(image4), title('Markers and object boundaries superimposed on original image (I4)')
figure(), imshow(Lrgb), title('Colored watershed label matrix (Lrgb)')
figure(), imshow(image), hold on
himage = imshow(Lrgb);
set(himage, 'AlphaData', 0.3);
title('Lrgb superimposed transparently on original image')





