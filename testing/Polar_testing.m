%% Polar image testing
close all
clear all
clc

warning('off', 'images:initSize:adjustingMag');
image = imread('C:\Users\Ross\Google Drive\Capstone Images\Binh\Array Image-1.tif');
image = imfilter(rgb2gray(image), fspecial('average'));
cells = HaloImage(image);
cells.process();
cells.makeCells('apop')


%% cart2pol
close all
clc

im = cells.aCells{1}.image;
thresh = imquantize(im, multithresh(im, 2));
nuc = boolean(thresh > 2);
[row, col] = find(nuc);
cn = [mean(row), mean(col)];
[H, W] = size(im);
[x, y] = meshgrid(1:W, 1:H);

[theta, rho] = cart2pol(x - cn(1), y - cn(2));

% figure()
% subplot(1, 2, 1), imshow(uint8(10 * (theta + pi)), [])
% subplot(1, 2, 2), imshow(uint8(rho), [])

figure(), surf(rho, theta, double(im)), shading interp


%% radialDist
% Reconstruct missing parts of image using regression of pixel intensity at radial distance r. Not 
% necessary with updated segmentation scheme.
close all
clc

for n = 1:numel(cells.aCells)
    c = cells.aCells{n};
    r = regionprops(c.image, 'Centroid');
    
    % Find r and th for each point > 0
    [H, W] = size(c.image);
    stats = zeros(length(c.image > 0), 2);
    index = 1;
    for i = 1:H
        for j = 1:W
            if c.image(i, j) > 0
                dist = sqrt((i - r(1).Centroid(1))^2 + (j - r(1).Centroid(2))^2);
                stats(index, 1:2) = [dist, c.image(i, j)];
                index = index + 1;
            end
        end
    end

    % Find regression for estimation of points to fill.
    [~, slope, intercept] = regression(stats(:, 1)', stats(:, 2)');

    % Predict points using regression
    for i = 1:H
        for j = 1:W
            if c.image(i, j) == 0
                dist = sqrt((i - r(1).Centroid(1))^2 + (j - r(1).Centroid(2))^2);
                c.image(i, j) = slope * dist + intercept;
            end
        end
    end

    c.image = imfilter(c.image, fspecial('average', 5));

    ttl = sprintf('Fourier Transform of image %d', i);
    fourier = c.fft();
    figure(1), imshow(c.image, []), title(sprintf('Image %d', i))
    figure(2), imshow(fourier, [0, 1000]), title(ttl)
    pause()
end