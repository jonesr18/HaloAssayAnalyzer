%% Class testing (HaloImage & CellRegion)
close all
clear all
clc

warning('off', 'images:initSize:adjustingMag');

% image = imread('C:\Users\Ross\Google Drive\Capstone Images\Binh\Array Image-1.tif');
% image = imread('C:\Users\Ross\Google Drive\Capstone Images\Binh\X-ray for A172.tif');
% image = imread('C:\Users\Ross\Google Drive\Capstone Images\Binh\nanoparticle for HeLa.png');
% image = imread('C:\Users\Ross\Google Drive\Capstone Images\Binh\metal ion for HeLa.png');
% image = imread('C:\Users\Ross\Google Drive\Capstone Images\Binh\figure pic.png');
% image = imread('C:\Users\Ross\Google Drive\Capstone Images\87 - 12.26.13 - FHA T5 (EB)\Images\Image-51.tif');
image = imread('C:\Users\Ross\Google Drive\Capstone Images\97 - 05.11.14 200 D1 (EB)\Images\Image-5.tif');

image = image(1:454, 1:708); % Crop our microscope images

if size(image, 3) == 4
    image = cmyk2rgb(image);
end
if size(image, 3) == 3
    image = rgb2gray(image);
end
image1 = imfilter(image, fspecial('gaussian'));
image2 = imfilter(image1, fspecial('average'));   


%% HaloImage Damage
close all
clc

if exist('cells', 'var')
    cells.delete()
    clear cells
end
cells = HaloImage(image2);
cells.process();
tic
cells.makeCells('damage', 'bw', 'grad', true);
toc

cells.calcDamage('all');
names = fieldnames(cells.damage);
for f = names'
    disp(f{:})
    disp(cells.damage.(f{:}))
end

% pause()
cells.plot('gradmag', 'binary', 'labels')
figure(), imshow(cells.image .* uint8(cells.binary), [])


%% HaloImage Apoptosis
close all
clc

if exist('cells', 'var')
    cells.delete()
    clear cells
end
cells = HaloImage(image);
% cells.setThresh(10, 130);
cells.process();
cells.makeCells('apoptosis', 'bw', 'thresh', true)

% for i = 1:numel(cells.aCells)
%     ttl = sprintf('Fourier Transform of image %d', i);
%     fourier = cells.aCells{i}.fft();
%     figure(1), imshow(cells.aCells{i}.image, []), title(sprintf('Image %d', i))
%     figure(2), imshow(fourier, [0, 10000]), title(ttl)
%     pause()
% end

am = load('non_adj_H2O2_allMags.mat');
an = load('non_adj_H2O2_answers.mat');
cells.classify(am.allMags, an.answers, 'kmeans');
cells.plotTypes()


%% GUI Testing
close all
clc

if exist('gui', 'var')
    try
        delete(gui)
    catch
    end
    clear gui
end
gui = Initialization;


%% Imports/Exports
close all
clc


out = Exports();
out.writeCSV('random', cell(1));



