%% FHA
clear all
clc

% Image numbers
apop = [1 13 19 35 39];
necr = [9 13 21 49 57];
cont = [1 13 21 45 57];

% FFT data directories
apopDir = 'C:\Users\Ross\Google Drive\Capstone Images\87 - 10.7.13 - FHA 2 (EB)\Data\';
necrDir = 'C:\Users\Ross\Google Drive\Capstone Images\87 - 12.26.13 - FHA T5 (EB)\Data\';
contDir = 'C:\Users\Ross\Google Drive\Capstone Images\87 - 11.21.13 - FHA C1 (EB)\Data\';

% Initialize data matricies
apopDataMatrix = zeros(128, 128, 25);
necrDataMatrix = zeros(128, 128, 25);
contDataMatrix = zeros(128, 128, 25);
index = 0;

% Extract and organize data from FFT data
for i = 1:5
    for j = 1:5
        apopFile = sprintf('FFT %d of Image-%d.txt', i, apop(j));
        necrFile = sprintf('FFT %d of Image-%d.txt', i, necr(j));
        contFile = sprintf('FFT %d of Image-%d.txt', i, cont(j));
        apopData = load(strcat(apopDir, apopFile));
        necrData = load(strcat(necrDir, necrFile));
        contData = load(strcat(contDir, contFile));
        apopMatrix = ones(128);
        necrMatrix = ones(128);
        contMatrix = ones(128);
        for k = 1:length(necrData)
            x = necrData(k, 1) + 1;
            y = necrData(k, 2) + 1;
            necrMatrix(2 * x - 1 : 2 * x, 2 * y - 1 : 2 * y) = necrData(k, 3);
        end
        for k = 1:length(contData)
            x = contData(k, 1) + 1;
            y = contData(k, 2) + 1;
            contMatrix(2 * x - 1 : 2 * x, 2 * y - 1 : 2 * y) = contData(k, 3);
        end
        for k = 1:length(apopData)
            apopMatrix(apopData(k, 1) + 1, apopData(k, 2) + 1) = apopData(k, 3);
        end
        index = index + 1;
        apopDataMatrix(:, :, index) = mean(apopMatrix, 3);
        necrDataMatrix(:, :, index) = necrMatrix;
        contDataMatrix(:, :, index) = contMatrix;
    end
end

% Initialize ttest result arrays
pvals = zeros(128, 128, 3);
signif = zeros(128, 128, 3);

% Find significant differences between FFT patterns
for i = 1:128
    for j = 1:128
        [~, pvals(i, j, 1)] = ttest2(apopDataMatrix(i, j, :), contDataMatrix(i, j, :));
        [~, pvals(i, j, 2)] = ttest2(apopDataMatrix(i, j, :), necrDataMatrix(i, j, :));
        [~, pvals(i, j, 3)] = ttest2(necrDataMatrix(i, j, :), contDataMatrix(i, j, :));
        for k = 1:3
            if pvals(i, j, k) < .001
                signif(i, j, k) = 3;
            elseif pvals(i, j, k) < .01
                signif(i, j, k) = 2;
            elseif pvals(i, j, k) < .05
                signif(i, j, k) = 1;
            else
                signif(i, j, k) = 0;
            end
        end
    end
end

% Average all FFTs for each cell type
apopMean = mean(apopDataMatrix, 3);
necrMean = mean(necrDataMatrix, 3);
contMean = mean(contDataMatrix, 3);

% Find relative quantities (RQs) for comparison between each cell type FFT
caRQ = contMean ./ apopMean;
naRQ = necrMean ./ apopMean;
ncRQ = necrMean ./ contMean;

% Identify the indexes of RQ values which are at least 85% of the max -- better than 80% and 90% 
caTopRQ = find(caRQ > 0.85 * max(caRQ(:)));
naTopRQ = find(naRQ > 0.85 * max(naRQ(:)));
ncTopRQ = find(ncRQ > 0.85 * max(ncRQ(:)));

% Create binary image where top RQs are 1 and else are 0
caBinRQ = zeros(128); caBinRQ(caTopRQ) = 1;
naBinRQ = zeros(128); naBinRQ(naTopRQ) = 1;
ncBinRQ = zeros(128); ncBinRQ(ncTopRQ) = 1;


%% Plots
close all

% Surface plots
figure()
subplot(1, 3, 1), surf(apopMean), title('Apoptosis')
subplot(1, 3, 2), surf(necrMean), title('Necrosis')
subplot(1, 3, 3), surf(contMean), title('Control')

figure() % P-values
subplot(1, 3, 1), imshow(label2rgb(signif(:, :, 1))), title('CA Signif')
subplot(1, 3, 2), imshow(label2rgb(signif(:, :, 2))), title('NA Signif')
subplot(1, 3, 3), imshow(label2rgb(signif(:, :, 3))), title('NC Signif')

figure() % Relative quantitites
subplot(1, 3, 1), surf(caRQ), title('RQ Cont / Apop')
subplot(1, 3, 2), surf(naRQ), title('RQ Necr / Apop')
subplot(1, 3, 3), surf(ncRQ), title('RQ Necr / Cont')

figure() % Top relative quantitites
subplot(1, 3, 1), imshow(caBinRQ), title('Top Cont / Apop RQ')
subplot(1, 3, 2), imshow(naBinRQ), title('Top Necr / Apop RQ')
subplot(1, 3, 3), imshow(ncBinRQ), title('Top Necr / Cont RQ')


%% K-means of FFT data
close all
clc

indexes = sort([caTopRQ; naTopRQ; ncTopRQ]);
n = length(indexes);

% Create vectors of data to cluster with k-means, each vector is for one axis
data  = zeros(25, n);
aData = zeros(25, n);
nData = zeros(25, n);
cData = zeros(25, n);

for i = 1:n
    % Find x and y coordinates from indexes
    y = mod(indexes(i), 128);
    x = (indexes(i) - y) / 128;
    
    for j = 1:25
        data(j * 3 - 2 : j * 3, i) = [apopDataMatrix(y, x, j), necrDataMatrix(y, x, j), contDataMatrix(y, x, j)];
        aData(i, n) = apopDataMatrix(y, x, j);
        nData(i, n) = necrDataMatrix(y, x, j);
        cData(i, n) = contDataMatrix(y, x, j);
    end
end

% K-means
[IDX, centroids] = kmeans(data, 3, 'replicates', 10);

% After many replicates, these are always the same. 
apopDef = find(IDX(1:3:end) ~= mode(IDX(1:3:end))) % none
necrDef = find(IDX(2:3:end) ~= mode(IDX(2:3:end))) % 9, 20 -- look like control, perhaps are?
contDef = find(IDX(3:3:end) ~= mode(IDX(3:3:end))) % 22 -- looks like apoptosis?


%% Matlab FFTs
close all
clear all
clc

% Import image
apop = imread('C:\Users\Ross\Google Drive\Capstone Images\87 - 10.7.13 - FHA 2 (EB)\Images\Image-31.tif');
necr = imread('C:\Users\Ross\Google Drive\Capstone Images\87 - 12.26.13 - FHA T5 (EB)\Images\Image-31.tif');
cont = imread('C:\Users\Ross\Google Drive\Capstone Images\87 - 11.21.13 - FHA C1 (EB)\Images\Image-31.tif');

% Convert to grayscale
apop = rgb2gray(apop);
necr = rgb2gray(necr);
cont = rgb2gray(cont);

% Plot grayscale images
figure(), imshow(apop), title('Apop Image')
figure(), imshow(necr), title('Necr Image')
figure(), imshow(cont), title(' Cont Image')

% Plot magnitudes of image FFTs
figure(), imshow(abs(fftshift(fft2(double(apop)))), [0, 10000]), title('Apop Magnitude')
figure(), imshow(abs(fftshift(fft2(double(necr)))), [0, 10000]), title('Necr Magnitude')
figure(), imshow(abs(fftshift(fft2(double(cont)))), [0, 10000]), title('Cont Magnitude')


%% EM_GM Testing
close all
clear all
clc

X = zeros(600, 2); 
X(1:200, :) = normrnd(0, 1, 200, 2); 
X(201:400, :) = normrnd(0, 2, 200, 2); 
X(401:600, :) = normrnd(0, 3, 200, 2); 
[W, M, V, L] = EM_GM(X, 3, [], [], 1, []);

