%% Testing for different methods of cluster identification
clear all
clc

% Image numbers
apop = [1 13 19 35 39];
necr = [9 13 21 49 57];
cont = [1 13 21 45 57];

% Whole Cell or Nucleus choice
test = 'whole';
% test = 'nucleus';

% FFT data directories
switch test
    case {'whole'}
        apopDir = 'C:\Users\Ross\Google Drive\Capstone Images\87 - 10.7.13 - FHA 2 (EB)\Data\Whole Cell FFT\';
        necrDir = 'C:\Users\Ross\Google Drive\Capstone Images\87 - 12.26.13 - FHA T5 (EB)\Data\Whole Cell FFT\';
        contDir = 'C:\Users\Ross\Google Drive\Capstone Images\87 - 11.21.13 - FHA C1 (EB)\Data\Whole Cell FFT\';
    case {'nucleus'}
        apopDir = 'C:\Users\Ross\Google Drive\Capstone Images\87 - 10.7.13 - FHA 2 (EB)\Data\Nucleus FFT\';
        necrDir = 'C:\Users\Ross\Google Drive\Capstone Images\87 - 12.26.13 - FHA T5 (EB)\Data\Nucleus FFT\';
        contDir = 'C:\Users\Ross\Google Drive\Capstone Images\87 - 11.21.13 - FHA C1 (EB)\Data\Nucleus FFT\';
end

% Initialize data matricies
apopDataMatrix = zeros(128, 128, 25);
necrDataMatrix = zeros(128, 128, 25);
contDataMatrix = zeros(128, 128, 25);

% Extract and organize data from FFT data
index = 0;
for i = 1:5
    for j = 1:5
        % Load data
        apopFile = sprintf('%s FFT %d of Image-%d.txt', test, i, apop(j));
        necrFile = sprintf('%s FFT %d of Image-%d.txt', test, i, necr(j));
        contFile = sprintf('%s FFT %d of Image-%d.txt', test, i, cont(j));
        apopData = load(strcat(apopDir, apopFile));
        necrData = load(strcat(necrDir, necrFile));
        contData = load(strcat(contDir, contFile));
        
        % Extract data
        apopMatrix = ones(max(apopData(:, 1)) + 1);
        necrMatrix = ones(max(necrData(:, 1)) + 1);
        contMatrix = ones(max(contData(:, 1)) + 1);
        for k = 1:length(apopData)
            apopMatrix(apopData(k, 1) + 1, apopData(k, 2) + 1) = apopData(k, 3);
        end
        for k = 1:length(necrData)
            necrMatrix(necrData(k, 1) + 1, necrData(k, 2) + 1) = necrData(k, 3);
        end
        for k = 1:length(contData)
            contMatrix(contData(k, 1) + 1, contData(k, 2) + 1) = contData(k, 3);
        end
        
        % Fix matrices to be correct size
        apopMatrix = imresize(apopMatrix, [128, 128]);
        necrMatrix = imresize(necrMatrix, [128, 128]);
        contMatrix = imresize(contMatrix, [128, 128]);
        
        % Store data
        index = index + 1;
        apopDataMatrix(:, :, index) = apopMatrix;
        necrDataMatrix(:, :, index) = necrMatrix;
        contDataMatrix(:, :, index) = contMatrix;
    end
end


%% Bulk analysis

% Average squares of size 8x8 and create data matrix for PCA
bData = zeros(75, 256);
for i = 1:16
    for j = 1:16
       for k = 1:25
           rows = 8 * (i - 1) + 1 : 8 * i;
           cols = 8 * (j - 1) + 1 : 8 * j;
           index = 16 * (i - 1) + j;
           bData(k, index) = mean(mean(mean(apopDataMatrix(rows, cols, k))));
           bData(k + 25, index) = mean(mean(mean(necrDataMatrix(rows, cols, k))));
           bData(k + 50, index) = mean(mean(mean(contDataMatrix(rows, cols, k))));
       end
    end
end


%% Individual analysis
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

% Sort indexes from lowest to highest
indexes = sort([caTopRQ; naTopRQ; ncTopRQ]);
n = length(indexes);

% Create vectors of data to cluster with k-means, each vector is for one axis
data  = zeros(75, n);

for i = 1:n
    % Find x and y coordinates from indexes
    y = mod(indexes(i), 128);
    x = (indexes(i) - y) / 128;
    
    for j = 1:25
        data(j, i) = apopDataMatrix(y, x, j);
        data(j + 25, i) = necrDataMatrix(y, x, j);
        data(j + 50, i) = contDataMatrix(y, x, j);
    end
end


%% K-means
close all
clc

[IDX, centroids] = kmeans(data, 3, 'replicates', 10);

% After many replicates, these are always the same. 
apopDef = find(IDX(1:25) ~= mode(IDX(1:25)))     % none
necrDef = find(IDX(26:50) ~= mode(IDX(26:50)))   % 9, 20 -- look like control, perhaps are?
contDef = find(IDX(51:end) ~= mode(IDX(51:end))) % 22 -- looks like apoptosis?

% Bulk data -- gets them ALL correct
[IDX, centroids] = kmeans(bData, 3, 'replicates', 3);

apopDef = find(IDX(1:25) ~= mode(IDX(1:25))) % 
necrDef = find(IDX(26:50) ~= mode(IDX(26:50))) % 
contDef = find(IDX(51:end) ~= mode(IDX(51:end))) % 


%% Principal Components
% Works great for visualization!
close all
clc

% Normal -- fairly comparable, I think the weighted is slightly better
[~, score, ~, ~, vars] = pca(data);
disp(vars(1:5))
figure(), plot(score(1:25, 1), score(1:25, 2), 'or', score(26:50, 1), score(26:50, 2), 'ob', score(51:end, 1), score(51:end, 2), 'og')
title('normal')

% Weighted
weights = 1./var(data);
[~, wScore, ~, ~, wVars] = pca(data, 'VariableWeights', weights);
disp(wVars(1:5))
figure(), plot(wScore(1:25, 1), wScore(1:25, 2), 'or', wScore(26:50, 1), wScore(26:50, 2), 'ob', wScore(51:end, 1), wScore(51:end, 2), 'og')
title('PCA Analysis of 2 dimensions'), legend('Apoptotic', 'Necrotic', 'Control')
xlabel('Component 1'), ylabel('Component 2')


%% PCA of bulk analysis

% Normal
[~, score, ~, ~, vars] = pca(bData);
disp(vars(1:5))
figure(), plot(score(1:25, 1), score(1:25, 2), 'or', score(26:50, 1), score(26:50, 2), 'og', score(51:end, 1), score(51:end, 2), 'ob')
xlabel(sprintf('C1: %.2f%%', vars(1)), 'fontsize', 14)
ylabel(sprintf('C2: %.2f%%', vars(2)), 'fontsize', 14)
title('Bulk PCA of 2D Frequency Spectrum', 'fontsize', 14)

% Weighted -- not as good actually!
weights = 1./var(bData);
[~, wScore, ~, ~, wVars] = pca(bData, 'VariableWeights', weights);
disp(wVars(1:5))
figure(), plot(wScore(1:25, 1), wScore(1:25, 2), 'or', wScore(26:50, 1), wScore(26:50, 2), 'og', wScore(51:end, 1), wScore(51:end, 2), 'ob')
xlabel(sprintf('C1: %.2f%%', wVars(1)), 'fontsize', 14)
ylabel(sprintf('C2: %.2f%%', wVars(2)), 'fontsize', 14)
title('Weighted Bulk PCA of 2D Frequency Spectrum', 'fontsize', 14)



%% Factor Analysis
% Doesn't work with this data
close all
clc

[lambda, PSI, T, STATS, F] = factoran(data, 3);


%% Heirarchical
% Doesn't give much useful info
close all
clc

distance = pdist(data, 'euclidean');
clustTree = linkage(distance,'average');
cophCorr = cophenet(clustTree, distance)
figure(), dendrogram(clustTree, 0)

