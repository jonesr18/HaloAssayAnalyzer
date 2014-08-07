%% Machine Learning
close all
clear all
clc


%% Collect New Apop Data
close all
clc

imArray = Imports.images();
index = 0;
for i = 1:numel(imArray)
    image = imArray{i};
    cells = HaloImage(image);
    cells.crop(1, 454, 1, 708);
    try
        cells.process();
        cells.makeCells('apop', 'bw', 'grad', true);
        figure(1), imshow(cells.image, [])
        for j = 1:numel(cells.aCells)
            
            % Extract data, get classification
            aCell = cells.aCells{j};
            [fourier, magbins, innerbins] = aCell.fft();
            figure(2), imshow(aCell.image, [])
            figure(3), imshow(fourier, [0, 10000])
            answer = input('Cell type? (a/n/h/i)\n', 's');
            
            % Store data - reshaping changes data like so:
            %  [A1 A2 A3
            %   B1 B2 B3    --->   [A1 B1 C1 A2 B2 C2 A3 B3 C3]
            %   C1 C2 C3]
            index = index + 1;
            magSpec(index, :) = reshape(fourier, 1, []);        %#ok<SAGROW>
            allMags(index, :) = reshape(magbins, 1, []);        %#ok<SAGROW>
            innerMags(index, :) = reshape(innerbins, 1, []);    %#ok<SAGROW>
            switch answer
                case 'a'
                    answer = 'apoptotic';
                case 'n'
                    answer = 'necrotic';
                case 'h'
                    answer = 'healthy';
                case 'i'
                    answer = 'ignored';
            end
            answers{index} = answer; %#ok<SAGROW>
        end
    catch
        warning('Threshold not found for image %d', 2 * i)
    end
end

f = datestr(now);
f(f == ':') = ';';
f(f == '-') = ' ';
filename = [f, ' Learning Data.mat'];

save filename magSpec allMags innerMags answers

%% Compile Data
% Used previously - not in use anymore
clc

% -- Dataset non-adj camera H2O2 -- %

allMags = cell(0);
answers = '';
% Denaturing FHA data
for i = [0, 25, 50, 100, 200, 500]
    am = load(sprintf('%d_D1_allMags.mat', i));
    an = load(sprintf('%d_D1_answers.mat', i));
    allMags = [allMags, am.allMags]; %#ok<AGROW>
    answers = [answers, an.answers]; %#ok<AGROW>
end
% Non-denaturing FHA data
for i = [50, 200, 500]
    am = load(sprintf('%d_N1_allMags.mat', i));
    an = load(sprintf('%d_N1_answers.mat', i));
    allMags = [allMags, am.allMags]; %#ok<AGROW>
    answers = [answers, an.answers]; %#ok<AGROW>
end

save non_adj_H2O2_allMags.mat allMags
save non_adj_H2O2_answers.mat answers

% -- Dataset adj camera H2O2 -- %

allMags = cell(0);
answers = '';
for i = [0, 25, 50, 100, 500]
    am = load(sprintf('%d_D1_newAM.mat', i));
    an = load(sprintf('%d_D1_newAn.mat', i));
    allMags = [allMags, am.allMags]; %#ok<AGROW>
    answers = [answers, an.answers]; %#ok<AGROW>
end

save adj_H2O2_allMags.mat allMags
save adj_H2O2_answers.mat answers


%% Make Data for Saved Files
close all
clc

dataset = 'non_adj';
switch dataset
    case 'non_adj'
        am = load('non_adj_H2O2_allMags.mat');
        allMags = am.allMags;
        an = load('non_adj_H2O2_answers.mat');
        answers = an.answers;
    case 'adj'
        am = load('adj_H2O2_allMags.mat');
        allMags = am.allMags;
        an = load('adj_H2O2_answers.mat');
        answers = an.answers;
end

numCells = length(answers);

% Compile magnitude specturm data
index = 0;
compMags = zeros(10, 10, numCells);
for mags = allMags
    mag = mags{:};
    if ~isempty(mag)
        for i = 1:size(mag, 3)
            index = index + 1;
            compMags(:, :, index) = mag(:, :, i);
        end
    end
end

% Reshape data
shifted = shiftdim(compMags, 2);
data = reshape(shifted, size(compMags, 3), []);
a = answers == 'a';
n = answers == 'n';
c = answers == 'h' | answers == 'c';
d = answers == 'd';
i = answers == 'i' | answers == 'b';


%% Principal Component Analysis
close all
clc

% Unwieghted
[coeff, score, ~, ~, vars] = pca(data);
disp(vars(1:5))
figure()
plot3(score(a, 1), score(a, 2), score(a, 3), 'ro',...
      score(n, 1), score(n, 2), score(n, 3), 'go',...
      score(c, 1), score(c, 2), score(c, 3), 'bo',...
      score(d, 1), score(d, 2), score(d, 3), 'mo')
xlabel(sprintf('C1: %.2f%%', vars(1)), 'fontsize', 14)
ylabel(sprintf('C2: %.2f%%', vars(2)), 'fontsize', 14)
zlabel(sprintf('C3: %.2f%%', vars(3)), 'fontsize', 14)
title('Bulk PCA of 2D Frequency Spectrum', 'fontsize', 14)

% Weighted
weights = 1./var(data);
[~, wScore, ~, ~, wVars] = pca(data, 'VariableWeights', weights);
disp(wVars(1:5))
figure()
plot3(wScore(a, 1), wScore(a, 2), wScore(a, 3), 'ro',...
      wScore(n, 1), wScore(n, 2), wScore(n, 3), 'go',...
      wScore(c, 1), wScore(c, 2), wScore(c, 3), 'bo',...
      wScore(d, 1), wScore(d, 2), wScore(d, 3), 'mo')
xlabel(sprintf('C1: %.2f%%', wVars(1)), 'fontsize', 14)
ylabel(sprintf('C2: %.2f%%', wVars(2)), 'fontsize', 14)
zlabel(sprintf('C3: %.2f%%', wVars(3)), 'fontsize', 14)
title('Weighted Bulk PCA of 2D Frequency Spectrum', 'fontsize', 14)

%{
figure()
spot = [1, 2, 6, 7];
for i = 1:4
    subplot(2, 5, spot(i)), imshow(reshape(coeff(:, i), 10, 10), [-1, 1])
    title(sprintf('C%d: %.2f%%', i, vars(i)))
end
colormap(gray)
%}


%% Partial Correlation Coefficients
% close all
clc

isapop = boolean(answers == 'a');
isnecr = boolean(answers == 'n');
iscont = boolean(answers == 'c');

[coeffApop, pvalsApop] = partialcorr([data, isapop']);
[coeffNecr, pvalsNecr] = partialcorr([data, isnecr']);
[coeffCont, pvalsCont] = partialcorr([data, iscont']);

% Reshape correlation and pvals into image positions
coeffApop = reshape(coeffApop(end, 1 : end - 1), 10, 10);
coeffNecr = reshape(coeffNecr(end, 1 : end - 1), 10, 10);
coeffCont = reshape(coeffCont(end, 1 : end - 1), 10, 10);
pvalsApop = reshape(pvalsApop(end, 1 : end - 1), 10, 10);
pvalsNecr = reshape(pvalsNecr(end, 1 : end - 1), 10, 10);
pvalsCont = reshape(pvalsCont(end, 1 : end - 1), 10, 10);

subplot(2, 5, 3), imshow(coeffApop, [-1, 1]), title('Apop Corr')
subplot(2, 5, 4), imshow(coeffNecr, [-1, 1]), title('Necr Corr')
subplot(2, 5, 5), imshow(coeffCont, [-1, 1]), title('Cont Corr')
subplot(2, 5, 8), imshow(.05 ./ pvalsApop), title('Apop p-vals')
subplot(2, 5, 9), imshow(.05 ./ pvalsNecr), title('Necr p-vals')
subplot(2, 5, 10), imshow(.05 ./ pvalsCont), title('Cont p-vals')
colormap(gray)


%% Cell Type Classification
close all
clc

color = 'rgbk';
symb = 'x+^*';

% Split the sample set
indexes = randperm(size(data, 1));
data1 = data(indexes(1 : round(end / 2)), :);
data2 = data(indexes(round(end / 2) + 1 : end), :);

% Set new type indexes
type1 = answers(indexes(1 : round(end / 2)));
type2 = answers(indexes(round(end / 2) + 1 : end));
types1 = {find(type1 == 'a'), find(type1 == 'n'), find(type1 == 'c'), find(type1 == 'd')};
types2 = {find(type2 == 'a'), find(type2 == 'n'), find(type2 == 'c'), find(type2 == 'd')};


%% Unsupervised K-means
close all
clc

% Make classifier with first set
[idx, cn] = kmeans(data1, 3, 'replicates', 3);

% Check classifier against seed data
[~, score] = pca(data1);
clstr = {find(idx == 1), find(idx == 2), find(idx == 3)};
figure(), hold on
for i = 1:numel(clstr)
    for j = 1:numel(types1)
        k = ismember(clstr{i}, types1{j}); % In cluster and of correct type
        plot(score(clstr{i}(k), 1), score(clstr{i}(k), 2), [color(i), symb(j)])
    end
end
xlabel('Component 1', 'fontsize', 14)
ylabel('Component 2', 'fontsize', 14)
title('Classification of Seed Data with 3 K-means Clusters', 'fontsize', 14), hold off

% Check classifier against withheld data
tic
dists = dist([data2', cn']);
dists = dists(end - 2 : end, 1 : end - 3);
[~, idx] = min(dists);
toc
[~, score] = pca(data2);
clstr = {find(idx == 1), find(idx == 2), find(idx == 3)};
figure(), hold on
for i = 1:numel(clstr)
    for j = 1:numel(types2)
        k = ismember(clstr{i}, types2{j}); % In cluster and of correct type
        plot(score(clstr{i}(k), 1), score(clstr{i}(k), 2), [color(i), symb(j)])
    end
end
xlabel('Component 1', 'fontsize', 14)
ylabel('Component 2', 'fontsize', 14)
title('Classification of Withheld Data with 3 K-means Clusters', 'fontsize', 14)


%% Supervised K-means
close all
clc

% Make classifier with first set
[~, cn_a] = kmeans(data1(type1 == 'a', :), 3, 'replicates', 3);
[~, cn_n] = kmeans(data1(type1 == 'n', :), 3, 'replicates', 3);
[~, cn_c] = kmeans(data1(type1 == 'c', :), 3, 'replicates', 3);

% Check classifier against seed data
dists = dist([data1', cn_a', cn_n', cn_c']);
dists = dists(end - 8 : end, 1 : end - 9);
[~, idx] = min(dists);
[~, score] = pca(data1);
clstr = {find(idx <= 3), find(idx > 3 & idx <= 6), find(idx > 6)};
total1 = 0;
figure(), hold on
for i = 1:numel(clstr)
    for j = 1:numel(types1)
        k = ismember(clstr{i}, types1{j}); % In cluster and of correct type
        total1 = total1 + sum(k);
        plot(score(clstr{i}(k), 1), score(clstr{i}(k), 2), [color(i), symb(j)])
    end
end
xlabel('Component 1', 'fontsize', 14)
ylabel('Component 2', 'fontsize', 14)
title('Classification of Seed Data with 4x3 K-means Clusters', 'fontsize', 14), hold off

% Determine fraction correct
found = sum([ismember(clstr{1}, types1{1}),...
             ismember(clstr{2}, types1{2}),...
             ismember(clstr{3}, types1{3})]);
frac = found / total1 * 100

% Check classifier against withheld data
tic
dists = dist([data2', cn_a', cn_n', cn_c']);
dists = dists(end - 8 : end, 1 : end - 9);
[~, idx] = min(dists);
toc
[~, score] = pca(data2);
clstr = {find(idx <= 3), find(idx > 3 & idx <= 6), find(idx > 6)};
total2 = 0;
figure(), hold on
for i = 1:numel(clstr)
    for j = 1:numel(types2)
        k = ismember(clstr{i}, types2{j}); % In cluster and of correct type
        total2 = total2 + sum(k);
        plot(score(clstr{i}(k), 1), score(clstr{i}(k), 2), [color(i), symb(j)])
    end
end
xlabel('Component 1', 'fontsize', 14)
ylabel('Component 2', 'fontsize', 14)
title('Classification of Withheld Data with 4x3 K-means Clusters', 'fontsize', 14)

% Determine fraction correct
found = sum([ismember(clstr{1}, types2{1}),...
             ismember(clstr{2}, types2{2}),...
             ismember(clstr{3}, types2{3})]);
frac = found / total2 * 100


%% Classification Tree
close all
clc

% Build classifier
tree = ClassificationTree.fit(data1, type1');
% tree.view('mode', 'graph')

% Check classifier against seed data
labels = tree.predict(data1);
[~, score] = pca(data1);
clstr = {find(labels == 'a'), find(labels == 'n'), find(labels == 'c')};
total1 = 0;
figure(), hold on
for i = 1:numel(clstr)
    for j = 1:numel(types1)
        k = ismember(clstr{i}, types1{j}); % In cluster and of correct type
        total1 = total1 + sum(k);
        plot(score(clstr{i}(k), 1), score(clstr{i}(k), 2), [color(i), symb(j)])
    end
end
xlabel('Component 1', 'fontsize', 14)
ylabel('Component 2', 'fontsize', 14)
title('Classification of Seed Data with Tree', 'fontsize', 14), hold off

% Determine fraction correct
found = sum([ismember(clstr{1}, types1{1});...
             ismember(clstr{2}, types1{2});...
             ismember(clstr{3}, types1{3})]);
frac = found / total1 * 100

% Check classifier against withheld data
tic
labels = tree.predict(data2);
toc
[~, score] = pca(data2);
clstr = {find(labels == 'a'), find(labels == 'n'), find(labels == 'c')};
total2 = 0;
figure(), hold on
for i = 1:numel(clstr)
    for j = 1:numel(types2)
        k = ismember(clstr{i}, types2{j}); % In cluster and of correct type
        total2 = total2 + sum(k);
        plot(score(clstr{i}(k), 1), score(clstr{i}(k), 2), [color(i), symb(j)])
    end
end
xlabel('Component 1', 'fontsize', 14)
ylabel('Component 2', 'fontsize', 14)
title('Classification of Withheld Data with Tree', 'fontsize', 14), hold off

% Determine fraction correct
found = sum([ismember(clstr{1}, types2{1});...
             ismember(clstr{2}, types2{2});...
             ismember(clstr{3}, types2{3})]);
frac = found / total2 * 100


%% K-Nearest Neighbors
close all
clc

% Build classifier
knn = ClassificationKNN.fit(data1, type1', 'DistanceWeight', 'squaredinverse');

% Check classifier against seed data
labels = knn.predict(data1);
[~, score] = pca(data1);
clstr = {find(labels == 'a'), find(labels == 'n'), find(labels == 'c')};
total1 = 0;
figure(), hold on
for i = 1:numel(clstr)
    for j = 1:numel(types1)
        k = ismember(clstr{i}, types1{j}); % In cluster and of correct type
        total1 = total1 + sum(k);
        plot(score(clstr{i}(k), 1), score(clstr{i}(k), 2), [color(i), symb(j)])
    end
end
xlabel('Component 1', 'fontsize', 14)
ylabel('Component 2', 'fontsize', 14)
title('Classification of Seed Data by K-NN', 'fontsize', 14), hold off

% Determine fraction correct
found = sum([ismember(clstr{1}, types1{1});...
             ismember(clstr{2}, types1{2});...
             ismember(clstr{3}, types1{3})]);
frac = found / total1 * 100

% Check classifier against withheld data
tic
labels = knn.predict(data2);
toc
[~, score] = pca(data2);
clstr = {find(labels == 'a'), find(labels == 'n'), find(labels == 'c')};
total2 = 0;
figure(), hold on
for i = 1:numel(clstr)
    for j = 1:numel(types2)
        k = ismember(clstr{i}, types2{j}); % In cluster and of correct type
        total2 = total2 + sum(k);
        plot(score(clstr{i}(k), 1), score(clstr{i}(k), 2), [color(i), symb(j)])
    end
end
xlabel('Component 1', 'fontsize', 14)
ylabel('Component 2', 'fontsize', 14)
title('Classification of Withheld Data by K-NN', 'fontsize', 14), hold off

% Determine fraction correct
found = sum([ismember(clstr{1}, types2{1});...
             ismember(clstr{2}, types2{2});...
             ismember(clstr{3}, types2{3})]);
frac = found / total2 * 100


%% Expectation Maximization - Gaussian Mixtrure
close all
clc

% Build classifier using first 5 principal components (using just data1 didn't work)
[~, score] = pca(data1);
emgm = gmdistribution.fit(score(:, 1:5), 3, 'replicates', 3);

% Check classifier against seed data
idx = emgm.cluster(score(:, 1:5));
clstr = {find(idx == 1), find(idx == 2), find(idx == 3)};
figure(), hold on
for i = 1:numel(clstr)
    for j = 1:numel(types1)
        k = ismember(clstr{i}, types1{j}); % In cluster and of correct type
        plot(score(clstr{i}(k), 1), score(clstr{i}(k), 2), [color(i), symb(j)])
    end
end
xlabel('Component 1', 'fontsize', 14)
ylabel('Component 2', 'fontsize', 14)
title('Classification of Seed Data by EM-GM', 'fontsize', 14), hold off

% Check classifier against withheld data
[~, score] = pca(data2);
tic
idx = emgm.cluster(score(:, 1:5));
toc
clstr = {find(idx == 1), find(idx == 2), find(idx == 3)};
figure(), hold on
for i = 1:numel(clstr)
    for j = 1:numel(types2)
        k = ismember(clstr{i}, types2{j}); % In cluster and of correct type
        plot(score(clstr{i}(k), 1), score(clstr{i}(k), 2), [color(i), symb(j)])
    end
end
xlabel('Component 1', 'fontsize', 14)
ylabel('Component 2', 'fontsize', 14)
title('Classification of Withheld Data by EM-GM', 'fontsize', 14), hold off



