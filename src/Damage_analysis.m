%% Comet Assay Data Analysis
close all
clear all
clc

% Import data -- some things are changed for each experiment
doses = [0, 25, 50, 100]; % Set doses used
numGroups = length(doses);
data = cell(1, numGroups);
stats = cell(1, numGroups);
in = Imports();
for i = 1:numGroups
    % Set filename regexp
    filename1 = sprintf('07.30.14 %d rads C1 (EB)', doses(i));
    filename2 = sprintf('07.30.14 %d rads C2 (EB)', doses(i));
    filename3 = sprintf('07.31.14 %d rads C1 (EB)', doses(i));
    filename4 = sprintf('07.31.14 %d rads C2 (EB)', doses(i));
    if i == 1
        [data1, ~, labels, statLabels] = in.cometData(filename1);
    else
        [data1, ~] = in.cometData(filename1);
    end
    [data2, ~] = in.cometData(filename2);
    [data3, ~] = in.cometData(filename3);
    [data4, ~] = in.cometData(filename4);
    combData = [data1; data2; data3; data4];
    data{i} = combData;
    stats{i} = [min(combData); max(combData); mean(combData); std(combData); median(combData)];
end
tests = {'Comet extent', 'Comet total intensity', 'Comet total area', 'Tail DNA[%]', 'Tail extent',...
         'Tail distributed moment (DNA migration)', 'Tail extent moment', 'Tail area (Singh)',...
         'Tail (Olive) moment /100'};

% Format data into struct
cometData = struct();
for i = 1:numGroups
    for j = 1:numel(labels)
        % Remove illegal struct name characters
        L = labels{j};
        illegal = ' []()/';
        for il = illegal
            L(L == il) = '_';
        end
        L(L == '%') = 'P'; % for Tail DNA[%]
        
        % Build struct
        D = sprintf('d%d', doses(i));
        cometData.(D).(L).data = data{i}(:, j);
        for k = 1:numel(statLabels)
            cometData.(D).(L).(statLabels{k}) = stats{i}(k, j);
        end
    end
end


%% ACA Histograms
close all
clc

warning('off', 'stats:lillietest:OutOfRangePLow');
warning('off', 'stats:lillietest:OutOfRangePHigh');

% Plot histograms with normal theory curve superimposed
pvals = zeros(numGroups, numel(tests));
for i = 1:numGroups
    figure()
    k = 0;
    for j = 1:numel(labels)
        if any(strcmpi(labels{j}, tests))
            k = k + 1;
            subplot(3, 3, k), histfit(data{i}(:, j))
            title([num2str(doses(i)), ' \muM H_2O_2 ', labels{j}])
            [~, p] = lillietest(data{i}(:, j));
            pvals(i, k) = p;
        end
    end
end

%% ACA Dose Response
close all
clc

% Extract stats
warning('off', 'MATLAB:catenate:DimensionMismatch');
S = zeros(0, numel(tests));
for i = 1:numGroups
    S = [S; stats{i}]; %#ok<AGROW>
end

% Plot results
figure()
k = 0;
numStats = numel(statLabels);
N = 400; % number of cells scored per dose
for i = 1:numel(labels)
    if any(strcmpi(labels{i}, tests))
        k = k + 1;
        subplot(3, 3, k)
        errorbar(doses, S(3:numStats:end, i), S(4:numStats:end, i) / sqrt(N)) % Bars are SEM
        title(labels{i})
    end
end


%% Halo Assay Data Analysis
close all
clc

% Import data -- some things are changed for each experiment
doses = [0, 25, 50, 100]; % Set doses used
numGroups = length(doses);
data = cell(1, numGroups);
stats = cell(1, numGroups);
in = Imports();
for i = 1:numGroups
    filename1 = sprintf('07.30.14 %d rads D1 (EB) damage', doses(i));
    filename2 = sprintf('07.30.14 %d rads D2 (EB) damage', doses(i));
    filename3 = sprintf('07.31.14 %d rads D1 (EB) damage', doses(i));
    filename4 = sprintf('07.31.14 %d rads D2 (EB) damage', doses(i));
    if i == 1
        [data1, ~, labels, statLabels] = in.haloData(filename1);
    else
        [data1, ~] = in.haloData(filename1);
    end
    [data2, ~] = in.haloData(filename2);
    [data3, ~] = in.haloData(filename3);
    [data4, ~] = in.haloData(filename4);
    combData = [data1; data2; data3; data4];
    data{i} = combData;
    N = size(combData, 1);
    stats{i} = [mean(combData); std(combData); std(combData) / sqrt(N)]; % Last is SEM
    fprintf(1, 'Dose: %d\nN: %d\n', doses(i), N);
end

% Format data into struct
haloData = struct();
for i = 1:numGroups
    for j = 1:numel(labels)
        % Remove illegal struct name characters
        L = labels{j};
        illegal = ' []()/';
        for il = illegal
            L(L == il) = '_';
        end
        L(L == '%') = 'P'; % for Tail DNA[%]
        
        % Build struct
        D = sprintf('d%d', doses(i));
        haloData.(D).(L).data = data{i}(:, j);
        for k = 1:numel(statLabels)
            haloData.(D).(L).(statLabels{k}) = stats{i}(k, j);
        end
    end
end


%% FHA Histograms
close all
clc

warning('off', 'stats:lillietest:OutOfRangePLow');
warning('off', 'stats:lillietest:OutOfRangePHigh');

% Plot histograms with normal theory curve superimposed
pvals = zeros(numGroups, numel(labels));
for i = 1:numGroups
    figure()
    for j = 1:numel(labels)
        subplot(3, 3, j), histfit(data{i}(:, j))
        title([num2str(doses(i)), ' \muM H_2O_2 ', labels{j}])
        [~, pvals(i, j)] = lillietest(data{i}(:, j));
    end
end


%% FHA Dose Response
close all
clc

% Extract stats
warning('off', 'MATLAB:catenate:DimensionMismatch');
numStats = 3;
S = zeros(0, numStats);
for i = 1:numGroups
    S = [S; stats{i}]; %#ok<AGROW>
end

% Plot results
figure()
for i = 1:numel(labels)
    subplot(3, 3, i)
    errorbar(doses, S(1:numStats:end, i), S(3:numStats:end, i)) % Bars are SEM
    title(labels{i})
end



%% Compare FHA and ACA Results
% For this comparison, I am matching comet assay parameters to those I have created for the FHA in
% order to assess their comparability. 
%
% Comparisons:
%   
%   NDF -------> Comet Total Area / (Comet Total Area - Tail Area)
%   HD --------> Comet Extent / (Comet Extent - Tail Extent)
%   HM, AHM ---> Tail Distributed Moment, Tail Extent Moment
%   IHI, AHI --> Comet Total Intensity  
%   RHI -------> Tail DNA[%], Tail (Olive) Moment / 100
close all
clc

% Create new ACA metrics
for i = 1:numGroups
    
    D = sprintf('d%d', doses(i));
    dose = cometData.(D);
    acaFHA = dose.Comet_total_area.data ./ (dose.Comet_total_area.data - dose.Tail_area__Singh_.data);
    acaHD = dose.Comet_extent.data ./ (dose.Comet_extent.data - dose.Tail_extent.data);
    
    dose.acaFHA.data = acaFHA;
    dose.acaFHA.Mean = mean(acaFHA);
    dose.acaFHA.StdDev = std(acaFHA);
    dose.acaFHA.Median = median(acaFHA);
    
    dose.acaHD.data = acaHD;
    dose.acaHD.Mean = mean(acaHD);
    dose.acaHD.StdDev = std(acaHD);
    dose.acaHD.Median = median(acaHD);
    cometData.(D) = dose;
end

% Make comparison plot
figure()
labFHA = {'NDF', 'HD', 'HM', 'AHM', 'IHI', 'AHI', 'RHI', 'RHI'};
labACA = {'acaFHA', 'acaHD', 'Tail_distributed_moment__DNA_migration_', 'Tail_extent_moment',...
          'Comet_total_intensity', 'Comet_total_intensity', 'Tail_DNA_P_', 'Tail__Olive__moment__100'};
nameACA = {'Comet Area Ratio', 'Comet Distance', 'Tail Distributed Moment', 'Tail Distributed Moment',...
           'Comet Total Intensity', 'Comet Total Intensity', 'Tail DNA %', 'Olive Moment'};

for n = 1:numel(labFHA)
    % Get ACA metric mean and SEM
    meansACA = zeros(1, numGroups);
    semsACA = zeros(1, numGroups);
    for i = 1:numGroups
        meansACA(i) = cometData.(sprintf('d%d', doses(i))).(labACA{n}).Mean;
        semsACA(i) = cometData.(sprintf('d%d', doses(i))).(labACA{n}).StdDev / 10;
    end
    
    % Get FHA metric mean and SEM
    meansFHA = zeros(1, numGroups);
    semsFHA = zeros(1, numGroups);
    for i = 1:numGroups
        meansFHA(i) = haloData.(sprintf('d%d', doses(i))).(labFHA{n}).mean;
        semsFHA(i) = haloData.(sprintf('d%d', doses(i))).(labFHA{n}).SEM;
    end
    
    % Normalization factor to make plots more comparable
    meanACA = mean(meansACA);
    meanFHA = mean(meansFHA);
    
    % Plot normalized dose response together
    subplot(4, 4, 2 * n - 1)
    errorbar(doses, meansACA / meanACA, semsACA / meanACA, 'b'), hold on
    errorbar(doses, meansFHA / meanFHA, semsFHA / meanFHA, 'r')
    title(sprintf('%s vs %s', labFHA{n}, nameACA{n}))
    xlabel('H_2O_2 Dose (/muM)')
    ylabel('Normalized Metric Score')
    legend(nameACA{n}, labFHA{n})
    
    % Plot correlation between metrics
    subplot(4, 4, 2 * n)
    [R, slope, intercept] = regression(meansFHA, meansACA);
    scatter(meansFHA, meansACA, 'k'), hold on
    plot(meansFHA, slope * meansFHA + intercept, 'k')
    title(sprintf('Correlation between %s and %s. R^2 = %.3f', labFHA{n}, nameACA{n}, R^2))
    xlabel(labFHA{n}), ylabel(nameACA{n})
end


%% Complete Comparison
% Compare every metric to each other and plot a heat map of R^2 values. In case interesting
% parameters have correlations that I would not guess.
close all
clc

% Make comparison plot
labFHA = {'NDF', 'HD', 'HM', 'AHM', 'IHI', 'AHI', 'RHI'};
labACA = {'Comet_extent', 'Comet_total_intensity', 'Comet_total_area', 'Tail_DNA_P_', 'Tail_extent',...
          'Tail_distributed_moment__DNA_migration_', 'Tail_extent_moment', 'Tail__Olive__moment'...
           'Tail_integrated_intensity__Singh_', 'Tail_integrated_intensity_ratio__Singh_',...
           'Tail_area__Singh_', 'Tail_break_Number__Singh_'};
nameACA = {'Comet Extent', 'Comet Intensity', 'Comet Area', 'Tail DNA %', 'Tail Extent',...
           'Tail Distributed Moment', 'Tail Extent Moment', 'Tail (Olive) Moment',...
           'Tail Integrated Intensity', 'TII Ratio', 'Tail Area', 'Tail Break #'};

corrs = zeros(numel(labFHA), numel(labACA));
for n = 1:numel(labFHA)
    for m = 1:numel(labACA)
        % Get FHA metric mean
        meansFHA = zeros(1, numGroups);
        for i = 1:numGroups
            meansFHA(i) = haloData.(sprintf('d%d', doses(i))).(labFHA{n}).mean;
        end
        
        % Get ACA metric mean
        meansACA = zeros(1, numGroups);
        for i = 1:numGroups
            meansACA(i) = cometData.(sprintf('d%d', doses(i))).(labACA{m}).Mean;
        end
        
        R = regression(meansFHA, meansACA);
        corrs(n, m) = R * abs(R);
    end
end

h = HeatMap(corrs);
h.ColumnLabels = nameACA;
h.RowLabels = labFHA;


%% Write corrs to file
clc

fileID = fopen(strcat(pwd, '\Data\Correlations.csv'), 'wt');
fprintf(fileID, ',');
fprintf(fileID, '%s,', nameACA{1:end-1});
fprintf(fileID, '%s\n', nameACA{end});
for i = 1:size(corrs, 1)
    fprintf(fileID, '%s,', labFHA{i});
    fprintf(fileID, '%d,', corrs(i, 1:end-1));
    fprintf(fileID, '%d\n', corrs(i, end));
end
fclose(fileID);


%% Cell Counts
close all
clc

doses = [0, 25, 50, 100, 200, 500];
hours = [0, 24, 48];

countData = struct();
data = xlsread('Counting Data.xlsx');
for i = 2:2:size(data, 1)
    for j = 1:size(data, 2)
        D = sprintf('d%d', doses(i / 2));   % Dose
        H = sprintf('hr%d', hours(j));      % Time point
        countData.(D).(H).live = data(i - 1, j);
        countData.(D).(H).dead = data(i, j);
    end
end


%% CC Dose Response
close all
clc

colors = 'mrgcbk';

figure(), hold on
for i = 2:2:size(data, 1)
    plot(hours, data(i - 1, :) / data(i - 1, 1), colors(i / 2))
end
for i = 2:2:size(data, 1)
    plot(hours, data(i, :) / data(i, 1), [colors(i / 2), '--'])
end
xlabel('Hours')
ylabel('# of Cells Relative to 0 hr')
leg = cell(1, length(doses));
for i = 1:length(doses)
    leg{i} = sprintf('%d uM H_2O_2', doses(i));
end
legend(leg)


