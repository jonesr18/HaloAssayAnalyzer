% Ross Jones
% Singh Lab
% University of Washington
% DNA Damage Analysis using multithresh

function [] = DNA_damage_kmeans(files, filePath, numFiles, seqGray)
    
    % Preallocate arrays
    image     = imread(fullfile(filePath, files{1}));
    sequence  = zeros([size(image) numFiles], class(image));
    seqLabels = zeros([size(image(:, :, 1)) numFiles], class(image));
    
    % Threshold each image in seqence
    for i = 1:numFiles
        image_i = imread(fullfile(filePath, files{i})); 
        sequence(:, :, :, i) = image_i;
        
        % Convert to a*b* space
        cform = makecform('srgb2xyz');
        Lab = applycform(image_i, cform);
        ab = double(Lab(:, :, 1:3));
        nrows = size(ab, 1);
        ncols = size(ab, 2);
        image_i = reshape(ab, nrows*ncols, 3);
        
        % Perform clustering (3 times to avoid local minima)
        try
            IDX = kmeans(image_i, 3, 'distance', 'sqEuclidean', 'Replicates', 3);
        catch me
            disp(me)
        end
        
        % Label pixels with k-means results
        seqLabels(:, :, i) = reshape(IDX, nrows, ncols);
    end
    
    % Surface area comparison

    % seqThresh is in form [1 2 3] for [background, halo, nucleus]. In order to
    % extract each part, we must use subtractions to set certain threshold groups to
    % zero and multiply to remove such regions. Then the matrix must be adjusted so
    % that the desired region is set to 1, and that regions is summed to get SA
    
    seqLabels   = uint8(seqLabels);
    nucThresh   = seqLabels - 2;
    nucSurface  = sum(sum(sum(nucThresh)));
    haloThresh  = (seqLabels - 1) .* uint8(abs(int8(seqLabels) - 3));
    haloSurface = sum(sum(sum(haloThresh)));
    surfDamage  = haloSurface / nucSurface;

    % Relative intensity comparison

    % Use the nuc and halo thesholds to isolate desired regions and compare relative
    % intensities. This method takes into account the actual intensity, which may be
    % eschewed by the dying and imaging of each slide, but could be more accurate
    % otherwise because it takes into account DNA density in areas.
    
    nucInten  = sum(sum(sum(nucThresh .* seqGray)));
    haloInten = sum(sum(sum(haloThresh .* seqGray)));
    intDamage = haloInten / nucInten;

    fprintf('Surface area comparison damage: %.4f\n', surfDamage)
    fprintf('Relative intensity comparison damage: %.4f\n', intDamage)
end







