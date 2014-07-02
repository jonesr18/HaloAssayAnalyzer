% Ross Jones
% Singh Lab
% University of Washington
% DNA Damage Analysis using multithresh

function [] = DNA_damage_multithresh(sequence, numFiles)
    
    % Threshold each image in seqence
    % PROBLEM: This is the command causing the upside down blue kid image to show up
    seqThresh = zeros(size(sequence), class(image));
    for i = 1:numFiles
        image_i = sequence(:, :, i);
        thresh = multithresh(image_i, 2);                   % Implement multithresh
        seqThresh(:, :, i) = imquantize(image_i, thresh);   % apply thresholds to obtain segmented image
    end
    
    % Surface area comparison

    % seqThresh is in form [1 2 3] for [background, halo, nucleus]. In order to
    % extract each part, we must use subtractions to set certain threshold groups to
    % zero and multiply to remove such regions. Then the matrix must be adjusted so
    % that the desired region is set to 1, and that regions is summed to get SA
    
    seqThresh   = uint8(seqThresh);
    nucThresh   = seqThresh - 2;
    nucSurface  = sum(sum(sum(nucThresh)));
    haloThresh  = (seqThresh - 1) .* uint8(abs(int8(seqThresh) - 3));
    haloSurface = sum(sum(sum(haloThresh)));
    surfDamage  = (haloSurface + nucSurface) / nucSurface;

    % Relative intensity comparison

    % Use the nuc and halo thesholds to isolate desired regions and compare relative
    % intensities. This method takes into account the actual intensity, which may be
    % eschewed by the dying and imaging of each slide, but could be more accurate
    % otherwise because it takes into account DNA density in areas.

    nucInten  = sum(sum(sum(nucThresh .* sequence)));
    haloInten = sum(sum(sum(haloThresh .* sequence)));
    intDamage = haloInten / nucInten;

    fprintf('Surface area comparison damage: %.4f\n', surfDamage)
    fprintf('Relative intensity comparison damage: %.4f\n', intDamage)
end







