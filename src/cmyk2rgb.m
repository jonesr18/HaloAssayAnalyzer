% Summary
%   Description
function rgb = cmyk2rgb(cmyk)
    validateattributes(cmyk, {'uint8'}, {'3d'}, mfilename, 'CMYK image', 1);
    dim3 = size(cmyk, 3);
    if dim3 ~= 4
        error('CMYK image not in correct format (N x M x 4)')
    end
    inprof = iccread('USSheetfedCoated.icc');
    outprof = iccread('sRGB.icm');
    rgb = applycform(cmyk, makecform('icc', inprof, outprof));
end