% Save_tif v2 (from HW)
function Save_tif(fname, imageStack, sliceLabels)
%%% This function writes a 3D array into a hyperstack .tif file.

%%% v2 - extends from max stack length of 5 to any size stack

    TIF = Tiff(fname, 'w');
    for k = 1:size(imageStack, 3)
        % Set TIFF tags for floating-point storage
        tagStruct.ImageLength = size(imageStack, 1);
        tagStruct.ImageWidth = size(imageStack, 2);
        tagStruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagStruct.BitsPerSample = 64; % 64-bit double precision
        tagStruct.SamplesPerPixel = 1;
        tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagStruct.SampleFormat = Tiff.SampleFormat.IEEEFP; % Floating-point format
        tagStruct.Software = 'MATLAB';
        % Set tags and write the slice
        if ~isempty(sliceLabels)
          tagStruct.PageName = sliceLabels{k};
        end
        TIF.setTag(tagStruct);
        TIF.write(imageStack(:, :, k));
        if k < size(imageStack, 3)
           TIF.writeDirectory(); % Create a new directory for the next slice
        end
    end
    TIF.close();
    return
end
