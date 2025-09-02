% Stitch together BonFIRE spectra
function [XSTITCH,YSTITCH] = bfstitch(XVAL1,YVAL1,XVAL2,YVAL2)
%%% This function stitches together two parts of a BonFIRE spectrum.
% (leaves YVAL1 intensities unchanged)

    % Ensure column vectors
    XVAL1 = ensurecolumn(XVAL1); XVAL2 = ensurecolumn(XVAL2);
    YVAL1 = ensurecolumn(YVAL1); YVAL2 = ensurecolumn(YVAL2);

    % Make sure lengths match
    [XVAL1,YVAL1] = lengthmatch(XVAL1,YVAL1);
    [XVAL2,YVAL2] = lengthmatch(XVAL2,YVAL2);

    % Use knnsearch to find stitch point
    XVAL1index = mode(knnsearch(XVAL1,XVAL2)); % freq in XVAL1 closest to the freqs in XVAL2
    XVAL2index = mode(knnsearch(XVAL2,XVAL1));
    normfactor = YVAL1(XVAL1index)./YVAL2(XVAL2index);
    
    YVAL2 = YVAL2.*normfactor;
    
    % Remove stitch point (double-counted)
    YVAL2(XVAL2index) = [];
    XVAL2(XVAL2index) = [];
    
    % Stitch together
    XSTITCH = [XVAL1; XVAL2];
    YSTITCH = [YVAL1; YVAL2];

    [XSTITCH,sortinds] = sort(XSTITCH);
    YSTITCH = YSTITCH(sortinds);
    return
end
