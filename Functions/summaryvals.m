% Summary values
function [SIGPEAK,PEAKSIGN,PEAKSD,INTEGSIG,INTEGSD,NOISE,SNR,BASESD,SNRV2,SNRSWEEP,INTSNR] = ...
    summaryvals(CORRSIG,CORRSIGBASE,SIGSDS)
%%% This function calculates summary statistics/values for data after 
% baseline correction.

    % Determine if signal peak is negative or positive
    abscorrsig = abs(CORRSIG);
    flipcorrsig = -1*CORRSIG;
    corrabsdist = sum((CORRSIG-abscorrsig).^2);
    flipabsdist = sum((flipcorrsig-abscorrsig).^2);
    if corrabsdist < flipabsdist
        [SIGPEAK,maxindex] = max(CORRSIG); % corrsig is closer to abscorrsig
        PEAKSIGN = 1; % peak is positive
    else
        [SIGPEAK,maxindex] = min(CORRSIG); % flipcorrsig is closer to abscorrsig
        PEAKSIGN = -1; % peak is negative
    end
    % Calculate summary values
    PEAKSD = SIGSDS(maxindex);
    INTEGSIG = sum(CORRSIG);
    INTEGSD = mean(SIGSDS,'all');
    NOISE = range(CORRSIGBASE);
    BASESD = std(CORRSIGBASE);
    SNR = SIGPEAK/NOISE;
    SNRV2 = SIGPEAK/BASESD;
    SNRSWEEP = CORRSIG/NOISE;
    INTSNR = sum(SNRSWEEP);
    return
end
