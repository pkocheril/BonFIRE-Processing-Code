% Normalize a signal
function NORMSIG = normsweep(SIGNAL)
%%% This function background-subtracts and normalizes a vector.
    
    NORMSIG = (SIGNAL-min(SIGNAL))./(max(SIGNAL)-min(SIGNAL));
    return
end