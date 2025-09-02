% Make sure lengths match
function [X,Y] = lengthmatch(X,Y)
%%% This function ensures that two vectors are the same length (longer
% vector is trimmed if there's a mismatch).

    if length(X) ~= length(Y)
        if length(X) < length(Y) % Y is longer
            Y = Y(1:length(X)); % Y trimmed to length(X)
        else % X is longer
            X = X(1:length(Y)); % X trimmed to length(Y)
        end
    end
    
    return
end
