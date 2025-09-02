% Ensure column vectors
function X = ensurecolumn(X)
% This function ensures that an input vector is a column vector.

    if length(X) == width(X)
        X = X.';
    end
    
    return
end