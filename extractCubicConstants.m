function extractCubicConstants(logFilePath, outputFilePath)
    % Extract cubic force constants from a Gaussian log file
    % Args:
    %   logFilePath (str): Path to the Gaussian log file
    %   outputFilePath (str): Path to save the extracted cubic force constants

    % Open the log file
    fid = fopen(logFilePath, 'r');
    if fid == -1
        error('Error opening the file: %s', logFilePath);
    end
    
    % Initialize variables
    cubicConstants = [];
    inCubicSection = false;
    
    % Loop through the file line by line
    while ~feof(fid)
        line = fgetl(fid); % Read a line from the file
        
        % Detect the start of the cubic force constants section
        if contains(line, 'Cubic Force Constants')
            inCubicSection = true;
            continue;
        end
        
        % Exit the cubic section if a blank line is encountered
        if inCubicSection && isempty(strtrim(line))
            break;
        end
        
        % Parse the cubic force constants if in the correct section
        if inCubicSection
            % Use textscan to parse lines with the expected format
            formatSpec = '%d %d %d %d %f'; % Mode, i, j, k, Value
            data = textscan(line, formatSpec);
            
            % Check if parsing was successful
            if ~isempty(data{1})
                cubicConstants = [cubicConstants; ...
                                  data{1}, data{2}, data{3}, data{4}, data{5}];
            end
        end
    end
    
    % Close the file
    fclose(fid);
    
    % Write the extracted data to the output file
    outputFile = fopen(outputFilePath, 'w');
    if outputFile == -1
        error('Error creating the output file: %s', outputFilePath);
    end
    
    % Write header
    fprintf(outputFile, 'Mode\t i\t j\t k\t Value\n');
    % Write data
    for i = 1:size(cubicConstants, 1)
        fprintf(outputFile, '%d\t %d\t %d\t %d\t %.6f\n', cubicConstants(i, :));
    end
    
    % Close the output file
    fclose(outputFile);
    
    fprintf('Extracted %d cubic constants. Saved to %s\n', size(cubicConstants, 1), outputFilePath);
end