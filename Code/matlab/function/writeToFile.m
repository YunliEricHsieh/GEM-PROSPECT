function writeToFile(fileName, rowData)
    % Convert all elements to strings
    formattedRow = cellfun(@convertToString, rowData, 'UniformOutput', false);

    % Join the row as a comma-separated string and append it to the file
    fid = fopen(fileName, 'a');  % Open file in append mode
    fprintf(fid, '%s\n', strjoin(formattedRow, ','));  % Write the formatted row
    fclose(fid);  % Close the file
end
