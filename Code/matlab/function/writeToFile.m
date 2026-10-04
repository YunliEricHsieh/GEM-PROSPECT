function writeToFile(fileName, rowData)
%WRITETOFILE Append a cell row to an existing CSV output without quoting.
%   fileName is a writable path; rowData contains identifiers and flux cells.
%   Uses convertToString; identifiers must not contain commas or newlines.
    % Convert all elements to strings
    formattedRow = cellfun(@convertToString, rowData, 'UniformOutput', false);

    % Join the row as a comma-separated string and append it to the file
    fid = fopen(fileName, 'a');  % Open file in append mode
    fprintf(fid, '%s\n', strjoin(formattedRow, ','));  % Write the formatted row
    fclose(fid);  % Close the file
end
