function str = convertToString(element)
    % Convert different types of elements to strings
    if iscell(element) && numel(element) == 1  % Single-element cell
        element = element{1};
    end
    if isnumeric(element)
        str = num2str(element);  % Convert numeric values to strings
    elseif ischar(element) || isstring(element)
        str = char(element);  % Ensure it's a character array
    else
        str = 'NA';  % Default for unsupported types
    end
end