function [mean_value, final_error] = meanAndError(values, errors)
    n = numel(values);
    mean_value = mean(values);
    sumsqErr = 0;
    for i=1:length(errors)
        sumsqErr = sumsqErr + errors(i)*errors(i);
    end
    final_error = sqrt(sumsqErr/(n*n)+std(values)*std(values));
end