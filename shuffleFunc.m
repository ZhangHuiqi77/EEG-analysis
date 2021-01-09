function xs = shuffleFunc(x,ns)
    % Huiqi
    % 02/01/2021
    % This function generate shuffled x based on randperm(irreplaceable)
    % Input:
        % x: data, can be many trials, each line is a trial
        % ns: the number you want to shuffle for every trialed line
    % Output:
        % xs: shuffled x, size is [ns*size(x,1),size(x,2)]
    xs = zeros(ns*size(x,1),size(x,2));			%Vector to hold shuffled x results
    for m = 1:size(x,1)
        for n = 1:ns				%For each surrogate,
            xs(n+ns*(m-1),:) = x(m,randperm(size(x,2)));	%Resample by y-index
        end
    end
end