% please input a line/column of data, and the function calbs will return
% baseline of the data with 98% datapoint included.


function [bs]= calbs(data)

    data = sort(data);
    N = length(data);
    bs_min = data(N * 0.02);
    bs_max = data(N * 0.98);
    bs = bs_max - bs_min;

end