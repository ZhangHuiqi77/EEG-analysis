function [xm,xmin,xmax] = calstat(x)
    xm = mean(x);
    xmin = min(x);
    xmax = max(x);
end