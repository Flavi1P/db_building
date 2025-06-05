function [new_x, new_y, new_n, bin_ixs] = binmedian_edges(old_x,old_y,x_edges)
% [new_x, new_y, new_n, bin_ixs] = binmedian_edges(old_x,old_y,x_edges)
%
% bins y data into medians, grouped by x data bounds defined in x_edges

rowvec = isrow(old_x);

if rowvec
    old_x = old_x';
    old_y = old_y';
end

n = length(x_edges)-1;

new_x = nan(n,1);
new_y = nan(n,size(old_y,2));
new_n = nan(n,size(old_y,2));

bin_ixs = nan(size(old_x));

for this_bin = 1:n
    good = old_x >= x_edges(this_bin) & old_x < x_edges(this_bin+1);
    if any(good)
        bin_ixs(good) = this_bin;
        new_x(this_bin) = nanmedian(old_x(good),1);
        new_y(this_bin,:) = nanmedian(old_y(good,:),1);
        new_n(this_bin,:) = sum(~isnan(old_y(good,:)),1);
    end
end
    