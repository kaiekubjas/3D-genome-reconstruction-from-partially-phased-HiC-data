function [ua_pairs,ambig_pairs] = partition(n,m,options)

arguments
    n
    m
    options.distribution
end

if ~isfield(options,'distribution') || strcmp(options.distribution,'random')
    ua_pairs = randperm(n);
    ua_pairs = sort(ua_pairs(1:n-m));
elseif strcmp(options.distribution,'even')
    ua_pairs = 1:floor(n/(n-m)):n;
    ua_pairs = ua_pairs(1:n-m);
elseif strcmp(options.distribution,'block')
    start = randi([1,m+1]);
    ua_pairs = start:(start+(n-m)-1);
end

ambig_pairs = setdiff(1:n,ua_pairs);


end