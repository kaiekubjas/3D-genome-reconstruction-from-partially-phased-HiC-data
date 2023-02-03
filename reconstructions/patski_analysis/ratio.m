function ratio_vector = ratio(v)

% Function that divides a vector by the sum of the absolute value of all non-nan entries

if nansum(abs(v)) == 0
    ratio_vector = ones(size(v))/length(v);
else
    ratio_vector = v*(1/nansum(abs(v)));
end

end
