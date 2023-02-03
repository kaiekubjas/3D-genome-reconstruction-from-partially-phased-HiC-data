function varargout = unpack_vector(v)
    for i = 1:nargout
        varargout(i) = {v(i)}; %#ok<AGROW>
    end
end
