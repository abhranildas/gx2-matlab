function gx2_check_method(method,w,s)
% GX2_CHECK_METHOD Validates that the parameters are compatible with the
% requested method. Ruben's method and the ellipse approximation both
% require all weights w to be the same sign and no normal term (s=0).
% Throws an error otherwise. Used by gx2cdf and gx2pdf.

if s || ~(all(w>0)||all(w<0))
    if strcmpi(method,'ruben')
        error("Ruben's method can only be used when all w are the same sign and s=0.")
    elseif strcmpi(method,'ellipse')
        error("The ellipse approximation can only be used when all w are the same sign and s=0.")
    end
end
end
