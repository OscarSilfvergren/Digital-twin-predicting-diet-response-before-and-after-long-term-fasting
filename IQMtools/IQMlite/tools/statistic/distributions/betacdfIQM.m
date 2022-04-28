function cdf_result = betacdfIQM(x, a, b)
% Cumulative density function of the Beta distribution
%
% USAGE:
% ======
% cdf_result = betacdfIQM(x, a, b)
%
% For each element of 'x', returns the CDF at 'x' of the beta
% distribution with parameters 'a' and 'b', i.e.,
% PROB (beta ('a', 'b') <= 'x').

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

  if (nargin ~= 3)
    error('Incorrect number of input arguments'); 
  end

  if (~isscalar (a) || ~isscalar(b))
    [retval, x, a, b] = common_size (x, a, b);
    if (retval > 0)
      error ('x, a and b must be of common size or scalar');
    end
  end

  sz = size(x);
  cdf_result = zeros (sz);

  k = find (~(a > 0) | ~(b > 0) | isnan (x));
  if (any (k))
    cdf_result (k) = NaN;
  end

  k = find ((x >= 1) & (a > 0) & (b > 0));
  if (any (k))
    cdf_result (k) = 1;
  end

  k = find ((x > 0) & (x < 1) & (a > 0) & (b > 0));
  if (any (k))
    if (isscalar (a) && isscalar(b))
      cdf_result (k) = betainc(x(k), a, b);
    else
      cdf_result (k) = betainc(x(k), a(k), b(k));
    end
  end

end
