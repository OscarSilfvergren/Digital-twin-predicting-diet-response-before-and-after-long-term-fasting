function inv = norminvIQM(x, m, s)
% Quantile function of the normal distribution
%
% USAGE:
% ======
% pdf = norminvIQM(x, m, s)
%
% For each element of 'x', compute the quantile (the inverse of the
% CDF) at 'x' of the normal distribution with mean 'm' and
% standard deviation 's'.
%
% Default values are 'm' = 0, 's' = 1.

% <<<COPYRIGHTSTATEMENT - IQM TOOLS LITE>>>

  if (nargin ~= 1 && nargin ~= 3)
    error('Incorrect number of input arguments'); 
  end

  if (nargin == 1)
    m = 0;
    s = 1;
  end

  if (~isscalar (m) || ~isscalar (s))
    [retval, x, m, s] = common_size (x, m, s);
    if (retval > 0)
      error ('x, m and s must be of common size or scalars');
    end
  end

  sz = size (x);
  inv = zeros (sz);

  if (isscalar (m) && isscalar (s))
    if (find (isinf (m) | isnan (m) | ~(s > 0) | ~(s < Inf)))
      inv = NaN * ones (sz);
    else
      inv =  m + s .* stdnormalinvIQM(x);
    end
  else
    k = find (isinf (m) | isnan (m) | ~(s > 0) | ~(s < Inf));
    if (any (k))
      inv(k) = NaN;
    end

    k = find (~isinf (m) & ~isnan (m) & (s > 0) & (s < Inf));
    if (any (k))
      inv(k) = m(k) + s(k) .* stdnormalinvIQM (x(k));
    end
  end

  k = find ((s == 0) & (x > 0) & (x < 1));
  if (any (k))
    inv(k) = m(k);
  end

  inv((s == 0) & (x == 0)) = -Inf;
  inv((s == 0) & (x == 1)) = Inf;

return
