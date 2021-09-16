function result = logsumexp(x,y,dim)

  if (nargin == 1)
    y = [];
    dim = 1;
  elseif (nargin == 2)
    if ~all(size(x)==size(y))
      error('Sizes must match when passing two arguments');
    end
  elseif (nargin == 3)
    if ~isempty(y)
      error('Second argument must be empty when passing three arguments');
    end
  end
  
  if isempty(y)
    maxval = max(x,[],dim);
    result = log(sum(exp(x-maxval),dim)) + maxval;
  else
    maxval = max(x,y);
    result = log(exp(x-maxval)+exp(y-maxval)) + maxval;
  end
  
  notfinite = ~isfinite(maxval);
  result(notfinite&(maxval>0)) = inf;
  result(notfinite&(maxval<0)) = -inf;
  
end
