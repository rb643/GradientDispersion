function [s,n] = sumsqr(x)
%SUMSQR Sum of squared elements of a matrix or matrices.
%
%  [S,N] = <a href="matlab:doc sumsqr">sumsqr</a>(X) returns S the sum of all squared finite elements in M,
%  and the number of finite elements N.  M may be a numeric matrix or a
%  cell array of numeric matrices.
%
%  For example:
%
%    s = <a href="matlab:doc sumsqr">sumsqr</a>([1 2;3 4])
%    [s,n] = <a href="matlab:doc sumsqr">sumsqr</a>({[1 2; NaN 4], [4 5; 2 3]})
%
%  See also MEANSQR, SUMABS, MEANABS.

% Mark Beale, 1-31-92
% Copyright 1992-2010 The MathWorks, Inc.

if nargin < 1,error(message('nnet:Args:NotEnough'));end

if isreal(x)
  notFinite = find(~isfinite(x));
  x(notFinite) = 0;
  s = x.*x;
  s = sum(s(:));
  n = numel(x) - length(notFinite);
elseif iscell(x)
  numx = numel(x);
  s = zeros(1,numx);
  n = zeros(1,numx);
  for i=1:numx
    [s(i),n(i)] = sumsqr(x{i});
  end
  s = sum(s);
  n = sum(n);
else
  error(message('nnet:NNData:NotNumOrCellNum'))
end