%
% Solver for tri-diagonal matrices.
%
% T must be in compact form: the sub-diagonal elements in the first
% column, the diagonal in the second column, and the super-diagonal
% in the third column.
%
% Note that T(1,1) and T(n,3) are never accessed, i.e. the
% subdiagonal entries start on the second row, and the
% superdiagonal elements end on the (n-1)st row.
%
function [x]=trisolve(T,b)

if nargin ~= 2
  error('trisolve:: two input arguments required.')
end
if nargout ~= 1
  error('trisolve: one output argument required.')
end

[n,m] = size(T);
if m ~= 3
  error('trisolve: matrix must have three columns.')
end

[N,M] = size(b);

if N ~= n
  error('trisolve: rhs and matrix must have same number of rows.');
end

work = zeros(n,1);

work(1) = T(1,2);
x(1,:)  = b(1,:);

% Forward sweep.
for i=2:n
  beta    = T(i,1)/work(i-1);
  x(i,:)  = b(i,:) - beta*x(i-1,:);
  work(i) = T(i,2) - beta*T(i-1,3);
end

x(n,:) = x(n,:)/work(n);

% Backward sweep.
for i=n-1:-1:1
  x(i,:) = (x(i,:) - T(i,3)*x(i+1,:)) / work(i);
end

