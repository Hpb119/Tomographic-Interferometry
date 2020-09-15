% Input: sparse system matrix A, data b.
% Output: SIRT reconstruction x.

tic
x = zeros(N * N, 1);
[rows cols] = size(A);
C = sparse(1 : cols, 1 : cols, 1 ./ sum(A));
R = sparse(1 : rows, 1 : rows, 1 ./ sum(A'));
CATR = C * A' * R;
for i = 1 : 100
  x = x + CATR * (B - A * x);
  i
end
imagesc(reshape(x,N,N)), colorbar
toc