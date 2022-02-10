function [v,x0,r] = fitCylinder (xyz)
  x = xyz(:,1);
  y = xyz(:,2);
  z = xyz(:,3);
  A = [x.^2 x.*y x.*z y.^2 y.*z z.^2 x y z];
  p = A \ ones(size(x));
%  sum( ( A*p - ones(size(x)) ).^2 )
  P = zeros(3);
  P([1 2 3 5 6 9]) = p(1:6);
  P = (P+P')/2;
  [V,e] = eig(P);
  e = diag(e);
  [r,n] = sort(e);
  v = V(:,n(1));
  B = -(p(7:9)'*V)'./e/2;
  r = r/(1+e' * B.^2);
  r = r(2:3).^-.5;
  x0 = V(:,n(2:3))*B(n(2:3));
end
