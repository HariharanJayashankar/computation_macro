function x=pentle(a,b,c,d,e,r);
%PENTLE solves the following pentadiagonal linear system of equations:
%
%  [c(1) d(1) e(1)   0    0    0  ...   0      0      0     0   ] 
%  [b(2) c(2) d(2) e(2)   0    0  ...   0      0      0     0   ] 
%  [a(3) b(3) c(3) d(3) e(3)   0  ...   0      0      0     0   ] 
%  [  0  a(4) b(4) c(4) d(4) e(4) ...   0      0      0     0   ] 
%  [             .                   .                 .        ]  x = r
%  [               .                   .                 .      ] 
%  [                 .                   .                 .    ] 
%  [  0    0    0    0    0    0 ... a(n-1) b(n-1) c(n-1) d(n-1)] 
%  [  0    0    0    0    0    0 ...   0    a(n)   b(n)   c(n)  ] 
%
%

%  Ellen McGrattan

n = size(c);
if min([min(abs( [a(3:n),e(1:n-2)] )),min(abs( [b(2:n),d(1:n-1)] )), ...
        min(abs( c ))])<1e-20;
  error('The five diagonals of A must have nonzero elements');
end; 
%
%  First row reduction:
%
for i=n-2:-1:1;
  ed=e(i)/d(i+1);
  b(i)=b(i)-a(i+1)*ed;
  c(i)=c(i)-b(i+1)*ed;
  d(i)=d(i)-c(i+1)*ed;
  r(i,:)=r(i,:)-r(i+1,:)*ed;
end;
%
%  Second row reduction
%
for i=3:n;
  ab=a(i)/b(i-1);
  b(i)=b(i)-c(i-1)*ab;
  c(i)=c(i)-d(i-1)*ab;
  r(i,:)=r(i,:)-r(i-1,:)*ab;
end;
x=trile(b,c,d,r);
