function x=trile(a,b,c,r);
%TRILE solves the following tridiagonal linear system(s) of equations:
%
%  [b(1) c(1)   0    0    0    0 ...   0      0      0   ] 
%  [a(2) b(2) c(2)   0    0    0 ...   0      0      0   ] 
%  [  0  a(3) b(3) c(3)   0    0 ...   0      0      0   ] 
%  [           .                     .                   ]  
%  [             .                     .                 ] x = r
%  [               .                     .               ] 
%  [  0    0    0    0    0    0 ... a(n-1) b(n-1) c(n-1)] 
%  [  0    0    0    0    0    0 ...   0    a(n)   b(n)  ] 
%
%  or A x = r,  where A: m x n, r: rm x rn, x: rm x rn
%

%
%  Row reduction:
%
n=length(r);
for i=n-1:-1:1;
  cb=c(i)/b(i+1);
  b(i)=b(i)-a(i+1)*cb;
  r(i)=r(i)-r(i+1)*cb;
end;
x(1,1)=r(1)/b(1);
for i=2:n;
  x(i,1)=(r(i)-a(i)*x(i-1))/b(i);
end;
