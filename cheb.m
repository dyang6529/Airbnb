% Could also write chebyshev recursively. 
function cc = cheb(xx,nn)
%xx: root for Chebyshev nodes
%nn: degree of polynomials. 
%cc: cheby poly for each root
cc = real(cos(kron(acos(xx),nn)));
end
