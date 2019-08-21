%finding roots for poly. 
function rr = rcheb(nn)
mod = nn-floor(nn/2)*2;
n1 = floor(nn/2);
k = 1:n1;
r1 = cos((2*k-1)*pi/(2*nn));
r1 = [r1 -r1];
if mod == 1
  r1 = [r1 0];
end
rr = real(sort(r1'));
end
