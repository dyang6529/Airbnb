
function array = cartesianproduct5(x,y,z,c,d)
[X,Y,Z,C,D] = ndgrid(x,y,z,c,d);
array = [X(:) Y(:) Z(:) C(:) D(:)];
end