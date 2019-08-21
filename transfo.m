%From [min,max] back to [-1,1] for any x

function z = transfo(x,xmin,xmax)
z = (2*(x-xmin)/(xmax-xmin))-1;
end
