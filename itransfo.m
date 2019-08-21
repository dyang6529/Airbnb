%from [-1,1] to [min,max] for any x. 
function z = itransfo(x,xmin,xmax)
z = 0.5*(x+1)*(xmax-xmin)+xmin;
end