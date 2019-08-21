function set = chebyextrema(i)
if i == 1
    m = 1;
    set(i) = 0;
else
    m = 2^(i-1)+1;
    for j = 1:m
        set(j) = -cos((j-1)*pi/(m-1));
    end
end
end