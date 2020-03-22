function m = OneDim2TwoDimIndex(I, D)

i = 0;
m = zeros(1, D);
totalSum = 0;

while i < I
    j = 0;
    while m(2) < totalSum
        m(1) = totalSum - j;
        m(2) = j;
        
        j = j+1;
        i = i + 1;
        if i == I
            return
        end
    end
    m = zeros(1, D);
    totalSum = totalSum+1;
end




end

