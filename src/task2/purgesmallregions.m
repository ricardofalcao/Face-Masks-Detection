function Out = purgesmallregions(BW)
    [Reg, N] = bwlabel(BW);
    Area = zeros(N,1);

    for i = 1:N
       Area(i) = nnz(Reg == i);
    end

    RegAvg = mean(Area(:));
    Out = zeros(size(BW));

    for i = 1:N
       if Area(i) >= RegAvg 
           Zone = (Reg == i);
           Out = Out | Zone;
       end
    end
end