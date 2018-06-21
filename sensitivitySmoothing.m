function [senSm] = sensitivitySmoothing(sens, wSm, Cc)

    Nsmooth = floor(length(wSm)/2);
    
    wrapN = @(index, N) (1 + mod(index-1, N));
    Cc2 = Cc(1:end-1);
    Nc2 = length(Cc2);
    senSm = sens;
    
    for i = 1:length(Cc2)
        k = Cc2(i);
        senSm(k) = wSm(1)*sens(k);
        for j = 1:Nsmooth
            senSm(k) = senSm(k) + wSm(1+j) * ( senSm( Cc2( wrapN(i+j, Nc2) ) ) + ...
                senSm( Cc2( wrapN(i-j, Nc2) ) ) );
        end
    end
    
end