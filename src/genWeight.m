% generate tensor W by size and missing rate
function W=genWeight(S,mr,id)
    W = ones(S); 
    if id == 1
        for k = 1:S(3)
            Omega = randperm(S(1)); 
            Omega = Omega(1:floor(mr*S(1)));
            W(Omega,:,k) = 0;
        end
    elseif id == 2
        block=[16 16];
        num=[floor(S(1)/block(1)) floor(S(2)/block(2))];
        for k = 1:S(3)
            Omega = randperm(num(1)*num(2)); 
            Omega = Omega(1:floor(mr*num(1)*num(2)));
            for j = 1:numel(Omega)
                temp=floor(Omega(j)/num(1))+1;
                y=int16(temp);
                temp=mod(Omega(j),num(1))+1;
                x=int16(temp);
                W((x-1)*block(1)+1:min(x*block(1),S(1)), (y-1)*block(2)+1:min(y*block(2),S(2)), k) = 0;
            end
        end
    else
        Omega = randperm(prod(S)); 
        Omega = Omega(1:floor(mr*prod(S)));
        W(Omega) = 0;
    end

end