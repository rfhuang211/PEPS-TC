function PEPS=initPEPS(S,r)
% return  PEPS core tensors
% mode: left-down-right-up-physical 
% left2right: i=1:L1
% down2up: j=1:L2
    [L1,L2]=size(S);
    if (L1 <= 1||L2 <= 1)
        error('L1 <= 1 or L2 <= 1.\n');
    end
    
    PEPS=cell(L1,L2);
    
    %% down boundary
    j=1;
    
    i=1;
    PEPS{i,j}=randn(1,1,r,r,S(i,j));
    
    for i = 2:L1-1
        PEPS{i,j}=randn(r,1,r,r,S(i,j));
    end
    
    i=L1;
    PEPS{i,j}=randn(r,1,1,r,S(i,j));
    
    %% up boundary
    j=L2;
    
    i=1;
    PEPS{i,j}=randn(1,r,r,1,S(i,j));
    
    for i = 2:L1-1
        PEPS{i,j}=randn(r,r,r,1,S(i,j));
    end
    
    i=L1;
    j=L2;
    PEPS{i,j}=randn(r,r,1,1,S(i,j));
    
    %% left boundary
    i=1;
    for j = 2:L2-1
        PEPS{i,j}=randn(1,r,r,r,S(i,j));
    end
    
    %% right boundary
    i=L1;
    for j = 2:L2-1
        PEPS{i,j}=randn(r,r,1,r,S(i,j));
    end
    
    %% internal tensors
    for j=2:L2-1
        for i=2:L1-1
            PEPS{i,j}=randn(r,r,r,r,S(i,j));
        end
    end
    
    %% Normalize
    for j=1:L2
        for i=1:L1
            PEPS{i,j}=PEPS{i,j}./max(abs(PEPS{i,j}(:)));
        end
    end

end
