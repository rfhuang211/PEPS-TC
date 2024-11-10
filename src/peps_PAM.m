function [X, G, Out]=peps_PAM(F, Omega, opts)
    %% variable set
    if (nargin ~= 3)
        error('Not enough input.\n');
    end
    

    Tol = opts.Tol;
    Rank = opts.Rank;
    MaxIter=opts.MaxIter;
    Rho=opts.Rho;

    if opts.Initial_flag == 1
        G=opts.Initial_value;
    else
        G=initPEPS(opts.PEPS_size, Rank);
    end
    
    Out.RSE = [];
    Out.Iter = 0;
    
    X = F; 
    S=size(F);
    
    
    for k = 1:MaxIter
        Xold = X;
        Gold = G;%        size(Gold)
        %alpha = 0.1;
        %% Update G
        [G, GX] = update_G_PAM(Xold,Gold,S,Rank,Rho);
        
        %% Update X
        X = (GX + Rho*Xold)/(1+Rho);
        X(Omega) = F(Omega);
        
        % %% only for diff(true, curr) plot
        % rse=norm(X(:)-opts.Xtrue(:));  
        % Out.RSE = [Out.RSE,rse];

        %% only for obj fun plot
        X_temp=contractPEPS(G);
        rse=0.5*norm(X(:)-X_temp(:)); 
        Out.RSE = [Out.RSE,rse];

        
        %% check the convergence 
        rse=norm(X(:)-Xold(:));   
        rse=rse/norm(Xold(:));       
        
        if mod(k, 50) == 0  ||   k == 1 
            fprintf('PEPS-TC: iter = %d   RSE = %f   \n', k, rse);
        end
        
        if rse < Tol 
            fprintf('****end: iter = %d   RSE = %f   \n', k, rse);
            break;
        end
    end
    Out.Iter = k;
end