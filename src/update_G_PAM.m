function [Gout, GX]=update_G_PAM(X,G,S,r,rho)
    [L1,L2]=size(G);
    dimNew = 1;
    numelX = numel(X);
    legLinks = cell(1,2);
    Gout = cell(L1,L2);
    
    %%calc old part
    enOld = enviOldPEPS_PAM(G);

    %% update G && accumulate new part
    %%**************** the 1st line
    y=1;
    x=1;
    G1d2 =  enOld{x,y};
    G1d2M = reshape(G1d2, [r*r numel(G1d2)/(r*r)]);
    XM = reshape(X, [S(1) numelX/S(1)]);
    GijM = squeeze(G{x,y});
    GijM = reshape(GijM, [r*r S(1)]);
    GijM = permute(GijM, [2 1]);
    GijM = (XM*G1d2M' + rho*GijM)*pinv(G1d2M*G1d2M'+rho*eye(r*r));
    GijM = permute(GijM, [2 1]);
    Gout{x,y} = reshape(GijM, [1 1 r r S(1)]);
    enNew = squeeze(Gout{x,y});

    %%
    x=2;
    dimNew = dimNew*S(1);
    sequence = 1;        
    tensorList = {enNew, enOld{x,y}};
    legLinks{1,1} = [-1 1 -4];
    legLinks{1,2} = [-2 1 -3 -5:-1:-5-(L1*L2-2-1)];
    G1d2 = ncon(tensorList,legLinks,sequence);
    G1d2M = reshape(G1d2, [r*r*r numel(G1d2)/(r*r*r)]);
    dimOld = numelX/dimNew/S(2);
    XM = reshape(X, [dimNew S(2) dimOld]);
    XM = permute(XM,[2 1 3]);
    XM = reshape(XM, [S(2) dimNew*dimOld]);
    GijM = squeeze(G{x,y});
    GijM = reshape(GijM, [r*r*r S(x)]);
    GijM = permute(GijM, [2 1]);
    GijM = (XM*G1d2M' + rho*GijM)*pinv(G1d2M*G1d2M'+rho*eye(r*r*r));
    GijM = permute(GijM, [2 1]);
    Gout{x,y} = reshape(GijM, [r 1 r r S(x)]);

    %%
    sequence = 1;        
    tensorList = {squeeze(Gout{x, y}), enNew};
    legLinks{1,1} = [1 -1 -3-(x-1-1) -5-2*(x-1-1)];
    legLinks{1,2} = [1 -2:-1:-2-(x-1-1) -4-(x-1-1):-1:-4-2*(x-1-1)];
    enNew = ncon(tensorList,legLinks,sequence);
    
    %%*************** inner point (two step)
    for x=3:L1-1
        %%
        dimNew = dimNew*S(x-1);
        sequence = 1:x-1;        
        tensorList = {enNew, enOld{x,y}};
        legLinks{1,1} = [-1 sequence -4:-1:-4-(x-1-1)];
        legLinks{1,2} = [-2 sequence -3 -5-(x-1-1):-1:-5-(x-1-1)-(L1*L2-x-1)];
        G1d2 = ncon(tensorList,legLinks,sequence);
        G1d2M = reshape(G1d2, [r*r*r numel(G1d2)/(r*r*r)]);
        dimOld = numelX/dimNew/S(x);
        XM = reshape(X, [dimNew S(x) dimOld]);
        XM = permute(XM,[2 1 3]);
        XM = reshape(XM, [S(x) dimNew*dimOld]);
        GijM = squeeze(G{x,y});
        GijM = reshape(GijM, [r*r*r S(x)]);
        GijM = permute(GijM, [2 1]);
        GijM = (XM*G1d2M' + rho*GijM)*pinv(G1d2M*G1d2M'+rho*eye(r*r*r));
        GijM = permute(GijM, [2 1]);
        Gout{x,y} = reshape(GijM, [r 1 r r S(x)]);
        
        %%
        sequence = 1;        
        tensorList = {squeeze(Gout{x, y}), enNew};
        legLinks{1,1} = [1 -1 -3-(x-1-1) -5-2*(x-1-1)];
        legLinks{1,2} = [1 -2:-1:-2-(x-1-1) -4-(x-1-1):-1:-4-2*(x-1-1)];
        enNew = ncon(tensorList,legLinks,sequence);
    end

    %%
    x=L1;
    dimNew = dimNew*S(x-1);
    sequence = 1:x-1;        
    tensorList = {enNew, enOld{x,y}};
    legLinks{1,1} = [-1 sequence -3:-1:-3-(x-1-1)];
    legLinks{1,2} = [sequence -2 -4-(x-1-1):-1:-4-(x-1-1)-(L1*L2-x-1)];
    G1d2 = ncon(tensorList,legLinks,sequence);
    G1d2M = reshape(G1d2, [r*r numel(G1d2)/(r*r)]);
    dimOld = numelX/dimNew/S(x);
    XM = reshape(X, [dimNew S(x) dimOld]);
    XM = permute(XM,[2 1 3]);
    XM = reshape(XM, [S(x) dimNew*dimOld]);
    GijM = squeeze(G{x,y});
    GijM = reshape(GijM, [r*r S(x)]);
    GijM = permute(GijM, [2 1]);    
    GijM = (XM*G1d2M' + rho*GijM)*pinv(G1d2M*G1d2M'+rho*eye(r*r));
    GijM = permute(GijM, [2 1]);
    Gout{x,y} = reshape(GijM, [r 1 1 r S(x)]);

    %%
    sequence = 1;        
    tensorList = {squeeze(Gout{x,y}), enNew};
    legLinks{1,1} = [1 -2-(x-1-1) -4-2*(x-1-1)];
    legLinks{1,2} = [1 -1:-1:-1-(x-1-1) -3-(x-1-1):-1:-3-2*(x-1-1)];
    enNew = ncon(tensorList,legLinks,sequence);
    
    %%**************** other lines
    % for y=2:L2-1
    %     %%TODO
    %     error('TODO')
    % end
    
    %%**************** the last line
    y=2;
    if y == L2
        %%only for 2 lines
        %%enNew for other lines      
        x=1;
        dimNew = dimNew*S(L1);
        sequence = 1:L1-1;        
        tensorList = {enNew, enOld{x,y}};
        legLinks{1,1} = [-1 sequence -3:-1:-3-(L1-1)];
        legLinks{1,2} = [-2 sequence -4-(L1-1):-1:-4-(L1-1)-(L1-1-1)];
        G1d2 = ncon(tensorList,legLinks,sequence);
        G1d2M = reshape(G1d2, [r*r numel(G1d2)/(r*r)]);
        dimOld = numelX/dimNew/S(L1+x);
        XM = reshape(X, [dimNew S(L1+x) dimOld]);
        XM = permute(XM,[2 1 3]);
        XM = reshape(XM, [S(L1+x) dimNew*dimOld]);
        GijM = squeeze(G{x,y});
        GijM = reshape(GijM, [r*r S(L1+x)]);
        GijM = permute(GijM, [2 1]);
        GijM = (XM*G1d2M' + rho*GijM)*pinv(G1d2M*G1d2M'+rho*eye(r*r));
        GijM = permute(GijM, [2 1]);
        Gout{x,y} = reshape(GijM, [1 r r 1 S(L1+x)]);
        
        left = squeeze(Gout{x,y});
        %%*************** inner point (two step)
        for x=2:L1-1
            %%
            legLinks = cell(1,3);
            dimNew = dimNew*S(L1+x-1);
            sequence = 1:L1-1;        
            tensorList = {left, enNew, enOld{x,y}};
            legLinks{1,1} = [1:x-1 -1 -5-(L1-1):-1:-5-(L1-1)-(x-1-1)];
            legLinks{1,2} = [1:x-1 -2 x:L1-1 -4:-1:-4-(L1-1)];
            legLinks{1,3} = [-3 x:L1-1 -6-(L1-1)-(x-1-1):-1:-6-(L1-1)-(x-1-1)-(L1-x-1)];
            G1d2 = ncon(tensorList,legLinks,sequence);
             
            G1d2M = reshape(G1d2, [r*r*r numel(G1d2)/(r*r*r)]);
            dimOld = numelX/dimNew/S(L1+x);
            XM = reshape(X, [dimNew S(L1+x) dimOld]);
            XM = permute(XM,[2 1 3]);
            XM = reshape(XM, [S(L1+x) dimNew*dimOld]);
            GijM = squeeze(G{x,y});
            GijM = reshape(GijM, [r*r*r S(L1+x)]);
            GijM = permute(GijM, [2 1]);
            GijM = (XM*G1d2M' + rho*GijM)*pinv(G1d2M*G1d2M'+rho*eye(r*r*r));
            GijM = permute(GijM, [2 1]);
            Gout{x,y} = reshape(GijM, [r r r 1 S(L1+x)]);

            %% update left
            tensorList = {left, squeeze(Gout{x,y})};
            legLinks = cell(1,2);
            sequence=1;
            legLinks{1,1} = [-1:-1:-(x-1) 1 -(x-1)-3:-1:-(x-1)-3-(x-1-1)];
            legLinks{1,2} = [1 -(x-1)-1 -(x-1)-2 -(x-1)-4-(x-1-1)];
            left = ncon(tensorList,legLinks,sequence);
        end

        x = L1;
        legLinks = cell(1,2);
        dimNew = dimNew*S(L1+x-1);
        sequence = 1:L1-1;        
        tensorList = {left, enNew};
        legLinks{1,1} = [1:x-1 -1 -4-(L1-1):-1:-4-(L1-1)-(x-1-1)];
        legLinks{1,2} = [1:x-1 -2 -3:-1:-3-(L1-1)];
        G1d2 = ncon(tensorList,legLinks,sequence);
         
        G1d2M = reshape(G1d2, [r*r numel(G1d2)/(r*r)]);
        XM = reshape(X, [dimNew S(L1+x)]);
        XM = permute(XM,[2 1]);
        XM = reshape(XM, [S(L1+x) dimNew]);
        GijM = squeeze(G{x,y});
        GijM = reshape(GijM, [r*r S(L1+x)]);
        GijM = permute(GijM, [2 1]);
        GijM = (XM*G1d2M' + rho*GijM)*pinv(G1d2M*G1d2M'+rho*eye(r*r));
        GijM = permute(GijM, [2 1]);
        Gout{x,y} = reshape(GijM, [r r 1 1 S(L1+x)]);

    else
        %%TODO
        error('TODO')
    end


    %% calc GX for updating X (contract PEPS)
    legLinks = cell(1,3);
    sequence=1:L1+1;
    tensorList = {left, enNew, squeeze(Gout{L1,L2})};
    legLinks{1,1} = [1:L1-1 L1 -L1-1:-1:-L1-1-(L1-1-1)];
    legLinks{1,2} = [1:L1-1 L1+1 -1:-1:-L1];
    legLinks{1,3} = [L1 L1+1 -L1*L2];
    GX = ncon(tensorList,legLinks,sequence);

    % %%test 
    % temp = contractPEPS(Gout);
    % max(temp(:)-GX(:))
end