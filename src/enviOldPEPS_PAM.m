function enOld=enviOldPEPS_PAM(peps_core)
% input must be cell mode core tensors with size L1xL2
% the index with a length equal to 1 will be squeezed
% algorithm:
% %************* calc the L2th line (L1-1:1)
% %************* calc ...
% %************* calc the 1sd line (L1:1)
% %Note: en(k) = en(k+1) + peps_core{k}
    [L1,L2]=size(peps_core);
    if (L1 <= 1 || L2 <= 1)
        error('L1 <= 1 or L2 <= 1.\n');
    end
    enOld=cell(L1,L2);
    legLinks = cell(1,2);

    if L2 ~= 2
        error('L2 ~= 2')
    end

    %%**************** the L2th line
    %%*************** endpoint
    y=L2;
    x=L1;
    enOld{x,y}=[];
    x=L1-1;
    enOld{x,y}=squeeze(peps_core{x+1,y});
    %%*************** inner point
    for x = L1-2:-1:1
        sequence = 1;        
        tensorList = {squeeze(peps_core{x+1,y}), enOld{x+1,y}};
        k = L1-(x+1);
        legLinks{1,1} = [-1 -2 1 -4-(k-1)];
        legLinks{1,2} = [1 -3:-1:-3-(k-1) -5-(k-1):-1:-5-2*(k-1)];
        enOld{x,y} = ncon(tensorList,legLinks,sequence);
    end
    
    %%*************** other lines
    y=L2-1;    
    %if y == L2-1 %(TODO)
        sequence = 1;
        tensorList = {squeeze(peps_core{1,y+1}), enOld{1,y+1}};
        k = L1-1;
        legLinks{1,1} = [-1 1 -3-(k-1)];
        legLinks{1,2} = [1 -2:-1:-2-(k-1) -4-(k-1):-1:-4-2*(k-1)];
        up = ncon(tensorList,legLinks,sequence);
    %end

    % for y = L2-1:-1:2
    % 
    %     %%TODO
    %     error('TODO')
    % 
    % end
    
    
    %%**************** the 1st line
    if  y == 1
        %%*************** endpoint
        x=L1;
        enOld{x,y}=up;
        
        x=L1-1;
        right = squeeze(peps_core{L1,y});
        sequence = 1;
        tensorList = {right, up};
        k = L1-1;
        legLinks{1,1} = [-1 1 -3-(k-1)];
        legLinks{1,2} = [-2:-1:-2-(k-1) 1 -4-(k-1):-1:-4-(k-1)-(L1-1)];
        enOld{x,y} = ncon(tensorList,legLinks,sequence);
        
        %%*************** inner point (two step)
        for x = L1-2:-1:1
            sequence = 1;
            tensorList = {squeeze(peps_core{x+1,y}), right};
            k = L1-(x+1);
            legLinks{1,1} = [-1 1 -2 -4-(k-1)];
            legLinks{1,2} = [1 -3:-1:-3-(k-1) -5-(k-1):-1:-5-2*(k-1)];
            right = ncon(tensorList,legLinks,sequence);

            sequence = 1:L1-x;
            tensorList = {right, up};
            k = L1-x;
            legLinks{1,1} = [-1 sequence -3-(x-1):-1:-3-(x-1)-(k-1)];
            legLinks{1,2} = [-2:-1:-2-(x-1) sequence -4-(x-1)-(k-1):-1:-4-(x-1)-(k-1)-(L1*(L2-y)-1)];
            enOld{x,y} = ncon(tensorList,legLinks,sequence);
        end
    end

    % %%print size of cores
    % for y = L2:-1:1
    %     for x=1:L1
    %         size(enOld{x,y})
    %     end
    % end
end