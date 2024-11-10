function tensor=contractPEPS(peps_core)
% combine peps-cores back to tensor
% input must be cell mode core tensors with size L1xL2
% Algorithm:
% S1.contact first line; 
% S2.contact next line; 
% S3;contact the result of S1 and the result of S2; 
% if not the last line, goto S2;

[L1,L2]=size(peps_core);
if (L1 <= 1||L2 <= 1)
    error('L1 <= 1 or L2 <= 1.\n');
end
%% down boundary
%S1
j = 1;
P=peps_core{1,j}; 
dimP = size(P);
for i=1:L1-1
    order = length(dimP);
    index=1:order;
    index(order-2)=order-1;
    index(order-1)=order;
    index(order)=order-2;
    P=permute(P,index);
    L=reshape(P,[numel(P)/dimP(order-2), dimP(order-2)]);
    dimP(order-2) = 1;
    
    dim = size(peps_core{i+1,j});
    order = length(dim);
    R=reshape(peps_core{i+1,j},[dim(1),numel(peps_core{i+1,j})/dim(1)]);
    P=L*R;
    dim(1)=1;
    dimP=[dimP dim];
    assert(numel(P)==prod(dimP));
    P=reshape(P,dimP);
end

%% internal lines
for j = 2:L2
    %S2
    T=peps_core{1,j}; 
    dimT = size(T);
    for i=1:L1-1
        order = length(dimT);
        index=1:order;
        index(order-2)=order-1;
        index(order-1)=order;
        index(order)=order-2;
        T=permute(T,index);
        L=reshape(T,[numel(T)/dimT(order-2), dimT(order-2)]);
        dimT(order-2) = 1;

        dim = size(peps_core{i+1,j});
        order = length(dim);
        R=reshape(peps_core{i+1,j},[dim(1),numel(peps_core{i+1,j})/dim(1)]);
        T=L*R;
        dim(1)=1;
        dimT=[dimT dim];
        assert(numel(T)==prod(dimT));
        T=reshape(T,dimT);
    end
    
    %S3
    indexP = 1:(j-1)*L1*5;
    matPdim2 = 1;
    for k=1:L1
        indexP((j-2)*L1*5 + (k-1)*4 + 1) = (j-2)*L1*5 + (k-1)*5 + 1;
        indexP((j-2)*L1*5 + (k-1)*4 + 2) = (j-2)*L1*5 + (k-1)*5 + 2;
        indexP((j-2)*L1*5 + (k-1)*4 + 3) = (j-2)*L1*5 + (k-1)*5 + 3;
        indexP((j-2)*L1*5 + (k-1)*4 + 4) = (j-2)*L1*5 + (k-1)*5 + 5;
        
        indexP((j-2)*L1*5 + L1*4 + k) = (j-2)*L1*5 + (k-1)*5 + 4;
        
        matPdim2 = matPdim2*dimP((j-2)*L1*5 + (k-1)*5 + 4);
    end
    
    indexT = 1:L1*5;
    matTdim1 = 1;
    for k=1:L1
        indexT(L1 + (k-1)*4 + 1) = (k-1)*5 + 1;
        indexT(L1 + (k-1)*4 + 2) = (k-1)*5 + 3;
        indexT(L1 + (k-1)*4 + 3) = (k-1)*5 + 4;
        indexT(L1 + (k-1)*4 + 4) = (k-1)*5 + 5;
        
        indexT(k) = (k-1)*5 + 2;
        
        matTdim1 = matTdim1*dimT((k-1)*5 + 2);
    end
    assert(matPdim2==matTdim1);
    
    L= permute(P,indexP);
    L= reshape(L,[numel(P)/matPdim2 matPdim2]);
    
    R= permute(T,indexT);
    R= reshape(R,[matTdim1 numel(T)/matPdim2]);
    P=L*R;
    
    for k=1:L1
        dimP((j-2)*L1*5 + (k-1)*5 + 4) = 1;
        dimT((k-1)*5 + 2) = 1;
    end
    dimP=[dimP dimT];
    assert(numel(P)==prod(dimP));
    P=reshape(P,dimP);
end

%reshape P to original tensor size
S=[];
for j=1:L2
    for i = 1:L1
        S=[S size(peps_core{i,j},5)];
    end
end
assert(numel(P)==prod(S));
tensor=reshape(P,S);
end