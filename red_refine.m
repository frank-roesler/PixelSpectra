function [c4nNew,n4eNew,DbNew,NbNew,P0,P1] = red_refine(c4n,n4e,Db,Nb)
% Refines triangulation by decomposing each element into 4 smaller ones
% This function is taken from [Bartels S. Numerical approximation of partial 
% differential equations. Springer; 2016]
[nC,d] = size(c4n); nE = size(n4e,1); 
nDb = size(Db,1); nNb = size(Nb,1); 
K = sparse(1:nC,1:nC,1:nC);
c4nNew = c4n; n4eNew = zeros(nE*2^d,d+1);
DbNew = zeros(nDb*2^(d-1),d); NbNew = zeros(nNb*2^(d-1),d);
P0 = sparse(2^d*nE,nE);
nr_nodes = nC;
tic
for j = 1:nE
    if mod(j,10000)==0
        disp(num2str(100*j/nE))
        toc
        tic
    end
    for k = 1:d+1
        for m = 1:d+1
            if K(n4e(j,k),n4e(j,m))==0
                nr_nodes = nr_nodes+1;
                K(n4e(j,k),n4e(j,m)) = nr_nodes;
                K(n4e(j,m),n4e(j,k)) = nr_nodes;
                I1(2*(nr_nodes-nC-1)+(1:2)) = [nr_nodes,nr_nodes];
                I2(2*(nr_nodes-nC-1)+(1:2)) = [n4e(j,k),n4e(j,m)];
                EE(2*(nr_nodes-nC-1)+(1:2)) = [1 1]/2;
                c4nNew(nr_nodes,:) = (c4n(n4e(j,k),:)+c4n(n4e(j,m),:))/2;
            end
        end
    end
    nodes = K(n4e(j,:),n4e(j,:));
    if d == 1
        n4eNew(2*(j-1)+(1:2),:) = ...
            [nodes(1,1),nodes(1,2);nodes(1,2),nodes(2,2)];
    elseif d == 2
        n4eNew(4*(j-1)+(1:4),:) = ...                                       
            [nodes(1,1),nodes(1,2),nodes(1,3);...                           
            nodes(1,2),nodes(2,3),nodes(1,3);...
            nodes(1,2),nodes(2,2),nodes(2,3);...
            nodes(1,3),nodes(2,3),nodes(3,3)];
    elseif d == 3
        n4eNew(8*(j-1)+(1:8),:) = ...
            [nodes(1,1),nodes(1,2),nodes(1,3),nodes(1,4);...
            nodes(1,2),nodes(2,2),nodes(2,3),nodes(2,4);...
            nodes(1,3),nodes(2,3),nodes(3,3),nodes(3,4);...
            nodes(1,4),nodes(2,4),nodes(3,4),nodes(4,4);...
            nodes(1,2),nodes(1,3),nodes(1,4),nodes(2,4);...
            nodes(2,3),nodes(1,3),nodes(1,2),nodes(2,4);...
            nodes(1,3),nodes(1,4),nodes(2,4),nodes(3,4);...
            nodes(1,3),nodes(2,4),nodes(2,3),nodes(3,4)];
    end
    P0(2^d*(j-1)+(1:2^d),j) = ones(2^d,1);
end
I1 = [I1,1:nC]; I2 = [I2,1:nC]; EE = [EE,ones(1,nC)];
P1 = sparse(I1,I2,EE,nr_nodes,nC);
for j = 1:nDb
    nodes = K(Db(j,:),Db(j,:));
    if d == 1
        DbNew = Db;
    elseif d == 2
        DbNew(2*(j-1)+(1:2),:) = ...
            [nodes(1,1),nodes(1,2);...
            nodes(1,2),nodes(2,2)];
    elseif d == 3
        DbNew(4*(j-1)+(1:4),:) = ...
            [nodes(1,1),nodes(1,2),nodes(1,3);...
            nodes(1,2),nodes(2,2),nodes(2,3);...
            nodes(1,2),nodes(2,3),nodes(1,3);...
            nodes(1,3),nodes(2,3),nodes(3,3)];
    end
end
for j = 1:nNb
    nodes = K(Nb(j,:),Nb(j,:));
    if d == 1
        NbNew = Nb;
    elseif d == 2
        NbNew(2*(j-1)+(1:2),:) = ...
            [nodes(1,1),nodes(1,2);...
            nodes(1,2),nodes(2,2)];
    elseif d == 3
        NbNew(4*(j-1)+(1:4),:) = ...
            [nodes(1,1),nodes(1,2),nodes(1,3);...
            nodes(1,2),nodes(2,2),nodes(2,3);...
            nodes(1,2),nodes(2,3),nodes(1,3);...
            nodes(1,3),nodes(2,3),nodes(3,3)];
    end
end
end



