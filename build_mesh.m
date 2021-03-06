function [c4n, n4e] = build_mesh(h)
% Computes uniform triangulation of size h for the filled Julia set.
    L = build_lattice(-2-1i, 2+1i, h);
    M = in_Julia(L, 1000, 2);
    J = L(M);    
    nodes = zeros(4,length(J));
    T1 = zeros(length(J(:)), 3);
    T2 = zeros(length(J(:)), 3);
    disp('identified domain')
    disp(['length J: ',num2str(length(J))])
    for i=1:length(J)
        z = J(i);
        corners = [z-h/2-h/2*1i, z+h/2+h/2*1i, z+h/2-h/2*1i, z-h/2+h/2*1i];
        T1(i,:) = [corners(3), corners(2), corners(1)];
        T2(i,:) = [corners(1), corners(2), corners(4)];
        nodes(:,i) = corners;
    end
    nodes = unique(nodes(nodes~=0));
    disp('nodes done')
    c4n = uniquetol([real(nodes), imag(nodes)],h/10,'ByRows',true)
    nodes = (c4n(:,1) + 1i*c4n(:,2)).';
    
    n4e_1 = zeros(length(J),3);
    n4e_2 = zeros(length(J),3);
    tol = h/10;
    parfor i=1:length(J)
        t1 = [find(abs(nodes-T1(i,1))<tol), find(abs(nodes-T1(i,2))<tol), find(abs(nodes-T1(i,3))<tol)];
        t2 = [find(abs(nodes-T2(i,1))<tol), find(abs(nodes-T2(i,2))<tol), find(abs(nodes-T2(i,3))<tol)];
        n4e_1(i,:) = t1;
        n4e_2(i,:) = t2;
    end
    n4e = [n4e_1; n4e_2];
    disp('n4e done')
    
    %% Plot and save results:
%     save(['Julia',num2str(h),'.mat'], 'L', 'h', 'M', 'c4n', 'n4e');
    figure
    patch('vertices',c4n,'faces',n4e,'edgecol','k','LineWidth',0.1,'facecol',[.8,.9,1]);
    drawnow
end


function tf = in_Julia(z, iters, escape_radius)
% Returns "True" if z is in Julia set.
    phi = (1+sqrt(5))/2;
    c = (1-phi);
    f = z.^2+c;
    for n=1:iters
        f = f.^2+c;
    end
    tf = abs(f) < escape_radius;
end
    
    
    
    
    
    
    
    