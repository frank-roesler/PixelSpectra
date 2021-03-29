function [c4n,n4e,Db,Nb] = refine_mesh(path)
% Returns a refined mesh of half the mesh size as the original one.
    load(path);
    TR = triangulation(n4e,c4n);
    Db = freeBoundary(TR);
    Nb = [];
    figure
    triplot(TR)
    title('Before:')

    [c4n,n4e,Db,Nb,P0,P1] = red_refine(c4n,n4e,Db,Nb);
    TR = triangulation(n4e,c4n);
    figure
    triplot(TR)
    title('After:')
    
    save([path(1:end-4),'_refined.mat'])
end
