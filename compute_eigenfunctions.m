% Computes FEM approximations of the N_eigs lowest eigenfunctions of the Julia
% set.

% ------------------------------------------------------------------------
% HOW THE TRINGULATION WORKS:
% c4n contains the coordinates of all nodes
% n4e contains triplets of indices identifying the 3 corners of a triangle. 
%
% Example: n4e(123,:)  = [5794, 5854, 5826]
%
%          c4n(5794,:) = [1.2743, -0.2334]  coordinates of first corner
%          c4n(5854,:) = [1.3065, -0.2346]  coordinates of second corner
%          c4n(5826,:) = [1.2905, -0.2082]  coordinates of third corner
% ------------------------------------------------------------------------

%% Load triangulation:
clear;
load('Julia0.025.mat')
TR = triangulation(n4e,c4n);
Db = freeBoundary(TR);
[nC,d]  = size(c4n);            % number of nodes
nE      = size(n4e,1);          % number of elements
dNodes  = unique(Db);           % Dirichlet boundary
fNodes  = setdiff(1:nC,dNodes); % free nodes
%% Get mass and stiffness matrices and compute eigenfunctions:
% s     = stiffness matrix,
% m     = mass matrix
% vol_T = volumes of elements
% mp_T  = midpoints of elements
N_eigs = 16;
[s,m,vol_T,mp_T] = fe_matrices(c4n,n4e);
S = s(fNodes,fNodes);
M = m(fNodes,fNodes);
[V,D] = eigs(full(S),full(M), N_eigs, 'smallestabs');
% [V, Spectrum, iresult] = sptarn(S,M,0,10);
Spectrum = cast(diag(D),'like',1+1i);
figure('Position',[100,100,1400,300])
plot(Spectrum, '.', 'MarkerSize',15);
title('Spectrum')

%% Tile plots:
figure('Position',[100,100,1400,700])
for j=1:N_eigs
    e = zeros(nC,1);
    e(fNodes) = V(:,j);
    e = e/sqrt(e'*(m*e));
    subplot(4,4,j)
    patch('vertices',c4n,'faces',n4e,'FaceVertexCData',e,'FaceColor','interp','EdgeColor','none')
    axis off
    colorbar
    c = max(abs(e)); 
    caxis([-0.9*c,0.9*c])
    if j==1
        title(['1st Eigenfunction'])
    elseif j==2
        title(['2nd Eigenfunction'])
    elseif j==3
        title(['3rd Eigenfunction'])
    else
        title([num2str(j),'th Eigenfunction'])
    end
end








