% Computes the lowest eigenvalue of the pair (s,m) by minimizing 
% Rayleigh's quotient using gradient descent.

clear;
load('Julia0.025.mat')
disp('Computing starting point')

TR = triangulation(n4e,c4n);
triplot(TR)
Db = freeBoundary(TR);
[nC,d]  = size(c4n);            % number of nodes
nE      = size(n4e,1);          % number of elements
dNodes  = unique(Db);           % Dirichlet boundary
fNodes  = setdiff(1:nC,dNodes); % free nodes
[s,m,vol_T,mp_T] = fe_matrices(c4n,n4e);

x = zeros(nC,1);
S = s(fNodes,fNodes);
M = m(fNodes,fNodes);
[V,D] = eigs(S,1,'smallestabs'); % EV of s is good initial guess for GD
x(fNodes) = V;
V = x/sqrt(x'*m*x);
disp('done; starting GD...')

%%
stepsize = 0.1;
DRay = 1;
iters = 100;
Rays = zeros(1,iters);
for i=1:iters
    num   = V'*s*V;
    denom = V'*m*V;
    Ray   = num/denom;
    DRay  = ((s*V) - (m*V)*Ray)/denom;
    V     = V - stepsize*DRay;
    V(dNodes) = 0;
    Rays(i) = Ray;
    i
end
disp('Done!')

Ray
error = abs(Rays(end-1)-Rays(end))/stepsize
plot(1:iters,Rays)
figure
patch('vertices',c4n,'faces',n4e,'FaceVertexCData',V,'FaceColor','interp','EdgeColor','none')
caxis([-1.6,1.6])
colorbar














