# PixelSpectra
This Matlab package computes the Dirichlet eigenvalues and eigenfunctions of domains defined through a pixelation. The computation is done via Finite Element approximation on a uniform triangulation, according to a method studied in https://alexeistepa.github.io/rough_eig_preprint.pdf.  
The main building blocks of the routine are
* ```build_mesh(h)``` computes the triangulation (of a filled Julia set) for a given mesh size ```h```;
* ```refine_mesh(c4n,n4e)``` computes a refined mesh from a given triangulation ```(c4n,n4e)``` by splitting every triangle into 4 smaller ones;
* ```compute_eigenfunctions.m``` computes the lowest eigenvalues and eigenfunctions for a given triangulation;
* ```compute_lowest_eig.m``` computes a high-accuracy approximation of the lowest eigenvalue using the Rayleigh-Ritz method & gradient descent.
A pre-computed triangulation for a filled Julia set is included to make the code executable out-of-the-box, and the code in ```build_mesh``` can easily be adapted to arbitrary domains.

The Finite Element code used in this package is partially adapted from [Bartels S. *Numerical approximation of partial differential equations.* Springer; 2016]

Any comments or queries are welcome at https://frank-roesler.github.io/contact/
