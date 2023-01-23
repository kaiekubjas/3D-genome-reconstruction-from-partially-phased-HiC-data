using HomotopyContinuation, DynamicPolynomials, LinearAlgebra, Random, MATLAB, NBInclude

mat"""
close all
clear all
"""

# Fix a random number generator seed
Random.seed!(123)
mat"rng(123,'twister')"

# Import functions for the NAG-based estimation
@nbinclude("NAG_functions.ipynb")

# Solve the square system for a random choice of parameters
@var x[1:3];
@var y[1:3];
@var Z[1:3,1:6];
@var C[1:6];
F = System([make_polynomial(x,y,Z[:,i],C[i]) for i=1:6],parameters=vcat(Z[:],C))
Z₀ = randn(ComplexF64,3,6)
C₀ = randn(ComplexF64, 6)
param₀ = vcat(Z₀[:],C₀)
result₀ = solve(F, target_parameters = vcat(Z₀[:],C₀))

# Simulation of chromosomes
n = 60;
m = 30;

mat"""
    n = double($n);
    m = double($m);
    [X,Y]=simulate_chromosomes(n,separation=3,method='brownian_normalized');
    [ua_pairs,ambig_pairs] = partition(n,m,distribution='random');
    $ua_pairs = ua_pairs;
    $ambig_pairs = ambig_pairs;
"""
ua_pairs = map(i->floor(Int,i),ua_pairs);
ambig_pairs = map(i->floor(Int,i),ambig_pairs);

mat"""

eps_percent = 90;
[U,P,A]=generate_contacts(X,Y,ua_pairs,max_error=eps_percent/100);

"""

# Estimation of unambiguous loci
mat"""

[Xua,Yua] = estimate_disambiguated(U,ua_pairs,optimization_method='chromsde');
$Zua = [Xua,Yua];
prows = [ua_pairs,ua_pairs+n];
$Pred = P(prows,ambig_pairs);

"""

# Estimation of ambiguous loci
Xambig, Yambig, Thtpy, attempts, ua_choices, minimal_distances = estimate_ambig_htpy(Pred,Zua,5,
            param0=param₀,result0=result₀,number_of_contacts=20,imag_part_threshold=0.20,max_attempts=100);
Xhtpy = zeros(Float64,3,n);
Xhtpy[:,ambig_pairs]=Xambig;
Xhtpy[:,ua_pairs]=Zua[:,1:n-m];
Yhtpy = zeros(Float64,3,n);
Yhtpy[:,ambig_pairs]=Yambig;
Yhtpy[:,ua_pairs]=Zua[:,n-m+1:end];

display(attempts)
display(maximum(attempts))
display(Thtpy)

# Refine estimations with local optimization
mat"""

[Xesthtpy,Yesthtpy,T,fvalbest,exitflagbest,outputbest,gradbest] = estimate_ambig(...
    Xua,Yua,ua_pairs,P,num_initializations=1,initialization_factor=0.0,Xstart=$Xhtpy,Ystart=$Yhtpy,...
    clustering_step=true);

"""

# Compute RMSD and visualize the estimation together with the true configuration
mat"""

rmsd_esthtpy = comparison(X,Xesthtpy,Y,Yesthtpy,ua_pairs,title_string='Estimation after NGA and local optimization step',...
    new_tag='estimated',allow_scaling=true,mode='lines',only_ua=false,shift=0.01);


view(0,90)
set(gcf, 'Position',  [100, 100, 500, 550])

title('')
subtitle('')

xL = xlim;
yL = ylim;
text(0.99*xL(2),0.99*yL(2),['RMSD = ',num2str(rmsd_esthtpy)],'HorizontalAlignment','right','VerticalAlignment','top')

set(gcf,'renderer','Painters')
print(gcf,'-depsc',['example-noise-',num2str(eps_percent)]);

"""

# # Compute RMSD and visualize the estimation without the local optimization step
mat"""

[Xhtpy,Yhtpy] = unmix_chromosomes($Xhtpy,$Yhtpy,ua_pairs);
Rhtpy = compute_rmsd(X,Xhtpy,Y,Yhtpy,allow_scaling=true);


rmsd_htpy = comparison(X,Xhtpy,Y,Yhtpy,ua_pairs,title_string='Estimation after NGA step',...
    new_tag='estimated',allow_scaling=false,mode='lines',only_ua=true,shift=0.01);


view(0,90)
set(gcf, 'Position',  [100, 100, 500, 550])

title('')
subtitle('')

xL = xlim;
yL = ylim;
text(0.99*xL(2),0.99*yL(2),['RMSD = ',num2str(rmsd_htpy)],'HorizontalAlignment','right','VerticalAlignment','top')

set(gcf,'renderer','Painters')
print(gcf,'-depsc',['NAG-noise-',num2str(eps_percent)]);

"""