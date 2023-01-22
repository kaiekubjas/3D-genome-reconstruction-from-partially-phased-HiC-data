function [X,Y] = simulate_chromosomes(n,options)

% This function simulates a chromose pair
%
% Each of the chromosomes consists of n beads, and the coordinates of these
% beads are recorded in 3-by-n matrices called X and Y, respectively.
%
% Will will think of the chromose described by X as the maternal
% chromosome, and the chromosome described by Y as the paternal chromosome.
%
% Input:
% n                      desired number of homologous bead pairs
% method (optional)      method of simulation, either 'brownian' (default) or
%                        'brownian_normalized'
% separation (optional)  distance between the first homologous bead pair


arguments
    n
    options.method
    options.separation
    options.plot_it
end

if isfield(options,'method')
    method = options.method;
else
    method = 'brownian';
end

if isfield(options,'separation')
    separation = options.separation;
else
    separation = 2;
end

if isfield(options,'plot_it')
    plot_it = options.plot_it;
else
    plot_it = false;
end

% Simulate the chromosomes
X=zeros(3,n);
X(:,1)=separation/2*[-1;-1;-1];
for i=1:n-1
    step = randn(3,1);
    if strcmp(method,'brownian_normalized')
        step = step/norm(step);
    end
    X(:,i+1)=X(:,i)+step;
end

Y=zeros(3,n);
Y(:,1)=separation/2*[1;1;1];
for i=1:n-1
    step = randn(3,1);
    if strcmp(method,'brownian_normalized')
        step = step/norm(step);
    end
    Y(:,i+1)=Y(:,i)+step;
end

% Plot the chromosome pair
if plot_it
    plot_true
end

end