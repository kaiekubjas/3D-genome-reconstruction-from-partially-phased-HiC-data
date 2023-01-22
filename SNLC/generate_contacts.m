function [U,P,A] = generate_PA_data(X,Y,ua_pairs,options)

arguments
    X
    Y
    ua_pairs
    options.alpha
    options.max_error
end

% This function computes ambiguous, partially ambiguous and unambiguous
% contact counts according to our model, for a given alpha, and 
% then multiplies each count by (1+eps) for eps drawn with uniform probability
% form the interval [-max_error,max_error].

% Input:
% X            3-by-n matrix of coordinates of the chromosome X
% Y            3-by-n matrix of coordinates of the chromosome Y
% ua_pairs     length-(n-m) list of the indicies of unambiguous loci
% alpha        the exponent in the distance-based model (defaults to -2)
% max_error    positive number determining radius of noise interval

if isfield(options,'max_error')
    max_error = options.max_error;
else
    max_error = 0;
end

if isfield(options,'alpha')
    alpha = options.alpha;
else
    alpha = -2;
end

n = size(X,2);
ambig_pairs = setdiff(1:n,ua_pairs);
m = length(ambig_pairs);


%% Compute contact frequences based on the theoretical model, and add noise

% Compute a matrix of fully ambiguous contact frquences
Atrue = zeros(n,n);
A = zeros(n,n);
for i = ambig_pairs
    for j = intersect(i+1:n,ambig_pairs)
        Atrue(i,j) = norm(X(:,i)-X(:,j))^alpha + norm(X(:,i)-Y(:,j))^alpha + norm(Y(:,i)-X(:,j))^alpha + norm(Y(:,i)-Y(:,j))^alpha;
        Atrue(j,i) = Atrue(i,j);
        A(i,j) = Atrue(i,j)*((1-max_error) + (1+max_error-(1-max_error)) * rand());
        A(j,i) = A(i,j);
    end
end

% Compute a matrix of partially ambiguous contact frequences
Ptrue = zeros(2*n,n);
P = zeros(2*n,n);
for i = ua_pairs
    for j = ambig_pairs
        Ptrue(i,j) = norm(X(:,i)-X(:,j))^alpha + norm(X(:,i)-Y(:,j))^alpha;
        P(i,j) = Ptrue(i,j)*((1-max_error) + (1+max_error-(1-max_error)) * rand());
        
        Ptrue(i+n,j) = norm(Y(:,i)-X(:,j))^alpha + norm(Y(:,i)-Y(:,j))^alpha;
        P(i+n,j) = Ptrue(i+n,j)*((1-max_error) + (1+max_error-(1-max_error)) * rand());
    end
end

% Compute a matrix of unambiguous contact frequences
Utrue = zeros(2*n,2*n);
U = zeros(2*n,2*n);
for i = ua_pairs
    for j = intersect(ua_pairs,i:n)
        Utrue(i,j) = norm(X(:,i)-X(:,j))^alpha;
        Utrue(j,i) = Utrue(i,j);
        U(i,j) = Utrue(i,j)*((1-max_error) + (1+max_error-(1-max_error)) * rand());
        U(j,i) = U(i,j);
        
        Utrue(i+n,j) = norm(Y(:,i)-X(:,j))^alpha;
        Utrue(j,i+n) = Utrue(i+n,j);
        U(i+n,j) = Utrue(i+n,j)*((1-max_error) + (1+max_error-(1-max_error)) * rand());
        U(j,i+n) = U(i+n,j);
        
        Utrue(i,j+n) = norm(X(:,i)-Y(:,j))^alpha;
        Utrue(j+n,i) = Utrue(i,j+n);
        U(i,j+n) = Utrue(i,j+n)*((1-max_error) + (1+max_error-(1-max_error)) * rand());
        U(j+n,i) = U(i,j+n);
        
        Utrue(i+n,j+n) = norm(Y(:,i)-Y(:,j))^alpha;
        Utrue(j+n,i+n) = Utrue(i+n,j+n);
        U(i+n,j+n) = Utrue(i+n,j+n)*((1-max_error) + (1+max_error-(1-max_error)) * rand());
        U(j+n,i+n) = U(i+n,j+n);
    end
end



