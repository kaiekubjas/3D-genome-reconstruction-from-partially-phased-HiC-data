function [Xua,Yua,T] = estimate_disambiguated(U,ua_pairs,options)

arguments
    U
    ua_pairs
    options.alpha
    options.optimization_method
end

if isfield(options,'optimization_method')
    optimization_method = options.optimization_method;
else
    optimization_method = 'mdscale';
end

if isfield(options,'alpha')
    alpha = options.alpha;
else
    alpha = -2;
end


% This script estimate the coordinates of the disambiguated beads from a
% matrix U of unambiguous contact frequencies
% These are the z's in the (old version of the) pdf file
% This corresponds to solving a noisy Euclidean distance problem

n = size(U,2)/2;
m = n-length(ua_pairs);

% Collapse U to an 2*(n-m)-by-2*(n-m) matrix
Uua = [U(ua_pairs,ua_pairs),U(ua_pairs,ua_pairs+n);...
    U(ua_pairs+n,ua_pairs),U(ua_pairs+n,ua_pairs+n)];

% Convert contact frequencies to distances
Dua = Uua.^(1/alpha);

if strcmp(optimization_method,'chromsde')
    tic;
    addpath(genpath('./ChromSDE_program2.2/program/'));
    [~,Pest2] = ChromSDE_knownAlpha(Uua,-1/alpha,1);
 
    
    Xua = Pest2(:,1:n-m);
    Yua = Pest2(:,n-m+1:end);
    
    %Compensate for a normalization that happens in ChromSDE
    Uua_new = NaN(2*(n-m));
    for p = 1:(n-m)
        for q = 1:(n-m)
            Uua_new(p,q) = norm(Xua(:,p)-Xua(:,q))^alpha;
            Uua_new(p+n-m,q) = norm(Yua(:,p)-Xua(:,q))^alpha;
            Uua_new(p,q+n-m) = norm(Xua(:,p)-Yua(:,q))^alpha;
            Uua_new(p+n-m,q+n-m) = norm(Yua(:,p)-Yua(:,q))^alpha;
        end
    end
    factor = (mean(Uua(isfinite(Uua)))/mean(Uua_new(isfinite(Uua_new))))^(1/alpha);
    Xua = factor*Xua;
    Yua = factor*Yua;
    T = toc;
    
%     for p = 1:(n-m)
%         for q = 1:(n-m)
%             Uua_new(p,q) = 1/norm(Xua(:,p)-Xua(:,q))^2;
%             Uua_new(p+n-m,q)=1/norm(Yua(:,p)-Xua(:,q))^2;
%             Uua_new(p,q+n-m) =1/norm(Xua(:,p)-Yua(:,q))^2;
%             Uua_new(p+n-m,q+n-m) = 1/norm(Yua(:,p)-Yua(:,q))^2;
%         end
%     end
%     factor = sqrt(mean(Uua(isfinite(Uua)))/mean(Uua_new(isfinite(Uua_new))));
%     Xua = Xua/factor;
%     Yua = Yua/factor;
%     
    rmpath('./ChromSDE_program2.2/program/helperfunctions')
end

if strcmp(optimization_method,'mdscale')
    tic;
    Dua(find(Dua==Inf))=NaN;
    for i = 1:size(Dua,2)
        Dua(i,i)=0;
    end
    
    output = mdscale(Dua,3,'Start','random');
    
    Xua = output(1:n-m,:)';
    Yua = output(n-m+1:end,:)';
    T = toc;
end

end


