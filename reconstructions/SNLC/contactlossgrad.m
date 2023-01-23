function gradient = contactlossgrad(v,Xua,Yua,P,ambig_pairs,options)

arguments
    v
    Xua
    Yua
    P
    ambig_pairs
    options.alpha
end

if isfield(options,'alpha')
    alpha = options.alpha;
else
    alpha = -2;
end

n = size(P,2);
ua_pairs = setdiff(1:n,ambig_pairs);
m = length(ambig_pairs);


X = NaN(3,n);
X(:,ua_pairs)=Xua;
X(:,ambig_pairs)=reshape(v(1:3*m),[3,m]);

Y = NaN(3,n);
Y(:,ua_pairs)=Yua;
Y(:,ambig_pairs)=reshape(v(3*m+1:end),[3,m]);


% Compute the gradient

Xgrad = zeros(3,n);
Ygrad = zeros(3,n);

for i=ua_pairs
    for j=ambig_pairs
        if ~isnan(P(i,j)) && ~isinf(P(i,j))
    
            Xgrad(:,j) = Xgrad(:,j) + 2*(P(i,j) - norm(X(:,i)-X(:,j))^alpha - norm(X(:,i)-Y(:,j))^alpha)*...
                (alpha*norm(X(:,i)-X(:,j))^(alpha-2))*(X(:,i)-X(:,j));
            
            Ygrad(:,j) = Ygrad(:,j) + 2*(P(i,j) - norm(X(:,i)-X(:,j))^alpha - norm(X(:,i)-Y(:,j))^alpha)*...
                (alpha*norm(X(:,i)-Y(:,j))^(alpha-2))*(X(:,i)-Y(:,j));
            
            Xgrad(:,i) = Xgrad(:,i) - 2*(P(i,j) - norm(X(:,i)-X(:,j))^alpha - norm(X(:,i)-Y(:,j))^alpha)*...
                alpha*((norm(X(:,i)-X(:,j))^(alpha-2))*(X(:,i)-X(:,j)) + (norm(X(:,i)-Y(:,j))^(alpha-2))*(X(:,i)-Y(:,j)) );
    
        end
        if ~isnan(P(i+n,j)) && ~isinf(P(i+n,j))
            
            Xgrad(:,j) = Xgrad(:,j) + 2*(P(i+n,j) - norm(Y(:,i)-X(:,j))^alpha - norm(Y(:,i)-Y(:,j))^alpha)*...
                (alpha*norm(Y(:,i)-X(:,j))^(alpha-2))*(Y(:,i)-X(:,j));
            
            Ygrad(:,j) = Ygrad(:,j) + 2*(P(i+n,j) - norm(Y(:,i)-X(:,j))^alpha - norm(Y(:,i)-Y(:,j))^alpha)*...
                (alpha*norm(Y(:,i)-Y(:,j))^(alpha-2))*(Y(:,i)-Y(:,j));
            
            Ygrad(:,i) = Ygrad(:,j) - 2*(P(i+n,j) - norm(Y(:,i)-X(:,j))^alpha - norm(Y(:,i)-Y(:,j))^alpha)*...
                alpha*((norm(Y(:,i)-Y(:,j))^(alpha-2))*(Y(:,i)-Y(:,j))+(norm(Y(:,i)-X(:,j))^(alpha-2))*(Y(:,i)-X(:,j)));
            
        end
    end
end

gradient = [Xgrad(:,ambig_pairs),Ygrad(:,ambig_pairs)];
gradient = gradient(:);

end