function loss = contactloss(v,Xua,Yua,P,ambig_pairs,options)

% Loss function that for a given choice of bead coordinates computes the
% corresponding partially ambiguous contact frequences and compute the
% sum squared difference to known partially ambiguous contact frequences.

% Inputs:
% v                            : vector of coordinates of the ambiguous
%                                beads of length 6*m
% Xua                          : 3 x (n-m) array of coordinates for the
%                                disambiguated maternal beads.
% Yua                          : 3 x (n-m) array of coordinates for the
%                                disambiguated maternal beads.
% P                            : 2*n x n array of partially ambiguous
%                                contact frequences, where P(i,j) is the
%                                contact counts between the ith bead in the
%                                sequence [X,Y] and X(:,j) and Y(:,j).
% ambig_pairs                  : list of indices for ambiguous beads of
%                                length m
% Outputs:
% loss                         : total loss

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


% Number of pairs of beads
n = size(P,2);
ua_pairs = setdiff(1:n,ambig_pairs);
m = length(ambig_pairs);

% Assemble v(1:3*m) and Xua into an 3-by-n matrix X
X = NaN(3,n);
X(:,ua_pairs)=Xua;
X(:,ambig_pairs)=reshape(v(1:3*m),[3,m]);

% Assemble v(3*m+1:end) and Yua into an 3-by-n matrix Y
Y = NaN(3,n);
Y(:,ua_pairs)=Yua;
Y(:,ambig_pairs)=reshape(v(3*m+1:end),[3,m]);

% Compute the loss
loss = 0;

for i=ua_pairs
    for j=ambig_pairs
        if ~isnan(P(i,j)) && ~isinf(P(i,j))
            loss = loss + ( -P(i,j) + norm(X(:,i)-X(:,j))^alpha + norm(X(:,i)-Y(:,j))^alpha )^2;
        end
        if ~isnan(P(i+n,j)) && ~isinf(P(i+n,j))
            loss = loss + ( -P(i+n,j) + norm(Y(:,i)-X(:,j))^alpha + norm(Y(:,i)-Y(:,j))^alpha )^2;
        end
    end
end



end

