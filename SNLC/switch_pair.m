function [X,Y] = switch_pair(X,Y,indices)

X_old = X;
Y_old = Y;

for i = indices
    
    X(:,i) = Y_old(:,i);
    Y(:,i) = X_old(:,i);
    
end

end