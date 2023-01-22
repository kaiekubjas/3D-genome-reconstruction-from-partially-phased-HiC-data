function [Unew,Pnew,Anew] = preprocess_contacts(U,P,A,ambig_pairs)

% Function that ensures that adapts the contact count matrices to the
% simplifying assumption that the loci can be partitioned into ambiguous
% and unambiguous loci.

n = size(A,1);

ua_pairs = setdiff(1:n,ambig_pairs);

Unew = U;
Pnew = P;
Anew = A;

for i=ua_pairs
    for j=ua_pairs
        
        % Split P(i,j) into contributions to U(i,j) and U(i,j+n).
        % (Since U is symmetric, we similarly split P(i,j) into a contribution to 
        % U(j,i) and U(j+n,i).)
        
        [Unew(i,j),Unew(i,j+n)] = unpack_vector(...
            [Unew(i,j),Unew(i,j+n)] + ...
            P(i,j)*ratio([U(i,j),U(i,j+n)]) );
        
        [Unew(j,i),Unew(j+n,i)] = unpack_vector(...
            [Unew(j,i),Unew(j+n,i)] + ...
            P(i,j)*ratio([U(j,i),U(j+n,i)]) );
        
        Pnew(i,j) = 0;
        
        % We split P(i+n,j) into contributions to U(i+n,j) and U(i+n,j+n), 
        % as well as U(j,i+n) and U(j+n,i+n).
        
        [Unew(i+n,j),Unew(i+n,j+n)] = unpack_vector(...
            [Unew(i+n,j),Unew(i+n,j+n)] + ...
            P(i+n,j)*ratio([U(i+n,j),U(i+n,j+n)]) );
        
        [Unew(j,i+n),Unew(j+n,i+n)] = unpack_vector(...
            [Unew(j,i+n),Unew(j+n,i+n)]+ ...
            P(i+n,j)*ratio([U(j,i+n),U(j+n,i+n)]) );
        
        Pnew(i+n,j) = 0;
        
        % Split A(i,j) into contributions to U(i,j),U(i,j+n),U(i+n,j),U(i+n,j+n).
        % (Since A is symmetric, we will split A(j,i) analogously in a
        % another iteration, which ensures that U remains symmetric.)
        
        [Unew(i,j),Unew(i,j+n),Unew(i+n,j),Unew(i+n,j+n)] = unpack_vector(...
            [Unew(i,j),Unew(i,j+n),Unew(i+n,j),Unew(i+n,j+n)] + ...
            A(i,j)*ratio([U(i,j),U(i,j+n),U(i+n,j),U(i+n,j+n)]) );
        
        Anew(i,j) = 0;
        
    end    
end

for i=ambig_pairs
    for j=ambig_pairs
        
        % Combine U(i,j), U(i+n,j), U(i,j+n) and U(i+n,j+n) into a 
        % contribution to Anew(i,j)
        
        Anew(i,j) = Anew(i,j) + U(i,j) + U(i+n,j) + U(i,j+n) + U(i+n,j+n);
        Unew(i,j) = 0;
        Unew(i+n,j) = 0;
        Unew(i,j+n) = 0;
        Unew(i+n,j+n) = 0;
        
        % Combine P(i,j), P(i+n,j), P(j,i) and P(j+n,i) into a contribution
        % to Anew(i,j). (In another iteration the same four terms will combine
        % into a contribution to Anew(j,i), which ensures that Anew remains
        % symmetric.)
        
        Anew(i,j) = Anew(i,j) + P(i,j) + P(i+n,j) + P(j,i) + P(j+n,i);
        Pnew(i,j) = 0;
        Pnew(i+n,j) = 0;
        Pnew(j,i) = 0;
        Pnew(j+n,i) = 0;
        
    end
end

for i=ua_pairs
    for j=ambig_pairs
        
        % Combine U(i,j) and U(i,j+n) into a contribution to Pnew(i,j)
        Pnew(i,j) = Pnew(i,j) + U(i,j) + U(i,j+n);
        
        Unew(i,j) = 0;
        Unew(i,j+n) = 0;
        
        Unew(j,i) = 0;
        Unew(j+n,i) = 0;
        
        % Combine U(i+n,j) and U(i+n,j+n) into a contribution to Pnew(i+n,j)
        Pnew(i+n,j) = Pnew(i+n,j) + U(i+n,j) + U(i+n,j+n);
        
        Unew(i+n,j) = 0;
        Unew(i+n,j+n) = 0;
        
        Unew(j,i+n) = 0;
        Unew(j+n,i+n) = 0;
        
        % Split A(i,j) into contributions to Pnew(i,j) and Pnew(i+n,j)
        [Pnew(i,j),Pnew(i+n,j)] = unpack_vector( [Pnew(i,j),Pnew(i+n,j)] + A(i,j)*ratio([P(i,j),P(i+n,j)]) );
        
        Anew(i,j) = 0;
        Anew(j,i) = 0;
        
        % Split P(j,i) and P(j+n,i) into contributions to Pnew(i,j) and Pnew(i+n,j)
        [Pnew(i,j),Pnew(i+n,j)] = unpack_vector( [Pnew(i,j),Pnew(i+n,j)] + P(j,i)*ratio([U(j,i),U(j,i+n)])  );     
        [Pnew(i,j),Pnew(i+n,j)] = unpack_vector( [Pnew(i,j),Pnew(i+n,j)] + P(j+n,i)*ratio([U(j+n,i),U(j+n,i+n)])  );     
        
        Pnew(j,i) = 0;
        Pnew(j+n,i) = 0;
        
    end
end


end