function [Xstart,Ystart,v_start] = starting_point(Xua,Yua,n,ua_pairs,options)

arguments
    Xua
    Yua
    n
    ua_pairs
    options.strategy
    options.factor
end

% This script generates a starting point that can be used for a local optimization.

% Concretely, it randomly distributes some points around the center of
% mass of Xua and Yua, respectively, using a normal distribution with
% standard estimation estimated from Xua and Yua.

% It relies on data previously computed by the scripts 'simulate.m'
% and 'generate_PA_data.m'.

% Note: This script currently assumes that X and Y correspond to a single
% chromosome each.

if isfield(options,'factor')
    factor = options.factor;
else
    factor = 1;
end

if isfield(options,'strategy')
    strategy = options.strategy;
else
    strategy = 'between_neighbors';
end

ambig_pairs = setdiff(1:n,ua_pairs);
m = length(ambig_pairs);


%% Compute the starting point

Xstart = NaN(3,n);
Xstart(:,ua_pairs) = Xua;

Ystart = NaN(3,n);
Ystart(:,ua_pairs) = Yua;

if strcmp(strategy,'around_origin')
    Xstart(:,ambig_pairs) = factor*std(Xua,0,2).*randn(3,m);
    Ystart(:,ambig_pairs) = factor*std(Xua,0,2).*randn(3,m);
end

if strcmp(strategy,'around_centers_of_mass')
    if size(Xua,2)>1
        Xstart(:,ambig_pairs)=mean(Xua,2) + factor*std(Xua,0,2).*randn(3,m);
    else
        Xstart(:,ambig_pairs)=mean(Xua,2) + randn(3,m);
    end
    
    if size(Yua,2)>1
        Ystart(:,ambig_pairs)=mean(Yua,2) + factor*std(Yua,0,2).*randn(3,m);
    else
        Ystart(:,ambig_pairs)=mean(Yua,2) + randn(3,m);
    end
end

if strcmp(strategy,'between_neighbors')
    
    
    for i = 1:ua_pairs(1)-1
        Xstart(:,i) = Xua(:,1)+(ua_pairs(1)-i)*(Xua(:,1)-Xua(:,2))/(ua_pairs(2)-ua_pairs(1));%+randn(3,1);
        Ystart(:,i) = Yua(:,1)+(ua_pairs(1)-i)*(Yua(:,1)-Yua(:,2))/(ua_pairs(2)-ua_pairs(1));%+randn(3,1);
    end
    
    for i = ua_pairs(end)+1:n
        Xstart(:,i) = Xua(:,end)+(i-ua_pairs(end))*(Xua(:,end)-Xua(:,end-1))/(ua_pairs(end)-ua_pairs(end-1));%+randn(3,1);
        Ystart(:,i) = Yua(:,end)+(i-ua_pairs(end))*(Yua(:,end)-Yua(:,end-1))/(ua_pairs(end)-ua_pairs(end-1));%+randn(3,1);
    end
    
    for k=1:length(ua_pairs)-1
        for i=ua_pairs(k)+1:ua_pairs(k+1)-1
            Xstart(:,i) = Xstart(:,i-1) + (Xua(:,k+1)-Xua(:,k))/(ua_pairs(k+1)-ua_pairs(k));% + (ua_pairs(k+1)-ua_pairs(k))*factor*randn(3,1);
            Ystart(:,i) = Ystart(:,i-1) + (Yua(:,k+1)-Yua(:,k))/(ua_pairs(k+1)-ua_pairs(k));% + (ua_pairs(k+1)-ua_pairs(k))*factor*randn(3,1);
        end
    end
    
    
%     for k=1:length(ua_pairs)-1
%         for i=ua_pairs(k)+1:ua_pairs(k+1)-1
%             Xstart(:,i) = Xstart(:,i) + factor*(Xua(:,k+1)-Xua(:,k))/(ua_pairs(k+1)-ua_pairs(k)).*randn(3,1);
%             Ystart(:,i) = Ystart(:,i) + factor*(Yua(:,k+1)-Yua(:,k))/(ua_pairs(k+1)-ua_pairs(k)).*randn(3,1);
%         end
%     end
    
    % for k=1:length(ua_pairs)-1
    %     for i=ua_pairs(k)+1:ua_pairs(k+1)-1
    %         Xstart(:,i) = Xstart(:,i) + (ua_pairs(k+1)-ua_pairs(k))*factor*randn(3,1);
    %         Ystart(:,i) = Ystart(:,i) + (ua_pairs(k+1)-ua_pairs(k))*factor*randn(3,1);
    %     end
    % end
    
    
    Xstart(:,ambig_pairs) = Xstart(:,ambig_pairs) + factor*std(Xua,0,2).*randn(3,m);
    Ystart(:,ambig_pairs) = Ystart(:,ambig_pairs) + factor*std(Yua,0,2).*randn(3,m);
    
end

%% Collect the coordinates of the ambiguous beads into a single column
% vector of length 6*m
v_start = [reshape(Xstart(:,ambig_pairs),[3*m,1]);reshape(Ystart(:,ambig_pairs),[3*m,1])];

end