function [Xbest,Ybest,T,fvalbest,exitflagbest,outputbest,gradbest] = estimate_ambig(Xua,Yua,ua_pairs,P,options)

arguments
    Xua
    Yua
    ua_pairs
    P
    options.num_initializations
    options.initialization_method
    options.Xstart
    options.Ystart
    options.initialization_factor
    options.clustering_step
    options.alpha
end

% This function estimates the coordinates of the ambiguous beads

% The coordinates of the unknown beads are represented by a single vector v
% of length 6*m. The star values of these coordinates are stored in v_start,
% and the optimum found in the optimization is stored in v_optimum


if isfield(options,'num_initializations')
    num_initializations = options.num_initializations;
else
    num_initializations = 100;
end

if isfield(options,'initialization_method')
    initialization_method = options.initialization_method;
else
    initialization_method = 'between_neighbors';
end


if isfield(options,'initialization_factor')
    initialization_factor = options.initialization_factor;
else
    initialization_factor = 2;
end

if isfield(options,'clustering_step')
    clustering_step = options.clustering_step;
else
    clustering_step = true;
end

if isfield(options,'alpha')
    alpha = options.alpha;
else
    alpha = -2;
end

n = size(P,2);
ambig_pairs = setdiff(1:n,ua_pairs);
m = length(ambig_pairs);



tic;

fval = NaN(1,num_initializations);

for i = 1:num_initializations
    if num_initializations~=1
        disp(['Initialization: ',num2str(i),' of ',num2str(num_initializations),'.'])
    end
    
    if isfield(options,'Xstart') && isfield(options,'Ystart')
        Xstart = options.Xstart + initialization_factor*std(Xua,0,2).*randn(3,n);
        Ystart = options.Ystart + initialization_factor*std(Yua,0,2).*randn(3,n);
    else
        [Xstart,Ystart] = starting_point(Xua,Yua,n,ua_pairs,...
                    strategy=initialization_method,factor=initialization_factor);
    end

        v_start = [reshape(Xstart(:,ambig_pairs),[3*m,1]);reshape(Ystart(:,ambig_pairs),[3*m,1])];

    optimization_options =  optimoptions('fminunc', ...
        'SpecifyObjectiveGradient', true, ...
        'Algorithm','quasi-newton', ...
        'CheckGradient',false,...
        'Display','final',...
        'MaxIterations',1e10,...
        'MaxFunctionEvaluations',1e10,...
        'FunctionTolerance',1e-7,...
        'StepTolerance',1e-10);
    
    [v_optimum,fval(i),exitflag,output,grad]= fminunc({...
        @(v) contactloss(v,Xua,Yua,P,ambig_pairs,alpha=alpha),...
        @(v) contactlossgrad(v,Xua,Yua,P,ambig_pairs,alpha=alpha)},...
        v_start,optimization_options);


    % Combine the coordinates of the ambiguous and unambiguous beads into
    % new arrays Xnew and Ynew

    Xest = NaN(3,n);
    Xest(:,ua_pairs) = Xua;
    Xest(:,ambig_pairs)=reshape(v_optimum(1:3*m),[3,m]);

    Yest = NaN(3,n);
    Yest(:,ua_pairs) = Yua;
    Yest(:,ambig_pairs)=reshape(v_optimum(3*m+1:6*m),[3,m]);
    
    if i==1 || fval(i) == min(fval)
        Xbest = Xest;
        Ybest = Yest;
        fvalbest=fval(i);
        exitflagbest=exitflag;
        outputbest=output;
        gradbest=grad;
    end
    
end

T = toc;

if clustering_step
    [Xbest,Ybest,switched_pairs] = unmix_chromosomes(Xest,Yest,ua_pairs);
end

end