function [rmsd,Xproc,Yproc,transformation] = compute_rmsd(X,Xest,Y,Yest,options)

arguments
    X
    Xest
    Y
    Yest
    options.allow_scaling
end

if isfield(options,'allow_scaling')
    allow_scaling = options.allow_scaling;
else
    allow_scaling = true;
end

n = size(X,2);

% Solve the Procrustes problem
if allow_scaling
    [~,proc_result,transformation] = procrustes([X,Y]',[Xest,Yest]');
else
    [~,proc_result,transformation] = procrustes([X,Y]',[Xest,Yest]','scaling',false);
end

proc_result = proc_result';
Xproc = proc_result(:,1:n);
Yproc = proc_result(:,n+1:end);

% Compute RMSD
rmsd = norm([X,Y]-[Xproc,Yproc],'fro') / sqrt(2*n);

end