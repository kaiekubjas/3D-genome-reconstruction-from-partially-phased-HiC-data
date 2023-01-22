function [rmsd,transform] = comparison(X,Xprime,Y,Yprime,ua_pairs,options)

arguments
    X
    Xprime
    Y
    Yprime
    ua_pairs
    options.title_string,
    options.tag,
    options.prime_tag,
    options.new_tag,
    options.allow_scaling,
    options.mode,
    options.time
    options.only_ua
    options.highlight
    options.shift
end

% This script plots the results of the optimization, together with
% information about how long the optimization took, and the RMSD value.
% It relies on data computed by the scripts 'simulate.m', 'generate_PA_data.m'
% and 'local_optimization.m'.

close all

% Handle defaults

if isfield(options,'title_string')
    title_string = options.title_string;
else
    title_string = 'Comparison';
end

if isfield(options,'tag')
    tag = options.tag;
else
    tag = 'true';
end

if isfield(options,'prime_tag')
    prime_tag = options.prime_tag;
else
    prime_tag = 'new';
end

if isfield(options,'new_tag')
    prime_tag = options.new_tag;
end

if isfield(options,'highlight')
    highlight = options.highlight;
else
    highlight = 'ua';
end


if isfield(options,'allow_scaling')
    allow_scaling = options.allow_scaling;
else
    allow_scaling = true;
end

if isfield(options,'mode')
    mode = options.mode;
else
    mode = 'lines';
end

if isfield(options,'time')
    T = options.time;
else
    T = NaN;
end

if isfield(options,'shift')
    shift = options.shift;
else
    shift = 0;
end


n = size(X,2);

ambig_pairs = setdiff(1:n,ua_pairs);


% Solve the Procrustes problem

if isfield(options,'only_ua') && options.only_ua
    
    
    if allow_scaling
        [~,~,transform] = procrustes([X(:,ua_pairs),Y(:,ua_pairs)]',[Xprime(:,ua_pairs),Yprime(:,ua_pairs)]');
    else
        [~,~,transform] = procrustes([X(:,ua_pairs),Y(:,ua_pairs)]',[Xprime(:,ua_pairs),Yprime(:,ua_pairs)]','scaling',false);
    end
    
    proc_result = transform.b*[Xprime,Yprime]'*transform.T + transform.c(1,:);
    
    
else
    
    if allow_scaling
        [~,proc_result,transform] = procrustes([X,Y]',[Xprime,Yprime]');
    else
        [~,proc_result,transform] = procrustes([X,Y]',[Xprime,Yprime]','scaling',false);
    end
    
end

proc_result = proc_result';
Xproc = proc_result(:,1:n);
Yproc = proc_result(:,n+1:end);

% Compute RMSD
rmsd = norm([X,Y]-[Xproc,Yproc],'fro') / sqrt(2*n);

% Plot the real and estimated data

if strcmp(mode,'lines')
    plot3(X(1,:),X(2,:),X(3,:),'color','b','marker','.','linewidth',0.7,'DisplayName',['X ',tag])
    hold on
    plot3(Xproc(1,:)+shift,Xproc(2,:)+shift,Xproc(3,:)+shift,'color','r','marker','.','linestyle','--','linewidth',0.7,'DisplayName',['X ',prime_tag])
    plot3(Y(1,:),Y(2,:),Y(3,:),'color','#7E2F8E','marker','.','linewidth',0.7,'DisplayName',['Y ',tag])
    plot3(Yproc(1,:)+shift,Yproc(2,:)+shift,Yproc(3,:)+shift,'linestyle','--','color','#EDB120','linewidth',0.7,'marker','.','DisplayName',['Y ',prime_tag])
elseif strcmp(mode,'points')
    plot3(X(1,:),X(2,:),X(3,:),'color','b','marker','.','markersize',12,'linestyle','none','DisplayName',['X ',tag])
    hold on
    plot3(Xproc(1,:)+shift,Xproc(2,:)+shift,Xproc(3,:)+shift,'color','r','marker','.','markersize',7,'linestyle','none','DisplayName',['X ',prime_tag])
    plot3(Y(1,:),Y(2,:),Y(3,:),'color','#7E2F8E','marker','.','markersize',12,'linestyle','none','DisplayName',['Y ',tag])
    plot3(Yproc(1,:)+shift,Yproc(2,:)+shift,Yproc(3,:)+shift,'marker','.','markersize',7,'linestyle','none','color','#EDB120','DisplayName',['Y ',prime_tag])
    alpha(0.1);
end


% Artificially create legend items for start beads and disambiguated beads
plot(NaN,NaN,'*k','DisplayName','Start')
if strcmp(highlight,'ua')
    plot(NaN,NaN,'ok','DisplayName','Unambiguous')
elseif strcmp(highlight,'ambig')
    plot(NaN,NaN,'ok','DisplayName','Ambiguous')
end

% Add legend to the figure
legend('AutoUpdate','off','Location','southoutside','NumColumns',3)

% Mark start points
plot3(X(1,1),X(2,1),X(3,1),'*','color','b','DisplayName','X start')
plot3(Xproc(1,1)+shift,Xproc(2,1)+shift,Xproc(3,1)+shift,'*','color','r','DisplayName','X start estimated')
plot3(Y(1,1),Y(2,1),Y(3,1),'*','color','#7E2F8E','DisplayName','Y start')
plot3(Yproc(1,1)+shift,Yproc(2,1)+shift,Yproc(3,1)+shift,'*','color','#EDB120','DisplayName','Y start estimated')

if strcmp(highlight,'ua')
    % Mark disambiguated beads
    plot3(X(1,ua_pairs),X(2,ua_pairs),X(3,ua_pairs),'o','color','b','DisplayName','X unambiguous')
    plot3(Xproc(1,ua_pairs),Xproc(2,ua_pairs),Xproc(3,ua_pairs),'o','color','r','DisplayName','X unambiguous estimated')
    plot3(Y(1,ua_pairs),Y(2,ua_pairs),Y(3,ua_pairs),'o','color','#7E2F8E','DisplayName','Y unambiguous')
    plot3(Yproc(1,ua_pairs),Yproc(2,ua_pairs),Yproc(3,ua_pairs),'o','color','#EDB120','DisplayName','Y unambiguous estimated')
elseif strcmp(highlight,'ambig')
    plot3(X(1,ambig_pairs),X(2,ambig_pairs),X(3,ambig_pairs),'o','color','b','DisplayName','X ambiguous')
    plot3(Xproc(1,ambig_pairs),Xproc(2,ambig_pairs),Xproc(3,ambig_pairs),'o','color','r','DisplayName','X ambiguous estimated')
    plot3(Y(1,ambig_pairs),Y(2,ambig_pairs),Y(3,ambig_pairs),'o','color','#7E2F8E','DisplayName','Y ambiguous')
    plot3(Yproc(1,ambig_pairs),Yproc(2,ambig_pairs),Yproc(3,ambig_pairs),'o','color','#EDB120','DisplayName','Y ambiguous estimated')
end
    

% Add a title and subtitle
title(title_string)
subtitle({['Indistinguishable pairs: ',num2str(n-length(ua_pairs)),'. ','Distinguishable pairs: ',num2str(length(ua_pairs)),'.'],['Time: ',num2str(T),' sec. ','RMSD: ',num2str(rmsd,'%.2e'),'.'],[]})
pbaspect([1 1 1])
hold off

end