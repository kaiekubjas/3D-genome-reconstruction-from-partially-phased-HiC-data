function [BI_X,BI_Y] = BI_diagram(X,Y,options)

arguments
    X
    Y
    options.alpha
    options.indices
    options.tag
    options.Xname
    options.Yname
end

if ~isfield(options,'alpha')
    alpha = -2;
else
    alpha = options.alpha;
end

if ~isfield(options,'indices')
    indices = 1:size(X,2)
else
    indices = options.indices;
end

if ~isfield(options,'tag')
    tag = [];
else
    tag = options.tag;
end

if ~isfield(options,'Xname')
    Xname = 'X';
else
    Xname = options.Xname;
end

if ~isfield(options,'Yname')
    Yname = 'Y';
else
    Yname = options.Yname;
end


%%
n = size(X,2);


%%


Uest_XX = NaN(n,n);

for i=1:n
    for j=1:n
        Uest_XX(i,j) = 1/norm(X(:,i)-X(:,j))^(-alpha);
    end
end

for i=1:n
   
    Uest_XX(i,i)=0;
    
end


BI_X = NaN(1,n);
numerator = NaN(1,n);
denominator = NaN(1,n);


for h=1:n
    numerator(h) = 0;
    denominator(h) = 0;
    
    for i=1:h
        for j=1:h
            numerator(h) = numerator(h) + Uest_XX(i,j)/h^2;
        end
    end
    
    for i=(h+1):n
        for j=h+1:n
            numerator(h) = numerator(h) + Uest_XX(i,j)/(n-h)^2;
        end
    end
    
    for i=1:h
        for j=(h+1):n
            denominator(h) = denominator(h) + 2*Uest_XX(i,j)/(h*(n-h));
        end
    end
    
    BI_X(h) = numerator(h)/denominator(h);
    
end

%%

Uest_YY = NaN(n,n);

for i=1:n
    for j=1:n
        Uest_YY(i,j) = 1/norm(Y(:,i)-Y(:,j))^(-alpha);
    end
end

for i=1:n
   
    Uest_YY(i,i)=0;
    
end


BI_Y = NaN(1,n);
numerator = NaN(1,n);
denominator = NaN(1,n);


for h=1:n
    numerator(h) = 0;
    denominator(h) = 0;
    
    for i=1:h
        for j=1:h
            numerator(h) = numerator(h) + Uest_YY(i,j)/h^2;
        end
    end
    
    for i=(h+1):n
        for j=h+1:n
            numerator(h) = numerator(h) + Uest_YY(i,j)/(n-h)^2;
        end
    end
    
    for i=1:h
        for j=(h+1):n
            denominator(h) = denominator(h) + 2*Uest_YY(i,j)/(h*(n-h));
        end
    end
    
    BI_Y(h) = numerator(h)/denominator(h);
    
end

%%

plot(indices(1:end),BI_X,'DisplayName',Xname)
hold on
plot(indices(1:end),BI_Y,'DisplayName',Yname)
hold off

xlabel('h')
ylabel('Bipartite index')

title([tag,' (\alpha = ',num2str(alpha),')'])
legend()
  