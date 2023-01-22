function [Xest,Yest,switched_pairs] = unmix_chromosomes(Xest,Yest,ua_pairs)

Xest_old = Xest;
Yest_old = Yest;

Xua = Xest(:,ua_pairs);
Yua = Yest(:,ua_pairs);

n = size(Xest,2);
ambig_pairs = setdiff(1:n,ua_pairs);

tvec = NaN(1,n);
tvec(ua_pairs)=1;

for i=flip(1:ua_pairs(1)-1)
    if norm(Xest(:,i)-Xest(:,i+1))^2 + norm(Yest(:,i)-Yest(:,i+1))^2  > norm(Xest(:,i)-Yest(:,i+1))^2 + norm(Yest(:,i)-Xest(:,i+1))^2
        tvec(i)=-tvec(i+1);
    else
        tvec(i)=tvec(i+1);
    end
end

for idx = 1:length(ua_pairs)-1
    
    ell=0;
    wvec = NaN(1,n);
    
    for i = ua_pairs(idx)+1:ua_pairs(idx+1)-1
        wvec(i) = norm(Xest(:,i)-Xest(:,i-1))^2 + norm(Yest(:,i)-Yest(:,i-1))^2  -( norm(Xest(:,i)-Yest(:,i-1))^2 + norm(Yest(:,i)-Xest(:,i-1))^2 );
        if wvec(i)>0
            tvec(i)=-tvec(i-1);
            ell = ell+1;
        else
            tvec(i)=tvec(i-1);
        end
    end
    
    if mod(ell,2)==1
        [~,k] = min(abs(wvec));
        wvec(k) = -wvec(k);
        for i = ua_pairs(idx)+1:ua_pairs(idx+1)-1
            wvec(i) = norm(Xest(:,i)-Xest(:,i-1))^2 + norm(Yest(:,i)-Yest(:,i-1))^2  -( norm(Xest(:,i)-Yest(:,i-1))^2 + norm(Yest(:,i)-Xest(:,i-1))^2 );
            if wvec(i)>0
                tvec(i)=-tvec(i-1);
                ell = ell+1;
            else
                tvec(i)=tvec(i-1);
            end
            
        end
    end
end



for i = ua_pairs(end)+1:n
    if norm(Xest(:,i)-Xest(:,i-1))^2 + norm(Yest(:,i)-Yest(:,i-1))^2  > norm(Xest(:,i)-Yest(:,i-1))^2 + norm(Yest(:,i)-Xest(:,i-1))^2;
        tvec(i) = -tvec(i-1);
    else
        tvec(i) = tvec(i-1);
    end
    
end


for i=1:n
    if tvec(i)==-1
        [Xest,Yest] = switch_pair(Xest,Yest,i);
    end
end

switched_pairs = find(tvec==-1);

if ~isempty(switched_pairs)
    disp(['The following pairs were switched: ',num2str(switched_pairs),' .'])
end

end