function [coeff,score,pc_percents,all_percents,data_new,rmse_all]=gps_vbpca(data_raw,k_num)


opts = struct( 'maxiters', 5000,...
               'algorithm', 'vb',...  %% ppca, map, vb 
               'uniquesv', 0,...
               'rmsstop',[1000 1e-8 1e-7 ],...
               'cfstop', [ 1000 0 0 ],...
               'minangle', 1e-14 );
          
[ coeff,score, mu, pcvarp, cv, hp, lc ] = pca_full(data_raw',k_num, opts ); 

score=score';mu=mu';
% compare raw and filtered time series
[no_epoch,no_site]=size(data_raw);
data_new = score*coeff'+ repmat(mu,no_epoch,1);

Var=zeros(no_site,1);
for i=1:1:no_site
    ok=~isnan(data_raw(:,i));
    d=data_raw(ok,i);
    dhat=data_new(ok,i);
    r=d-dhat;
   Var(i,1)=(d'*d-r'*r)/(d'*d);
end
all_percents=mean(Var);
rmse_all=Var;
disp(' Variacne reduction (%)'); mean(Var)

data1=score*coeff';
for i=1:k_num
    data2=score(:,i)*coeff(:,i)';
    for j=1:1:no_site
        d=data1(:,j);
        dhat=data2(:,j);
        r=d-dhat;
        Var1(i,j)=(d'*d-r'*r)/(d'*d);
    end
    pc_percents(i)=mean(Var1(i,:));
end
