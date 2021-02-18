function LW_PLS_OPTIMIZE(X_CAL,Y_CAL,X_test,Y_test,PLS,options)
% calculate residual matrix while varying PHI and num_samp
% haha
tic
% VIP=vip(PLS_model);
% VIP=[];
PHI=[0.1,0.5,1,10,50,100,1000,10000,100000];
num_samp=floor(10:(size(X_CAL,1)-10)/8:size(X_CAL,1));

for i=1:length(PHI)
    for r=1:length(num_samp) 
        
        OUTPUT=LW_PLS(X_CAL,Y_CAL,[],[],X_test,Y_test,PHI(i),PLS,num_samp(r),'score',options);
        y_pred(i,r)=OUTPUT.Y_VAL_PRED(end);
        res(i,r)=OUTPUT.RES;
    end
end
res
[a,b]=min(min(res));
for i=1:length(num_samp)
    [a(i),b(i)]=min(res(:,i));
end
[c,d]=min(a);
% lowest_res=res(b(d),d);
% 
% ind=zeros(length(PHI),length(num_samp));
% for i=1:length(num_samp)
%     [e,f]=find(res(:,i)<=lowest_res*1.1);  %the percentage
%     ind(e,i)=1;
% end
% ind
% [g,h]=find(ind~=0);
% temp_PHI=PHI(g);
% temp_num_samp=num_samp(h);
% 
% for k=1:length(g)
%     DIST(k)=sqrt((g(k)-9).^2+(h(k)-9).^2);
% end
% 
% [Y,I]=min(DIST);
% [x,y]=find(DIST==min(DIST));
% if length(y)>1
%     if h(I)>h(y)
%         I=I;
%     else
%         I=y;
%     end
% else
%     I=I;
% end     

disp(['ideal phi=' num2str(PHI(b(d)))])
disp(['ideal samp num=' num2str(num_samp(d))])
disp(['LWPLS pred value=' num2str(y_pred(b(d),d))])

% disp(['ideal phi=' num2str(PHI(g(I)))])
% disp(['ideal samp num=' num2str(num_samp(h(I)))])
% disp(['LWPLS pred value=' num2str(y_pred(g(I),h(I)))])
disp(['PLS pred value=' num2str(y_pred(end,end))])
disp(['ref value=' num2str(Y_test)])
toc
% pred=y_pred(b(d),d);
% PHI_optimal=PHI(b(d));
% samp=num_samp(d);
% pred=y_pred(g(I),h(I));
% PHI=PHI(g(I))
% samp=num_samp(h(I));

rmsep=abs(ones(length(PHI),length(num_samp))*Y_test-y_pred);

% figure
% imagesc(num_samp,log10(PHI),y_pred)
% xlabel('samp size')
% ylabel('log10(phi)')
% title(['reference value =' num2str(Y_test)])
% colorbar

figure
imagesc(num_samp,log10(PHI),rmsep)
xlabel('samp size')
ylabel('log10(phi)')
title(['lowest error =' num2str(min(min(rmsep)))])
colorbar
% 
figure
imagesc(num_samp,log10(PHI),res)
xlabel('samp size')
ylabel('log10(phi)')
title('spectral residual')
colorbar


