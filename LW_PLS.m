function OUTPUT=LW_PLS(X_CAL,Y_CAL,X_aug,Y_aug,X_VAL,Y_VAL,phi,PLS,num_samp,FLAG,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LW-PLS per International Journal of Pharmaceutics 421 (2011) 269– 274 
% with following modifications
%       1. using either vip or score distance to replace the 
%            original distance calculation 
%       2. adding prediction at each LV up to the LV defined in PLS via 
%           conducting svd across latent variable one at a time.
%
% Inputs: 
%   X_CAL,Y_CAL,X_VAL and Y_VAL: spectra and reference values of
%   calibration and test data
%   phi: the parameter
%   num_samp: the number of local samples used for local PLS model
%   PLS: original original PLS model between X_CAL and Y_CAL, without
%   adding X_aug and Y_aug. Be careful if there are samples removed in PLS
%   struct.
%   FLAG: either "vip" or "score" used to weigh distance between new samples
%   and X_CAL
%
% Written by Zhenqi (Pete) Shi
% Eli Lilly and Company
% 2020.4.2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clock
% h = waitbar(0,'Please wait...');

LV=size(PLS.loads{2,1},2);
t_cal=zeros(num_samp,LV);
p_cal=zeros(size(X_CAL,2),LV);
b_cal=zeros(1,LV);
t_val=zeros(size(X_VAL,1),LV);
options_1.display='off';
options_1.plots='none';
if isempty(X_aug)
    d_CAL=sqrt(sum(PLS.loads{1,1}(PLS.detail.includ{1,1},:).^2,2));
    X_CAL=X_CAL;
    Y_CAL=Y_CAL;
    t_CAL=PLS.loads{1,1}(PLS.detail.includ{1,1},:);
else
    pred=pls(X_aug,PLS,options_1);
    d_CAL=sqrt(sum([PLS.loads{1,1}(PLS.detail.includ{1,1},:);pred.loads{1,1}].^2,2));
    X_CAL=[X_CAL;X_aug];
    Y_CAL=[Y_CAL;Y_aug];
    t_CAL=[PLS.loads{1,1}(PLS.detail.includ{1,1},:);pred.loads{1,1}];
end

X_CAL = preprocess('calibrate',options,X_CAL);
X_VAL = preprocess('calibrate',options,X_VAL);
X_CAL=X_CAL.data;
X_VAL=X_VAL.data;
    
for i=1:size(X_VAL,1)
%     if FLAG_diff_pre==1
%         d=real(snv(X_CAL)-(ones(size(X_CAL,1),1)*snv(X_VAL(i,:))));   %log10 preprocessing can be customized
%     else

    X_CAL_1=[];
    Y_CAL_1=[];

    pred=pls(X_VAL(i,:),PLS,options_1);
    d_VAL=sqrt(sum(pred.loads{1,1}.^2,2));
    if max(d_VAL)>max(d_CAL)        %need to update this 2020/4/2
        OUTPUT.outlier=1;
    else
        OUTPUT.outlier=0;
    end
        
    if strcmp(FLAG,'vip')
        d=X_CAL-(ones(size(X_CAL,1),1)*X_VAL(i,:));
        d = d .* repmat(vip(PLS)',size(X_CAL,1),1) .* d;
        d = sqrt(sum(d,2));
    elseif strcmp(FLAG,'score')        
        d=(t_CAL-ones(size(t_CAL,1),1)*pred.loads{1,1}).^2;
        d = sqrt(sum(d,2));        
    end
    
    w=exp(-1*d/(std(d)*phi));   %the paper says standard deviation
    [B,I] = sort(w,'descend');
    w=w(I(1:num_samp));
    weight=diag(w);
    size(I);
    size(X_CAL);
    X_CAL_1=X_CAL(I(1:num_samp),:);
    Y_CAL_1=Y_CAL(I(1:num_samp),:);
    
    X_CAL_MEAN=sum(weight*X_CAL_1,1)/sum(w);
    Y_CAL_MEAN=sum(weight*Y_CAL_1,1)/sum(w);
    
    X_CAL_SCALE=X_CAL_1-ones(size(X_CAL_1,1),1)*X_CAL_MEAN;
    Y_CAL_SCALE=Y_CAL_1-ones(size(Y_CAL_1,1),1)*Y_CAL_MEAN;
    
    X_VAL_SCALE=X_VAL(i,:)-X_CAL_MEAN;
  
        
    for r=1:LV  
        
        [U,S,V]=svd(X_CAL_SCALE'*weight*Y_CAL_SCALE*Y_CAL_SCALE'*weight*X_CAL_SCALE);

        t_cal(:,r)=X_CAL_SCALE*V(:,1);
%         size(t_cal)
        p_cal(:,r)=(X_CAL_SCALE'*weight*t_cal(:,r))*inv(t_cal(:,r)'*weight*t_cal(:,r));
%         size(p_cal)
        b_cal(i,r)=(Y_CAL_SCALE'*weight*t_cal(:,r))*inv(t_cal(:,r)'*weight*t_cal(:,r));
%         size(b_cal)
        t_val(i,r)=X_VAL_SCALE*V(:,1);
%         size(t_val)
        Y_VAL_PRED(i,r)=Y_CAL_MEAN+t_val(i,1:r)*b_cal(i,1:r)';
        
        
        X_CAL_SCALE=X_CAL_SCALE-t_cal(:,r)*p_cal(:,r)';
        Y_CAL_SCALE=Y_CAL_SCALE-t_cal(:,r)*b_cal(i,r)';
        X_VAL_SCALE=X_VAL_SCALE-t_val(i,r)*p_cal(:,r)';
    end
    %     waitbar(i/size(X_VAL,1),h)
    OUTPUT.RES(i)=sum(X_VAL_SCALE.^2);
    OUTPUT.samp_ind(i,:)=I(1:num_samp);
end
% close (h)

if isempty(Y_VAL)
    RMSE=[];
elseif size(Y_VAL,1)>1
    RMSE=sqrt(sum((Y_VAL*ones(1,LV)-Y_VAL_PRED).^2)/size(Y_VAL,1));
else
    RMSE=sqrt((Y_VAL*ones(1,LV)-Y_VAL_PRED).^2/size(Y_VAL,1));
end

OUTPUT.Y_VAL_PRED=Y_VAL_PRED;
OUTPUT.RMSE=RMSE;


% OUTPUT.b_cal=b_cal;
% OUTPUT.p_cal=p_cal;
% OUTPUT.w=w;
% OUTPUT.d=d;
% OUTPUT.Y_CAL_MEAN=Y_CAL_MEAN;
% OUTPUT.t_val=t_val;        
% clock  