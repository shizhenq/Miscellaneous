function OUTPUT=LW_PLS_v0(X_CAL,Y_CAL,X_VAL,Y_VAL,phi,PLS,num_samp,options)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LW-PLS per International Journal of Pharmaceutics 421 (2011) 269– 274
% with following modifications
%       1. using either vip-weighting to replace the
%            original distance calculation
%       2. adding prediction at each LV up to the LV defined in PLS via
%           conducting svd across latent variable one at a time.
%
% Inputs:
%   X_CAL,Y_CAL,X_VAL and Y_VAL: spectra and reference values of
%   calibration and test data. It is recommended to run one query sample at
%   a time.
%   phi: the parameter
%   num_samp: the number of local samples used for local PLS model
%   PLS: PLS model between X_CAL and Y_CAL.Be careful if there are samples
%   removed in PLS struct.
%   options: contain preprocessing struct for X_CAL and X_VAL. Mncn is not 
%   included. Please use the following line of command to generate such a 
%   struct. s=preprocess;
%
% Outputs:
%    Y_VAL_PRED: predicted values with individual value per PC used all 
%       the way up to the number of PC used in PLS.
%   RMSE: prediction error if Y_VAL is provided.
%   RES: spectral residual
%   outlier: if the query sample is a X outlier in the score space of PLS
%   samp_ind: the subset of training set was used to form LW-PLS.
%
%
% Written by Zhenqi (Pete) Shi
% Eli Lilly and Company
% 2020.4.14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LV=size(PLS.loads{2,1},2);
t_cal=zeros(num_samp,LV);
p_cal=zeros(size(X_CAL,2),LV);
b_cal=zeros(1,LV);
t_val=zeros(size(X_VAL,1),LV);
options.display='off';
options.plots='none';

for i=1:size(X_VAL,1)
    d_CAL=sqrt(sum(PLS.loads{1,1}(PLS.detail.includ{1,1},:).^2,2));
    t_CAL=PLS.loads{1,1}(PLS.detail.includ{1,1},:);
    
    pred=pls(X_VAL(i,:),PLS,options);
    d_samp=(t_CAL-ones(size(t_CAL,1),1)*pred.loads{1,1}).^2;
    d_samp = sqrt(sum(d_samp,2));
    [B_samp,I_samp] = sort(d_samp,'ascend');
    d_VAL=sqrt(sum(pred.loads{1,1}.^2,2));
    if (d_VAL)>max(d_CAL)   %updated @ 2020/4/14 and changed it back @ 2020/9/14
        OUTPUT.outlier=1;
    else
        OUTPUT.outlier=0;
    end
    
%     X_CAL = preprocess('calibrate',options,X_CAL);
%     X_VAL = preprocess('calibrate',options,X_VAL);
%     X_CAL=X_CAL.data;
%     X_VAL=X_VAL.data;

%     X_VAL=snv(savgol(X_VAL,15,2,1));
%     X_CAL=snv(savgol(X_CAL,15,2,1));

    d=X_CAL-(ones(size(X_CAL,1),1)*X_VAL(i,:)); %only VIP was deployed @ 2020/5/15
    d = d .* repmat(vip(PLS)',size(X_CAL,1),1) .* d;
    d = sqrt(sum(d,2));
    
    w=exp(-1*d/(std(d)*phi));   %the paper says standard deviation
    [B,I] = sort(w,'descend');
    w=w(I(1:num_samp));
    weight=diag(w);
    size(I);
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
        
        p_cal(:,r)=(X_CAL_SCALE'*weight*t_cal(:,r))*inv(t_cal(:,r)'*weight*t_cal(:,r));
        
        b_cal(i,r)=(Y_CAL_SCALE'*weight*t_cal(:,r))*inv(t_cal(:,r)'*weight*t_cal(:,r));
        
        t_val(i,r)=X_VAL_SCALE*V(:,1);
        
        Y_VAL_PRED(i,r)=Y_CAL_MEAN+t_val(i,1:r)*b_cal(i,1:r)';
        
        
        X_CAL_SCALE=X_CAL_SCALE-t_cal(:,r)*p_cal(:,r)';
        Y_CAL_SCALE=Y_CAL_SCALE-t_cal(:,r)*b_cal(i,r)';
        X_VAL_SCALE=X_VAL_SCALE-t_val(i,r)*p_cal(:,r)';
    end
    OUTPUT.RES(i)=sum(X_VAL_SCALE.^2);
    OUTPUT.samp_ind(i,:)=I(1:num_samp);
end

if isempty(Y_VAL)
    RMSE=[];
elseif size(Y_VAL,1)>1
    RMSE=sqrt(sum((Y_VAL*ones(1,LV)-Y_VAL_PRED).^2)/size(Y_VAL,1));
else
    RMSE=sqrt((Y_VAL*ones(1,LV)-Y_VAL_PRED).^2/size(Y_VAL,1));
end

OUTPUT.Y_VAL_PRED=Y_VAL_PRED;
OUTPUT.RMSE=RMSE;

