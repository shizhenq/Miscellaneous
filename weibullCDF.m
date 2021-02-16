function [ANS,PDR_fit,RMSEP,R2]=weibullCDF(PDR,Time,lambda_initial,k_initial,A_inital)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%weibull function CDF fitting for INDIVIDUAL or REPLICATE dissolution profiles
%ANS: colume 1 is lambda - scale factor, colume 2 is k - shape factor, colume 3 is A - potency factor
%PDR_fit: fitted percentage drug released
%PDR: actual percentage drug released on individual or replicate samples
%(each disso profile is presented as a row vector)
%Time: time points (must be a row vector)
%R2: coefficient of determination
%RMSEP: fitted error
%
% From Yuxiang Zhao @ Duquesne University
% Modified by Zhenqi (Pete) Shi @ 2020.5.19
% Modified by Zhenqi (Pete) Shi @ 2020.11.18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Time=ones(size(PDR,1),1)*Time;
fun=@(x)(sum(sum((PDR-x(3)*(1-exp(-(Time./x(1)).^x(2)))).^2))); %cost function_least squares
x0=[lambda_initial,k_initial,A_inital];
x=[];

Aeq=[0,1,0];    %k is fixed at unity for Loxo305
Beq=1;
A=[0,0,1];  %remove inequality constraint on 11/18/2020 and bring it back on 12/6/2020
b=105;
lb=[0 1 0]; % realize the ub and lb was better way to do inequality constraint on 12/6/2020
ub=[60 1 105];
x=fmincon(fun,x0,[],[],Aeq,Beq,lb,ub); % fitting with constraint
ans=x;
pdr_fit=x(1,3)*(1-exp(-(Time./x(1,1)).^x(1,2)));


ANS=ans;
PDR_fit= pdr_fit;

RMSEP=sqrt(sum((PDR-PDR_fit).^2,2)/length(Time));
 
for m=1:size(PDR,1)
    corrcoef(PDR(m,:),PDR_fit(m,:));
    R2(m)=ans(1,2).^2;
end
  
  

    