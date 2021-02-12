function [LABEL,DATA,PP,TIMEAXIS]=DISSO_pivot
% no repeat data entries are allowed in the excel sheet
[FileName,PathName,FilterIndex] = uigetfile( '*.*','MultiSelect', 'on');

[ndata, text, alldata] = xlsread([PathName '\' FileName]);

ind_1=(strcmp(alldata(1,:),'DissoResults_SampleUID')); % these column heads may need to be updated when coming across different excel sheet download from CRAVE
ind_2=(strcmp(alldata(1,:),'DissoResults_ResultName'));
ind_3=strcmp(alldata(1,:),'WeightCorrectedDissoResult');
ind_4=strcmp(alldata(1,:),'DissoResults_DissoTimepoint');
ind_5=strcmp(alldata(1,:),'NIR Sample ID');
ind_6=strcmp(alldata(1,:),'CompV4_TabSF');
ind_7=strcmp(alldata(1,:),'Tablet Thickness mm');
ind_8=strcmp(alldata(1,:),'CompV4_MainCompForce');
ind_9=strcmp(alldata(1,:),'CompV4_TensileStrength');
ind_10=strcmp(alldata(1,:),'Tablet Weight mg');
ind_PP=logical(sum([ind_6; ind_10; ind_7; ind_8; ind_9])); %combine tickness/tablet weight/main force/TSF/tensile strength into PP and in that order   

ind=zeros(size(alldata,1)-1,1);
ind_count=1;
for i=1:size(alldata,1)-2
    if strcmp([cell2str(alldata(i+1,ind_1)) '-' cell2str(alldata(i+1,ind_5)) '-' cell2str(alldata(i+1,ind_2))],[cell2str(alldata(i+2,ind_1)) '-' cell2str(alldata(i+2,ind_5)) '-' cell2str(alldata(i+2,ind_2))])
        ind(i:i+1)=ind_count;
    else
       ind=ind;
       ind_count=ind_count+1;
    end
end
alldata=alldata(2:end,:);
TIMEAXIS=cell2mat(alldata(ind==1,ind_4));
LABEL={};
DATA=[];
PP=[];
size(alldata)
for r=1:max(ind)    
    temp=find(ind==r);    
    temp_1=cell2mat(alldata(temp,ind_3));
    alldata(temp(1),ind_PP)
    temp_2=cell2mat(alldata(temp(1),ind_PP));   
%     [cell2str(alldata(temp(1),ind_1)) '-' cell2str(alldata(temp(1),ind_5)) '-' cell2str(alldata(temp(1),ind_2))]    
    if length(temp)==8 | length(temp)==9    % [5;10;15;20;30;45;60;75] or [5;10;15;20;30;45;60;75;90] 
        disp('Scenario#1')
        LABEL{size(LABEL,1)+1,1}=[cell2str(alldata(temp(1),ind_1)) '-' cell2str(alldata(temp(1),ind_5)) '-' cell2str(alldata(temp(1),ind_2))];
        DATA(size(DATA,1)+1,1:length(temp))=cell2mat(alldata(temp,ind_3))';
        PP(size(PP,1)+1,:)=cell2mat(alldata(temp(1),ind_PP));
    elseif length(temp)==16 | length(temp)==18  %repeat samples
        if diff(temp_1(1:2))~=0
            disp('Scenario#2a')
            LABEL{size(LABEL,1)+1,1}=[cell2str(alldata(temp(1),ind_1)) '-' cell2str(alldata(temp(1),ind_5)) '-' cell2str(alldata(temp(1),ind_2))];
            DATA(size(DATA,1)+1,1:length(temp)/2)=cell2mat(alldata(temp(1:length(temp)/2),ind_3))';
            PP(size(PP,1)+1,:)=cell2mat(alldata(temp(1),ind_PP));
        else
            disp('Scenario#2b')
            LABEL{size(LABEL,1)+1,1}=[cell2str(alldata(temp(1),ind_1)) '-' cell2str(alldata(temp(1),ind_5)) '-' cell2str(alldata(temp(1),ind_2))];
            DATA(size(DATA,1)+1,1:length(temp)/2)=cell2mat(alldata(temp(1:2:length(temp)),ind_3))';
            PP(size(PP,1)+1,:)=cell2mat(alldata(temp(1),ind_PP));
        end
    elseif length(temp)==32     %repeat samples
        if diff(temp_1(1:2))~=0
            disp('Scenario#3a')
            LABEL{size(LABEL,1)+1,1}=[cell2str(alldata(temp(1),ind_1)) '-' cell2str(alldata(temp(1),ind_5)) '-' cell2str(alldata(temp(1),ind_2))];
            DATA(size(DATA,1)+1,1:length(temp)/4)=cell2mat(alldata(temp(1:length(temp)/4),ind_3))';
            PP(size(PP,1)+1,:)=cell2mat(alldata(temp(1),ind_PP));
        else
            disp('Scenario#3b')
            LABEL{size(LABEL,1)+1,1}=[cell2str(alldata(temp(1),ind_1)) '-' cell2str(alldata(temp(1),ind_5)) '-' cell2str(alldata(temp(1),ind_2))];
            DATA(size(DATA,1)+1,1:length(temp)/4)=cell2mat(alldata(temp(1:4:length(temp)),ind_3))';
            PP(size(PP,1)+1,:)=cell2mat(alldata(temp(1),ind_PP));
        end
    else
        [cell2str(alldata(temp(1),ind_1)) '-' cell2str(alldata(temp(1),ind_5)) '-' cell2str(alldata(temp(1),ind_2))]
        LABEL=LABEL;
        DATA=DATA;
        PP=PP;
    end
end

% 
% find(DATA(:,1)~=nan);
% DATA=DATA(a,:);

    