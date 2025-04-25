function [D1] = loadAllDataIntoTable

% This file assumes the existence of csv files that summarize block results. 
% Grab current folder
folDir=pwd;

%% Dataset 1

% WP Group
baseDir='\data';
nsub=[2 4:7 9:11 13:17];   

for j=1:length(nsub) 
    if nsub(j)<10
        specDir='\WP0';
    else
        specDir='\WP';
    end
    exit=false;
    blockNo=0;
    subDir=[folDir baseDir specDir num2str(nsub(j)) '\'];

    while ~exit
        blockNo=blockNo+1;
        try
            t=readtable([subDir 'Slow' num2str(blockNo) '.csv']);            
        catch
            exit=true;
            break
        end

        t=t(~isnan(t.pertSize),:);
        
        taux=table(j.*ones(size(t.pertSize,1),1), blockNo*ones(size(t.pertSize,1),1), [1; zeros(size(t.pertSize,1)-1,1)], [zeros(size(t.pertSize,1)-1,1); 1] ,'VariableNames',{'subID','blockNo','isFirstInBlock','isLastInBlock'});
        t=cat(2,t,taux);

        if blockNo==1
            superT=t;
        else
            superT=cat(1,superT,t);
        end

    end

    if nsub(j)==nsub(1)
        WP=superT;
    else
        WP=cat(1,WP,superT);
    end

end


% WPS Group
superT = [];
superSuperT = [];
id_idx=unique(WP.subID); id_idx = id_idx(end);
nsub=[1 3 5 6 9 11:13 15:19];   

for j=1:length(nsub) 
    if nsub(j)<10
        specDir='\WPS0';
    else
        specDir='\WPS';
    end
    exit=false;
    blockNo=0;
    subDir=[folDir baseDir specDir num2str(nsub(j)) '\'];

    while ~exit
        blockNo=blockNo+1;
        try
            t=readtable([subDir 'Slow' num2str(blockNo) '.csv']);            
        catch
            exit=true;
            break
        end

        t=t(~isnan(t.pertSize),:);
        
        taux=table((id_idx+j).*ones(size(t.pertSize,1),1), blockNo*ones(size(t.pertSize,1),1), [1; zeros(size(t.pertSize,1)-1,1)], [zeros(size(t.pertSize,1)-1,1); 1] ,'VariableNames',{'subID','blockNo','isFirstInBlock','isLastInBlock'});
        t=cat(2,t,taux);

        if blockNo==1
            superT=t;
        else
            superT=cat(1,superT,t);
        end

    end

    if nsub(j)==nsub(1)
        WPS=superT;
    else
        WPS=cat(1,WPS,superT);
    end

end


% WP Group 3
nsub=[18:19 21:22 24 26:29 31:33 36];  
id_idx=unique(WPS.subID); id_idx = id_idx(end);

for j=1:length(nsub) 
    if nsub(j)<10
        specDir='\WP0';
    else
        specDir='\WP';
    end
    exit=false;
    blockNo=0;
    subDir=[folDir baseDir specDir num2str(nsub(j)) '\'];

    while ~exit
        blockNo=blockNo+1;
        try
            t=readtable([subDir 'Slow' num2str(blockNo) '.csv']);            
        catch
            exit=true;
            break
        end

        t=t(~isnan(t.pertSize),:);
        
        taux=table((id_idx+j).*ones(size(t.pertSize,1),1), blockNo*ones(size(t.pertSize,1),1), [1; zeros(size(t.pertSize,1)-1,1)], [zeros(size(t.pertSize,1)-1,1); 1] ,'VariableNames',{'subID','blockNo','isFirstInBlock','isLastInBlock'});
        t=cat(2,t,taux);

        if blockNo==1
            superT=t;
        else
            superT=cat(1,superT,t);
        end

    end

    if nsub(j)==nsub(1)
        WPF=superT;
    else
        WPF=cat(1,WPF,superT);
    end

end


%% Concatenate

Dataset1=[WP;WPS; WPF]; %;WPF
    
clear WP;
clear WPS;
clear WPF;

%% Add missing columns of information to the concatenated datalog

% Define quantities of interest

%Adding prev perturbation to table:
Dataset1.prevSize=[0;Dataset1.pertSize(1:end-1)]; 
Dataset1.prevSize(Dataset1.isFirstInBlock~=0)=0; %Assigning NaN to previous perturbation for first trial in each block

%Creating binary response variable(s):
Dataset1.leftResponse=Dataset1.initialResponse==-1;
Dataset1.rightResponse=Dataset1.initialResponse==1;
Dataset1.noResponse=isnan(Dataset1.initialResponse);
Dataset1.nullTrials=Dataset1.pertSize==0;
Dataset1.correctResponses=Dataset1.initialResponse==-sign(Dataset1.pertSize) & ~Dataset1.nullTrials;
Dataset1.incorrectResponses=Dataset1.initialResponse==sign(Dataset1.pertSize) & ~Dataset1.nullTrials;

%%  Delete the first perturbation size in each block

% Get rid of the first stimulus size per block except for WP02. The first 3 subjects collected did not have the very first stimulus size of the block repeated at the end. 
idx=find(Dataset1.isLastInBlock==1);
idx=idx(Dataset1.pertSize(idx(1:end-1))==Dataset1.pertSize(idx(1:end-1)+1)); %This is the index of the last in block perturbation where the controller most likely failed so we want to remove it
Dataset1(idx,:) = []; %Delete the rows for the stimulus size where the controller failed

Dataset1.isFirstInBlock(Dataset1.subID==1 & Dataset1.isFirstInBlock==1)=0; % Get rid of the first stimulus size of each block
% Dataset1.isFirstInBlock(Dataset1.subID==33 &  Dataset1.blockNo==7 & Dataset1.isFirstInBlock==1)=0; % We do not want to remove this perturbation when the controller broke
% Dataset1.isFirstInBlock(Dataset1.subID==35 &  Dataset1.blockNo==7 & Dataset1.isFirstInBlock==1)=0; % We do not want to remove this perturbation when the controller broke

Dataset1=Dataset1(Dataset1.isFirstInBlock~=1,:);
Dataset1=removevars(Dataset1,{'isFirstInBlock', 'isLastInBlock'});
% G = groupcounts(Dataset1,{'pertSize','subID'});


% Change block numbers per participant 

% WP04, subID=2
Dataset1.blockNo(Dataset1.subID == 2 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 2 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 2 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP05, subID=3
Dataset1.blockNo(Dataset1.subID == 3 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2 | Dataset1.blockNo == 3 )) = 1;
Dataset1.blockNo(Dataset1.subID == 3 & (Dataset1.blockNo == 4 | Dataset1.blockNo == 5)) = 2;
Dataset1.blockNo(Dataset1.subID == 3 & (Dataset1.blockNo == 6 | Dataset1.blockNo == 7)) = 3;

% WP06, subId=4
Dataset1.blockNo(Dataset1.subID == 4 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 4 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 4 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP07, subId=5
Dataset1.blockNo(Dataset1.subID == 5 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 5 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 5 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP09, subId=6
Dataset1.blockNo(Dataset1.subID == 6 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 6 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4 | Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 2;
Dataset1.blockNo(Dataset1.subID == 6 & (Dataset1.blockNo == 7 | Dataset1.blockNo == 8)) = 3;

% WP10, subId=7
Dataset1.blockNo(Dataset1.subID == 7 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 7 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 7 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP11, subId=8
Dataset1.blockNo(Dataset1.subID == 8 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 8 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 8 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP13, subId=9
Dataset1.blockNo(Dataset1.subID == 9 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 9 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 9 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP14, subId=10
Dataset1.blockNo(Dataset1.subID == 10 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 10 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 10 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP15, subId=11
Dataset1.blockNo(Dataset1.subID == 11 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 11 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 11 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP16, subId=12
Dataset1.blockNo(Dataset1.subID == 12 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 12 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 12 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP17, subId=13
Dataset1.blockNo(Dataset1.subID == 13 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 13 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 13 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WPS01, subId=14
Dataset1.blockNo(Dataset1.subID == 14 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 14 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 14 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WPS03, subId=15
Dataset1.blockNo(Dataset1.subID == 15 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 15 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 15 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WPS05, subId=16
Dataset1.blockNo(Dataset1.subID == 16 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2 | Dataset1.blockNo == 3)) = 1;
Dataset1.blockNo(Dataset1.subID == 16 & (Dataset1.blockNo == 4 | Dataset1.blockNo == 5)) = 2;
Dataset1.blockNo(Dataset1.subID == 16 & (Dataset1.blockNo == 6 | Dataset1.blockNo == 7)) = 3;

% WPS06, subId=17 
Dataset1.blockNo(Dataset1.subID == 17 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 17 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 17 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WPS09, subId=18 
Dataset1.blockNo(Dataset1.subID == 18 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 18 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 18 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WPS11, subId=19 
Dataset1.blockNo(Dataset1.subID == 19 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 19 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 19 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WPS12, subId=20 
Dataset1.blockNo(Dataset1.subID == 20 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 20 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 20 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WPS13, subId=21 
Dataset1.blockNo(Dataset1.subID == 21 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 21 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 21 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WPS15, subId=22 
Dataset1.blockNo(Dataset1.subID == 22 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 22 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 22 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WPS16, subId=23 
Dataset1.blockNo(Dataset1.subID == 23 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 23 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 23 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WPS17, subId=24 
Dataset1.blockNo(Dataset1.subID == 24 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 24 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 24 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WPS18, subId=25 
Dataset1.blockNo(Dataset1.subID == 25 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 25 & (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 25 & (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WPS19, subId=26 
Dataset1.blockNo(Dataset1.subID == 26 & (Dataset1.blockNo == 1 | Dataset1.blockNo == 2 | Dataset1.blockNo == 3)) = 1;
Dataset1.blockNo(Dataset1.subID == 26 & (Dataset1.blockNo == 4 | Dataset1.blockNo == 5)) = 2;
Dataset1.blockNo(Dataset1.subID == 26 & (Dataset1.blockNo == 6 | Dataset1.blockNo == 7)) = 3;

% WP18, subId=27 
Dataset1.blockNo(Dataset1.subID == 27 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 27 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2; 
Dataset1.blockNo(Dataset1.subID == 27 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP19, subId=28
Dataset1.blockNo(Dataset1.subID == 28 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 28 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2; 
Dataset1.blockNo(Dataset1.subID == 28 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP21, subId=29
Dataset1.blockNo(Dataset1.subID == 29 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 29 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4 )) = 2; 
Dataset1.blockNo(Dataset1.subID == 29 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6 | Dataset1.blockNo == 7)) = 3;

% WP22, subId=30 
Dataset1.blockNo(Dataset1.subID == 30 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 30 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4 )) = 2; 
Dataset1.blockNo(Dataset1.subID == 30 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6 )) = 3;

% WP24, subId=31 
Dataset1.blockNo(Dataset1.subID == 31 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 31 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4 )) = 2; 
Dataset1.blockNo(Dataset1.subID == 31 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6 )) = 3;

% WP26 subID = 32, this participants responses will be flipped 
Dataset1.blockNo(Dataset1.subID == 32 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 32 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 32 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

condition = Dataset1.subID == 32 & Dataset1.nullTrials ~= 1 & Dataset1.noResponse ~= 1;
Dataset1.leftResponse(condition) = ~Dataset1.leftResponse(condition);
Dataset1.rightResponse(condition) = ~Dataset1.rightResponse(condition);
Dataset1.initialResponse(Dataset1.subID == 32) = -1.*Dataset1.initialResponse(Dataset1.subID == 32);
Dataset1.correctResponses(condition & sign(Dataset1.initialResponse) ~= sign(Dataset1.pertSize)) = 1;
Dataset1.correctResponses(condition & sign(Dataset1.initialResponse) == sign(Dataset1.pertSize)) = 0;
Dataset1.incorrectResponses(condition & sign(Dataset1.initialResponse) ~= sign(Dataset1.pertSize)) = 0;
Dataset1.incorrectResponses(condition & sign(Dataset1.initialResponse) == sign(Dataset1.pertSize)) = 1; 

% WP27, subId=33
Dataset1.blockNo(Dataset1.subID == 33 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 33 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4 )) = 2; 
Dataset1.blockNo(Dataset1.subID == 33 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6 | Dataset1.blockNo == 7)) = 3;

% WP28, subId=34
Dataset1.blockNo(Dataset1.subID == 34 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 34 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2; 
Dataset1.blockNo(Dataset1.subID == 34 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP29, subId=35
Dataset1.blockNo(Dataset1.subID == 35 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 35 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2; 
Dataset1.blockNo(Dataset1.subID == 35 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6| Dataset1.blockNo == 7)) = 3;

% WP31, subId=36
Dataset1.blockNo(Dataset1.subID == 36 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 36 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2; 
Dataset1.blockNo(Dataset1.subID == 36 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP32, subId=37
Dataset1.blockNo(Dataset1.subID == 37 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 37 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2; 
Dataset1.blockNo(Dataset1.subID == 37 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP33, subId=38
Dataset1.blockNo(Dataset1.subID == 38 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 38 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 38 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

% WP36, subId=39
Dataset1.blockNo(Dataset1.subID == 39 &  (Dataset1.blockNo == 1 | Dataset1.blockNo == 2)) = 1;
Dataset1.blockNo(Dataset1.subID == 39 &  (Dataset1.blockNo == 3 | Dataset1.blockNo == 4)) = 2;
Dataset1.blockNo(Dataset1.subID == 39 &  (Dataset1.blockNo == 5 | Dataset1.blockNo == 6)) = 3;

D1 = Dataset1;

end