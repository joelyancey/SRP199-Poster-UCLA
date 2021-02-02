% MATLAB code by Joel Yancey for Buonomano/Blair collaboration experiment at UCLA.
% Brief: This code uses classification algorithms to train and test whether information is available
% in recordings taken from a population of neurons in the rat auditory cortex.
%
% Citation
% Chih-Chung Chang and Chih-Jen Lin, LIBSVM : a library for support vector machines.
% ACM Transactions on Intelligent Systems and Technology, 2:27:1--27:27, 2011.
% Software available at http://www.csie.ntu.edu.tw/~cjlin/libsvm
%
% Original poster
% Decoding Stimulus Features From Cortical Population Responses (Yancey, J., Halladay, L., DeGuzman, R.,
% Blair, T., & Buonomano, D.) Poster presented at UCLA Neuroscience Undergraduate Poster Fair in
% California, Los Angeles.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%% CONTROL BOARD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

RUN_SVM = (1);
SAVE_ACTIVITY = (0);

pick_animal = (6);
TargetStims = [1 4 7];
%TargetStims = [1:9];

AuditoryCellCriteria  = 30; % based on peak-valley difference of convolved raster
conv_width=5; %convolve raster width


%%%%% SVM SETTINGS: SET # OF BLOCKS AND # of TEST TRIALS %%%%%

nPull = 10; % number of trials to be pulled out for 'svmpredict'
nBlocks = 500; % number of blocks for which SVM will execute

WINDOW_SELECTION = (2); % 1=1st window. 2=2nd window.  
%window = [15 80];
%window = [15 60];
window = [10 50];
%window = [15 29; 30 44; 45 59];
%window = [15 29; 30 44; 45 59;60 74;75 90];
%window = [15 60;80 100];
%window = [15 20; 21 25; 26 30; 31 40; 41 90];
%window = [15 19; 20 24; 25 29; 30 39; 40 90];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

antiwindow = [650 680]; % ANTI-WINDOW
TRAIN_WITH_ANTI = (0); 
TEST_WITH_ANTI = (0);


%%% ignore %%%
[numwindow ~]=size(window);
if (WINDOW_SELECTION == 1) window=window+100; else window=window+200; end

%%% SET DATA / TASK %%%
Rasters={}; minNumTrials=[];
% ORDER EXPERIMENT
isOrderTask=0;
if (pick_animal==5) load('AC5Rast.mat'); Rasters{end+1}=double(AC5Rast); minNumTrials = [minNumTrials min([AC5SessInfo_O.numTrial])]; 
    isOrderTask=1; ClusterName=AC5CellNames_O; end
if (pick_animal==6) load('AC6Rast.mat'); Rasters{end+1}=double(AC6Rast); minNumTrials = [minNumTrials min([AC6SessInfo_O.numTrial])]; 
    isOrderTask=1; ClusterName=AC6CellNames_O; end
if (pick_animal==7) load('AC7Rast.mat'); Rasters{end+1}=double(AC7Rast); minNumTrials = [minNumTrials min([AC7SessInfo_O.numTrial])]; 
    isOrderTask=1; ClusterName=AC7CellNames_O; TargetStims = [1 2 3 4 5 6 7 8 9 10 11 12]; CORRECT_STIM_NUMS=(0); end
if (pick_animal==8) load('AC8Rast.mat'); Rasters{end+1}=double(AC8Rast); minNumTrials = [minNumTrials min([AC8SessInfo_O.numTrial])]; 
    isOrderTask=1; ClusterName=AC8CellNames_O; TargetStims = [1 2 3 4 5 6 7 8 9 10 11 12]; CORRECT_STIM_NUMS=(0); end



% INTERVAL EXPERIMENT
isIntervalTask=0;
if (0) load('AC6Rast_interval.mat'); Rasters{end+1}=double(AC6Rast_interval); minNumTrials = [minNumTrials min([AC6SessInfo_I.numTrial])]; 
    ClusterName=AC6CellNames_I; isIntervalTask=1; window=[200,300,500,700]; TargetStims=[2 3 4 5]; end




nAnimals=numel(Rasters);
nTrials=min(minNumTrials);
numStim=length(TargetStims);

%%%% PICK AUD CELLS USING PEAK-VALLEY DIFFERENCE OF CONVOLVED RASTERS %%%%%
fprintf('Contructing Activity...  0%')
totalNumAud=0; ACTIVITY=[]; ANTI_ACTIVITY=[];
for j=1:nAnimals
    
    Raster=Rasters{j};
    
    StartTrial  = find(Raster(:,1)==99999);
    spikeindex  = setdiff(1:length(Raster),StartTrial);
    maxt        = ceil(max(Raster(spikeindex,3)))+1;
    pret        = 100;
    pret2       = 200;
    numCells    = max(Raster(spikeindex,2)); %%%%%

    
    % AC7 and AC8 have 12 stimuli. This will convert the stimulus numbers to the same as AC6
    if ((max(Raster(spikeindex,1))==12) && CORRECT_STIM_NUMS)
        for i=1:length(Raster)
            if (Raster(i,1)==1 || Raster(i,1)==5 || Raster(i,1)==9)  Raster(i,1)=0;
            elseif (Raster(i,1)==2 || Raster(i,1)==3 || Raster(i,1)==4)  Raster(i,1)=Raster(i,1)-1;
            elseif (Raster(i,1)==6 || Raster(i,1)==7 || Raster(i,1)==8)  Raster(i,1)=Raster(i,1)-2;
            elseif (Raster(i,1)==10 || Raster(i,1)==11 || Raster(i,1)==12)  Raster(i,1)=Raster(i,1)-3;
            end
        end
    end
    
    
    Conv_Rast = zeros(numStim,numCells,maxt+pret,'single');
    Delta=zeros(numStim,numCells);
    stimLists=cell(1,numStim);
    count=0;
    
    

    for i=TargetStims
        count=count+1;
        list = find(Raster(spikeindex,1)==i);
        conv_rast=sparse(Raster(spikeindex(list),2),floor(Raster(spikeindex(list),3)+100),1,numCells,maxt+pret);
        Conv_Rast(count,:,:)= CONV_RASTER(full(conv_rast),conv_width,'gauss');
        data = squeeze(Conv_Rast(count,:,:));
        rast_store(i,:,:)=full(conv_rast);
        
        if count==1
            [Y maxIndex]= max(data');
            [Y minIndex]= min(data');
        end
        
        if (max(data')>maxIndex)
        [Y maxIndex]= max(data'); %test
        end
        if (min(data')<minIndex)
        [Y minIndex]= min(data'); %test
        end
        
        Delta(count,:) = max(data')-min(data');
        
        stimLists{count}=[find(Raster(spikeindex,1)==i)]; % saves processing time when making Activity matrix]
    end
    

    % select auditory cells based on peak/valley analysis
    AudCells = find(max(Delta)>AuditoryCellCriteria);
    
    %dummy = find(
    
%     for c=[14,26,27,42,43,49,114,162,177,182,197,211,212,234,237,238,253,255,256,263,264,280,284,287,292,304,309,320,326,327,328,343]
%         AudCells=AudCells(AudCells~=c);
%     end
    
    
%         for c=[131,144,207,305,306,314,325,338]
%         AudCells=AudCells(AudCells~=c);
%     end
    
    
    numAud=length(AudCells);
    totalNumAud=totalNumAud+numAud; %%%%%% probably not correct for multiple cells

    
    %%%%%%%%%%%% MAKE ACTIVITY MATRIX: [T1(S1-SX) T2(S1-SX) ...] %%%%%%%%%%%%%%
    Activity=zeros(numStim*nTrials,numAud*numwindow);
    AntiActivity=zeros(numStim*nTrials,numAud);
    for trial=1:nTrials
        trialList=find((Raster(spikeindex,4)==trial));
        %for stim=1:numStim
        
        if (isIntervalTask)
            for stim=2:5
                List=intersect(trialList,stimLists{stim-1});
                raster=sparse(Raster(spikeindex(List),2),floor(Raster(spikeindex(List),3)+pret),1,numCells,maxt+pret);
                %Activity((trial-1)*numStim+stim,numAud*(w-1)+1:numAud*w) =sum(raster(AudCells,window(stim-1)+15:window(stim-1)+80),2);
                Activity((trial-1)*numStim+stim,1:numAud) =sum(raster(AudCells,window(stim-1)+15:window(stim-1)+80),2);
            end
        end
        
        if (isOrderTask)
            for stim=1:numStim
                List=intersect(trialList,stimLists{stim});
                raster=sparse(Raster(spikeindex(List),2),floor(Raster(spikeindex(List),3)+pret),1,numCells,maxt+pret);
                for (w=1:numwindow)
                    Activity((trial-1)*numStim+stim,numAud*(w-1)+1:numAud*w)=sum(raster(AudCells,window(w,1):window(w,2)),2);
                    %AntiActivity((trial-1)*numStim+stim,numAud*(w-1)+1:numAud*w)=sum(raster(AudCells,antiwindow(1):antiwindow(2)),2);
                    %AntiActivity((trial-1)*numStim+stim,numAud*(w-1)+1:numAud*w)=sum(raster(AudCells,500:1000),2);
                end
            end
        end
        
        %AntiActivity=AntiActivity/(window(2)-window(1));
        
        
        %if TRAIN_WITH_ANTI==1 | TEST_WITH_ANTI==1
            AntiActivity((trial-1)*numStim+stim,numAud*(w-1)+1:numAud*w)=sum(raster(AudCells,window(w,1):window(w,2)),2);
        %end

        fprintf('\b\b\b%2.0f%%',floor((trial/nTrials)*100))
   
        
    end
    
    ACTIVITY=[ACTIVITY Activity];
    ANTI_ACTIVITY=[ANTI_ACTIVITY AntiActivity];
end

for (d=1:27) fprintf('\b'); end




if (SAVE_ACTIVITY)
%save Activity.mat Activity Raster ClusterName AudCells window rast_store spikeindex nTrials Conv_Rast
if (exist('AC5Rast')) save AC5Activity.mat Activity Raster ClusterName AudCells window spikeindex nTrials Conv_Rast rast_store conv_width
end
if (exist('AC6Rast')) save AC6Activity.mat Activity Raster ClusterName AudCells window spikeindex nTrials Conv_Rast rast_store conv_width
end
if (exist('AC7Rast')) save AC7Activity.mat Activity Raster ClusterName AudCells window spikeindex nTrials Conv_Rast rast_store conv_width
end
if (exist('AC8Rast')) save AC8Activity.mat Activity Raster ClusterName AudCells window spikeindex nTrials Conv_Rast rast_store conv_width
end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% SUPPORT VECTOR MACHINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ACTIVITYflipd=ACTIVITY';

%ACTIVITY=ACTIVITY-AntiActivity;

%%% ATTEMPT DATA SCALING
% Maxs=max(ACTIVITY);
% Mins=min(ACTIVITY);
% MaxMinusMins=Maxs-Mins;
% 
% MinsRep=repmat(Mins,numStim*nTrials,1);
% MaxMinusMinsRep=repmat(MaxMinusMins,numStim*nTrials,1);
% ACTIVITY=(ACTIVITY-MinsRep)./MaxMinusMinsRep;


if RUN_SVM~=1 return; end
%fprintf('Testing SVM (%d blocks x %d test trials = %d tests per stimulus)...\n',nBlocks,nPull,nBlocks*nPull)

nOut=numStim; % # of output cells (= # of conditions being tested)

% MAKE LABELS
dummy = ones(nOut)*-1;
for i=1:nOut  dummy(i,i)=1;  end
trainLabels=repmat(dummy,nTrials-nPull,1);
predictLabels=repmat(dummy,nPull,1);
bothLabels=[trainLabels;predictLabels];


%%%%%%%%%%%%%%%%%%% RUN MANY TESTS OF SVM PERFORMANCE %%%%%%%%%%%%%%%%%%%%%
fprintf('Running SVM...  0%')
aSTORE=zeros(nPull*nOut,nBlocks,nOut);
for i=1:nBlocks
    
    % FOR EACH BLOCK, PULL RANDOM TRIALS FOR CROSS-VALIDATION
    pull=(randsample(nTrials,nPull)-1)*nOut;
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          PULLindex=[];
    for j=1:nPull  PULLindex=horzcat(PULLindex,pull(j)+1:pull(j)+nOut); end
    TRAINindex=setdiff(1:nOut*nTrials,PULLindex');
    
    %dummyMax=repmat(max(ACTIVITY(TRAINindex',:)),[size(ACTIVITY,1),1]);
    
    % select data from window or antiwindow for training/testing
    if (TRAIN_WITH_ANTI == 0) 
        trainDATA=ACTIVITY(TRAINindex',:);
        %trainDATA=ACTIVITY(TRAINindex',:)./(dummyMax(TRAINindex,:)+.0001);
    else trainDATA=ANTI_ACTIVITY(TRAINindex',:);
    end
    if (TEST_WITH_ANTI == 0) 
        testDATA=ACTIVITY(PULLindex',:);
        %testDATA=ACTIVITY(PULLindex',:)./(dummyMax(PULLindex,:)+.0001);
    else testDATA=ANTI_ACTIVITY(PULLindex',:);
    end
    
    for j=1:nOut
        %svm = svmtrain(trainLabels(:,j),trainDATA,'-t 0 -b 1 -c 100 -e 0.01 -q');
        %svm = svmtrain(trainLabels(:,j),trainDATA,'-t 1 -b 1 -c 100 -e 0.01 -q');
        %svm = svmtrain(trainLabels(:,j),trainDATA,'-t 2 -b 0 -c 1 -q');
        svm = svmtrain(trainLabels(:,j),trainDATA,'-t 0 -b 0 -c .1 -q'); %<----
        [l, a, out] = svmpredict(predictLabels(:,j), testDATA, svm, '-q');
        
        aSTORE(:,i,j)=out; % STORAGE MATRIX FOR OUTPUT VALUES
        
    end
    
    fprintf('\b\b\b%2.0f%%',floor((i/nBlocks)*100))
    
end

for (d=1:18) fprintf('\b'); end

%%%%%%%%%%%%%%%%%% CALCULATE ACCURACY/ SVM PERFORMANCE %%%%%%%%%%%%%%%%%%%%

% aSTORE has size nOut*nPull, nBlocks, nOut
% essentially... all responses to all trial data of all stims x 1 block x 1 svm  

[row col depth]=size(aSTORE);
for i=1:nBlocks
   for j=1:row
      o = squeeze(aSTORE(j,i,:));
      [dummy choice(i,j)] = max(o);
      correct(i,j) = mod(j-1,nOut)+1;
   end
end

% 'choice' and 'correct' size nBlocks, nOut*nPull. 

dummy1=choice-correct;
dummy1(find(dummy1~=0))=9999;
dummy1(find(dummy1==0))=1;
dummy1(find(dummy1~=1))=0;
% 'dummy1' has same size as 'choice'/'correct'. 1 if correct, else 0.

correctperblock = sum(dummy1,2)/size(dummy1,2); % percent correct per block
STD_ERROR = std(correctperblock)/sqrt(nBlocks);

dummy2 = find((choice(:)-correct(:))==0); %makes dummy2 one dimensional
for i=1:nOut
    accuracies(i) = length(intersect(dummy2,find(correct(:)==i)))/(nPull*nBlocks);
    std_errors(i) = std(mean(choice(:,i:nOut:end)'))/sqrt(nBlocks);
end

ACCURACY = length(dummy2==0)/length(correct(:));
ACCURACY_dean = mean(correctperblock); % Dean's method



%%%%%%%%%%%%%%%%%%%%%% PRINT RESULTS TO CONSOLE %%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n----------------SUMMARY----------------\n')
fprintf('SETTINGS:\n')
if (isIntervalTask) fprintf('task=INTERVAL') 
else fprintf('task=ORDER | Target Stims: '); end
for i=1:numStim fprintf('%d ', TargetStims(i)); end
fprintf(' | Animal(s) used: ')
if exist('AC5Rast') fprintf('AC5 '); end; if exist('AC6Rast') fprintf('AC6 '); end
if exist('AC7Rast') fprintf('AC7 '); end; if exist('AC8Rast') fprintf('AC8 '); end
fprintf('| # cells: %d (criteria: %d)',totalNumAud,AuditoryCellCriteria) 
%fprintf('\nAud Cells: %d (criteria: %d)  |  Trials: %d  |  Response Window: %d:%d\n',totalNumAud,AuditoryCellCriteria,nTrials,min(window),max(window));
fprintf('\nTrials: %d',nTrials);
fprintf(' | # Blocks: %d | # Cross-validations: %d\n',nBlocks,nPull)
fprintf('RESULTS:\n')
fprintf('Overall accuracy: %2.4f (std error: %2.4f)\n',ACCURACY,STD_ERROR)
if ((ACCURACY-ACCURACY_dean>.000001) | (ACCURACY_dean-ACCURACY>.000001)) fprintf('WARNING: WARNING: Joels accuracy ~= Deans accuracy!'); end
%fprintf('Overall standard error: %2.4f\n',STD_ERROR)
for i=1:nOut fprintf('Output accuracy for preferred stimulus %d: %2.4f (std error: %2.4f)\n',i,accuracies(i),std_errors(i)); end

