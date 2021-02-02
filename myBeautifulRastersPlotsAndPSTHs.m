clear

load('Activity.mat')

%%~~~~~~~~~~~~ CHOOSE STIMULI TO COMPARE ~~~~~~~~~~~~~%%%
COMPARE{1} = [1 4 7];   strings{1} = 'Low-second';
COMPARE{2} = [2 5 8];   strings{2} = 'High-second';
COMPARE{3} = [3 6 9];   strings{3} = 'Noise-second';
%stimLabels={'L,L','L,H','L,N','H,L','H,H','H,N','N,L','N,H','N,N'};

numStim=9;


%%%~~~~~~~~~~~~~~~ PSTH/PLOT PARAMETERS ~~~~~~~~~~~~~~%%%
dBin                     = 5; % bin size (ms)
range                    = 300; % range end (ms)
ctrlRange                = 400:1100;
nBins=floor(range/dBin);
centers=dBin/2:dBin:range-dBin;


%%%~~~~~~~~~~~~~ GET CELLS FROM USER INPUT ~~~~~~~~~~~~~%%%
userChoice1=input('View auditory cells? Enter [1]\nView choice cells? Enter [2]\n');
if userChoice1==1 cells=AudCells;
elseif userChoice1==2 cells=input('Enter cells (ex. [4 18 29]: ');
end

% userChoice1=input('View auditory cells? Enter [1]\nView choice cells? Enter [2]\n');
% if userChoice1==1 cells=AudCells;
% elseif userChoice1==2 cells=input('Enter cells (ex. [4 18 29]: ');
% end

spikes=zeros(numStim,length(cells),range); % - 'spikes' matrix stores spike times
spikes=rast_store(:,cells,1:range);


% baseLine response for each cell (across stims 1 through 9) in the range 'ctrlRange'
% responses/bin
%baseLine=sum(sum(squeeze(rast_store(:,cells,ctrlRange)),3),1)*dBin/(length(ctrlRange)*numStim);
% responses/sec
baseLine=sum(sum(squeeze(rast_store(:,cells,ctrlRange)),3),1)*1000/(length(ctrlRange)*numStim*nTrials);


%%%~~~~~~~~~~~~ INITIALIZE GRAPHIC OBJECTS ~~~~~~~~~~~~~%%%
hFig=figure;
pos=get(gcf,'Position').*[.5 .5 2.5 1.25]; % scale figure up in size
set(gcf,'position',pos,'doublebuffer','on','DefaultAxesColorOrder',[0 0 1;0 1 0;1 0 0]);

hRast=axes('position',[.04 .07 .63 .27]);
% - use 'Firing rate (Hz)'?

hPlots=[
    axes('position',[.04 .40 .30 .50]),...         % left plot
    axes('position',[.09+.80/3 .40 .30 .50]),...   % middle plot
    axes('position',[.14+1.6/3 .40 .30 .50])];   % right plot
hold on;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~ VIEW DATA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count=0;
for i=cells
    count=count+1;
    
    hold all
    
    unitID=sprintf('Unit #:  %d (aud cell#%d)      |     Session Name:  %s      |     Tetrode:  %s',i,find(AudCells==i),ClusterName{i}.session,ClusterName{i}.tetrode);
    fprintf(unitID);fprintf('\n');
    
    %%%~~~~~~~~~~~~~~~ STRINGS ~~~~~~~~~~~~~~~%%%
    hTitle=axes('Position',[0 0 1 1],'Visible','off');
    Title=text(.5,1,unitID,'Interpreter','none','VerticalAlignment','top',...
        'HorizontalAlignment','center','fontsize',18,'fontweight','bold','Units','normalized');
    
    string1= sprintf(['EXPERIMENTAL CONDITIONS:\n+ 0ms        - Trial start\n+ 100ms'... 
        '   - First tone onset\n+ 200ms   - Second tone onset\n+ ~2s          - Trial complete']);
    text1=text(.75,.26,string1,'units','normalized','fontsize',12);
    
%     string2=sprintf(['Avg. Response/Bin:  %0.2f\nNumber of Trials:  %d\nBin size:  %dms'],avgResponse(count),nTrials,dBin);
    string2=sprintf(['Baseline Frequency (spikes/s):  %0.2f\nNumber of Trials:  %d\nBin size:  %dms'],baseLine(count),nTrials,dBin);
    
    text2=text(.75,.10,string2,'units','normalized','fontsize',14);
    
    %%%~~~~~~~~~~~~~~~ RASTERS ~~~~~~~~~~~~~~~%%%
    axes(hRast);
    imagesc(squeeze(Conv_Rast(:,i,:)))
    colormap(gray)
    set(gca,'xlim',[0 800])
    line([100 100],get(gca,'ylim'),'Color','r')
    line([200 200],get(gca,'ylim'),'Color','r')
    ylabel('Stimulus Number','fontsize',15)
    xlabel('Time (ms)','fontsize',12)
    
    %%%~~~~~~~~~~~~~~~ PLOTS ~~~~~~~~~~~~~~~%%%
    temp1=[];
    for bin=1:dBin:range-dBin % assign spikes to bins
        temp1=[temp1 sum(squeeze(spikes(:,find(cells==i),bin:bin+dBin-1))')'];
    end
    
    temp1=temp1*1000/(dBin*nTrials); %change scaling to frequency (spikes/sec)
    
    axes(hPlots(1)); hold off
    plot(centers,temp1(COMPARE{1},:),'linewidth',1.5);
    title(strings{1},'fontsize',15);
    ylabel('Responses / Bin','fontsize',15)
    ylabel('Frequency (spikes/s)','fontsize',15)
    legend('L,L','H,L','N,L');
    set(gca,'YGrid','on')
    dummyA=get(gca,'ylim');
    
    axes(hPlots(2)); hold off
    plot(centers,temp1(COMPARE{2},:),'linewidth',1.5);
    title(strings{2},'fontsize',15);
    legend('L,H','H,H','N,H');
%     set(gca,'ytick',[])
    set(gca,'YGrid','on')
    dummyB=[dummyA get(gca,'ylim')];
    
    axes(hPlots(3)); hold off
    plot(centers,temp1(COMPARE{3},:),'linewidth',1.5);
    title(strings{3},'fontsize',15);
    legend('L,N','H,N','N,N');
%     set(gca,'ytick',[])
    set(gca,'YGrid','on')
    
    maxYlim=max([dummyB get(gca,'ylim')]); % maximum ylim of 3 plots

    for (j=1:3)   axes(hPlots(j));   % stuff to do to all the plots
        ylim([0 maxYlim]); % scales y-axis
%         text(.1,.1,sprintf('plot ylim: %d',maxYlim),'fontunits','normalized','position',[.85,.95],'units','normalized');
%         text(.1,.1,sprintf('avg. response/bin:\n %0.2f',avgResponse(count)),'fontunits','normalized','position',[.97,.95],'units','normalized','horizontalalignment','right');
%         text(.1,.1,string2,'fontunits','normalized','position',[.03,.95],'units','normalized','horizontalalignment','left');
        hold on
        line([100 100],get(gca,'ylim'),'Color','r')
        line([200 200],get(gca,'ylim'),'Color','r')
        hold on
    end

    waitforbuttonpress;
    cla(hTitle); % only way I could figure out to refresh title
end

close all;



