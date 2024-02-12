% 
%Main script for Time-Series Analyses
%Used in: 
% Pisanski, A., Prostebby, M., Dickson, D., Pagliardini, S. (2024).
% "Mapping responses to focal injections of bicuculline in the 
% lateral parafacial region identifies core regions for maximal generation 
% of active expiration" 

%Requires Functions: 
% LCdataload()
% datacensor()
% cycledefinev2()
% shadedErrorBar()

% Loads and analyses data from a lab chart .mat export. Outputs figures and
% saves the analyzed variables as a structure "annettesummary" in file:
% "Annette2023.mat"

%Written by Mitchell Prostebby, 2023

[DATA,INFO] = LCdataload();

%Data segmentation
artflag = 'y';
Overallbase = figure; 
while strcmp(artflag,'y')
    t = linspace(1,length(DATA(:,1)),length(DATA(:,1)));
    d = DATA(:,1);
    subplot(3,1,1)
    pp1 = plot(t,d); pp1.XDataSource = 't'; pp1.YDataSource = 'd';
    hold on;
    if isfield(INFO,'COMMENTTIMES')~= 0
        for i = 1:length(INFO.COMMENTTIMES)
            plot([1 1].*INFO.COMMENTTIMES(i,1),[min(d) max(d)],'--g');
        end
    end
    subplot(3,1,2); plot(t,DATA(:,2)); subplot(3,1,3); plot(t,DATA(:,3));
    artflag = input('Remove Artifacts?(y/n):','s');
    if strcmp(artflag,'y')
        [DATA(:,1),DATA(:,2),DATA(:,3)] = datacensor('nan',DATA(:,1),DATA(:,2),DATA(:,3));
        refreshdata(pp1,'base');
    end
end
    disp('select base and response period (ignoreNAN)')
    [o_base,~] = ginput(4);
    
basecorr = sum(isnan(DATA(1:o_base(2),1))); 
inj1corr = sum(isnan(DATA(o_base(2):o_base(3),1)));
inj2corr = sum(isnan(DATA(o_base(3):o_base(4),1)));
o_base(2:4) = o_base(2:4)-basecorr; o_base(3:4) = o_base(3:4)-inj1corr; o_base(4) = o_base(4)-inj2corr; 
DATA(isnan(DATA(:,1)),:) = []; %all channels 
close(Overallbase)
% DATA(1:o_base(1),:) = [];
data_norm = zeros(size(DATA));
% data_norm = zeros(round(o_base(4),0)-round(o_base(1),0)+1,3);

for i = 1:size(DATA,2)
    base = figure; plot(DATA(:,i));
    input('Zoom to desired area')
    [x,~] = ginput(2);
    close(base)
    % 
    % data_norm(:,i) = DATA(round(o_base(1),0):round(o_base(4),0),i)-mean(DATA(x(1):x(2),i)); %subtract the '0'-point for each signal
    data_norm(:,i) = DATA(:,i)-mean(DATA(x(1):x(2),i)); %subtract the '0'-point for each signal

    % data_norm(:,i) = data_norm(:,i)./max(data_norm(:,i));
    data_norm(:,i) = data_norm(:,i)./std(data_norm(:,i),[],"all");
    
    %PCA Norm
    PCAdata_norm(:,i) = DATA(:,i)-mean(DATA(o_base(1):o_base(2),i));
    PCAdata_norm(:,i) = PCAdata_norm(:,i)./max(PCAdata_norm(:,i));
    % data_norm(:,i) = DATA(:,i)-mean(DATA(:,i));
    % data_norm(:,i) = data_norm(:,i)./max(data_norm(:,i));
    % data_norm(:,i) = DATA(:,i)-mean(DATA(:,i));
    % data_norm(:,i) = data_norm(:,i)./(std(data_norm(:,i)).^2);

    % DATA(:,i) = DATA(:,i) - mean(DATA(x(1):x(2),i));
end

% %For quick reference to baseline and response period
% datacell = {data_norm(o_base(1):o_base(2),:),data_norm(o_base(3):o_base(4),:)};

% figure
% ax7 = subplot(3,1,1); plot(DATA(:,3),'k')
% ax8 = subplot(3,1,2); plot(DATA(:,2),'m')
% ax9 = subplot(3,1,3); plot(DATA(:,1),'r')
% linkaxes([ax7 ax8 ax9],'x')

groupcolor = 'k';
% groupcolor = "#D95319";
figure
ax7 = subplot(3,1,1); plot(linspace(0,length(data_norm(:,1))-1,length(data_norm(:,1)))./1000,data_norm(:,2),'Color',groupcolor,'LineWidth',0.5); ylabel('Airflow (STD)'); set(gca,'Box','off');
ax8 = subplot(3,1,2); plot(linspace(0,length(data_norm(:,1))-1,length(data_norm(:,1)))./1000,data_norm(:,3),'Color',groupcolor,'LineWidth',0.5); ylabel('\int Diaphragm (STD)'); set(gca,'Box','off');
ax9 = subplot(3,1,3); plot(linspace(0,length(data_norm(:,1))-1,length(data_norm(:,1)))./1000,data_norm(:,1),'Color',groupcolor,'LineWidth',0.5); ylabel('\int Abdominal (STD)'); set(gca,'Box','off');
hold on;  plot([1 1]*o_base(2)./1000,[min(data_norm(:,1)) max(data_norm(:,1))],'--g','LineWidth',2);
plot([1 1]*o_base(3)./1000,[min(data_norm(:,1)) max(data_norm(:,1))],'--g','LineWidth',2);
linkaxes([ax7 ax8 ax9],'x'); xlabel('Time (s)');

%% Average Cycle analysis
cutoff = 2; %HZ
[b,a] = butter(1,20./(1000/2));  

% 
strt_offset = 2000; %samples
% f_dia = filtfilt(b,a,DATA(:,3));
% [strt_ind,stimind,pkind,pkamp,pkwidth] = cycledefinev1(f_dia(strt_offset:end),[],cutoff,1000);
% strt_ind = strt_ind+strt_offset;
% pkind = pkind+strt_offset;


%Triggering on Airflow
f_air = filtfilt(b,a,data_norm(:,2));
% f_air = data_norm(:,2);
[ind,pkind2,pkamp2] = cycledefinev2(f_air(strt_offset:end),1);
strt_ind = ind+strt_offset;
pkind2 = pkind2+strt_offset;
strt_ind(strt_ind<o_base(1)) = [];
pkind2(pkind2<o_base(1)) = []; %Get rid of cycles before the beginning of the baseline period

preinj_period = find(strt_ind<o_base(2),1,'last');
half_triglen = 2*1000;%s * fs
if strt_ind(end)+half_triglen>length(DATA(:,1))
    strt_ind(end) = [];
end
timelag = linspace(-half_triglen./1000,half_triglen./1000,half_triglen.*2+1)';

for i = 1:length(strt_ind)
    TrigDia(:,i) = data_norm((strt_ind(i)-half_triglen):(strt_ind(i)+half_triglen),3);
    TrigAir(:,i) = data_norm((strt_ind(i)-half_triglen):(strt_ind(i)+half_triglen),2);
    TrigAbd(:,i) = data_norm((strt_ind(i)-half_triglen):(strt_ind(i)+half_triglen),1);
end
% zTrigDia = (TrigDia-nanmean(TrigDia,'all'))./std(TrigDia,[],"all",'omitnan'); %0centering and putting into STD units
% zTrigAir = (TrigAir-nanmean(TrigAir,'all'))./std(TrigAir,[],"all",'omitnan'); %0centering and putting into STD units
% zTrigAbd= (TrigAbd-nanmean(TrigAbd,'all'))./std(TrigAbd,[],"all",'omitnan'); %0centering and putting into STD units

%Baseline Averages
mTrigDia_base = mean(TrigDia(:,1:preinj_period),2); 
mTrigAir_base = mean(TrigAir(:,1:preinj_period),2); 
mTrigAbd_base = mean(TrigAbd(:,1:preinj_period),2); 

errTrigDia_base = std(TrigDia(:,1:preinj_period),[],2);
errTrigAir_base = std(TrigAir(:,1:preinj_period),[],2);
errTrigAbd_base = std(TrigAbd(:,1:preinj_period),[],2);

figure
histogram(diff(strt_ind))
% f1 = find(diff(strt_ind)<1100);
% fm1 = find(diff(strt_ind)<1100)-1;
% fp1 = find(diff(strt_ind)<1100)+1;
% pkamp2(f1) = NaN; pkamp2(fm1) = NaN; pkamp2(fp1) = NaN; pkamp2(isnan(pkamp2)) = [];
% pkind2(f1) = NaN; pkind2(fm1) = NaN; pkind2(fp1) = NaN; pkind2(isnan(pkind2)) = [];
% strt_ind(f1) = NaN; strt_ind(fm1) = NaN; strt_ind(fp1) = NaN; strt_ind(isnan(strt_ind)) = [];


%% Post injection 
postinj_period = o_base(3);%in samples
postinj_nBins = 10; %each 2 minutes long?
bindur = 2; %mintes 
%bindur*postinj_nBins SHOULD EQUAL 20MINUTES
m_resprate(1,1) = mean(diff(pkind2(1:preinj_period))./1000); err_resprate(1,1) = std(diff(pkind2(1:preinj_period))./1000)./sqrt(length(pkind2(1:preinj_period)));
ring_ind(1,1) = half_triglen-round(0.4*m_resprate(1,1)*1000,0);
ring_ind(1,2) = ring_ind(1,1)+round((m_resprate(1,1)*1000),0);
ring_sample(1) = ring_ind(1,2)-ring_ind(1,1)+1;

%Distance plots
TrigDia_minusbase = TrigDia-mTrigDia_base;
TrigAir_minusbase = TrigAir-mTrigAir_base;
TrigAbd_minusbase = TrigAbd-mTrigAbd_base;

%Euclidean Distance
    EucDistcycle = sqrt((TrigDia_minusbase.^2)+(TrigAir_minusbase.^2)+(TrigAbd_minusbase.^2));
    mEucDist_base = mean(EucDistcycle(:,1:preinj_period),2);
    % errEucDist_base = std(EucDistcycle(:,1:preinj_period,[],2)./sqrt(preinj_period);
    errEucDist_base = std(EucDistcycle(:,1:preinj_period),[],2);
    
    %Mahalanobis distance
    % MahDist = nan(size(TrigDia')); %cycles along the rows, timelags across the columns
    % for k = 1:length(timelag)
    %     MahDist(:,k) = mahal(horzcat(TrigDia(k,:)',TrigAir(k,:)',TrigAbd(k,:)'),horzcat(TrigDia(k,1:preinj_period)',TrigAir(k,1:preinj_period)',...
    %         TrigAbd(k,1:preinj_period)'));
    % end
    mMahDist_base = nan(length(timelag),1);
    for k = 1:length(timelag)
    mMahDist_base(k,1) = mahal([mTrigDia_base(k,1),mTrigAir_base(k,1),mTrigAbd_base(k,1)],horzcat(TrigDia(k,1:preinj_period)',TrigAir(k,1:preinj_period)',...
            TrigAbd(k,1:preinj_period)'));
    end

%Down sampleing for scatter plots
d_sample = 30000;
d_sample = d_sample*4/(postinj_nBins+1); %per segment
% persection = round(sqrt(d_sample/(postinj_nBins+1)),0);
random_mat_base = ones(length(timelag),preinj_period);
random_mat_base(ring_ind(1,1):ring_ind(1,2),:) = rand(ring_sample(1),preinj_period); randselection_base = false(length(timelag),preinj_period);
randselection_base(random_mat_base<(d_sample/(preinj_period*(ring_sample(1))))) = true;
sampled_base = sum(randselection_base,"all");
 
% colors = horzcat(linspace(1,0,postinj_nBins+1)',ones(postinj_nBins+1,1)*0.25,linspace(0,1,postinj_nBins+1)');
colors = turbo(postinj_nBins+1);
colors = vertcat([0,0,0],colors,[1 1 1]);
colors(2,:) = [0,0,0];
lwidth = 2.5;
malph1 = 0.13; msize = 36;
malph2 = 0.05;
malph4 = 0.3;

roundout = @(x) ceil(abs(x)).*sign(x); %function to round away from 0;
f1 = figure;
ax3d = axes();
s3 = scatter3(reshape(TrigAbd(randselection_base),sampled_base,1),reshape(TrigAir(randselection_base),sampled_base,1),reshape(TrigDia(randselection_base),...
    sampled_base,1),[],colors(2,:),'filled'); hold on;

% plot3(linspace(roundout(min(TrigAbd,[],'all')),roundout(max(TrigAbd,[],'all')),2),[0 0],[0 0],'--k','LineWidth',2); hold on;
% plot3([0 0],linspace(roundout(min(TrigAir,[],'all')),roundout(max(TrigAir,[],'all')),2),[0 0],'--k','LineWidth',2); hold on;
% plot3([0,0],[0,0],linspace(roundout(min(TrigDia,[],'all')),roundout(max(TrigDia,[],'all')),2),'--k','LineWidth',2); hold on;

plot3(linspace(-5,12,2),[0 0],[0 0],'--k','LineWidth',2); hold on;
plot3([0 0],linspace(-3,4,2),[0 0],'--k','LineWidth',2); hold on;
plot3([0,0],[0,0],linspace(-1,5,2),'--k','LineWidth',2); hold on;

% f2 = figure;
p1 = plot3(mTrigAbd_base(ring_ind(1,1):ring_ind(1,2),1),mTrigAir_base(ring_ind(1,1):ring_ind(1,2),1),mTrigDia_base(ring_ind(1,1):ring_ind(1,2),1),'Color',colors(2,:),'LineWidth',lwidth); hold on;
% colororder(f1,colors)

pp1 = figure;
% ax1 = subplot(1,2,1);

% 'linestyle','--','edgecolor','g',...
%     'linewidth',3,'edgealpha',0.2
% partS1b = scatter(TrigAir(randselection_base),TrigDia(randselection_base),[],colors(2,:),'filled'); hold on;
% partP1b = plot(mTrigAir_base(ring_ind(1,1):ring_ind(1,2),1),mTrigDia_base(ring_ind(1,1):ring_ind(1,2),1),'Color',colors(2,:),'LineWidth',lwidth); hold on;
% partS1b = scatter(TrigAir(randselection_base),TrigDia(randselection_base),[],colors(2,:),'filled'); hold on;
% partP1b = patchline(mTrigAir_base(ring_ind(1,1):ring_ind(1,2),1),mTrigDia_base(ring_ind(1,1):ring_ind(1,2),1),...
%     'linestyle','-','edgecolor',colors(2,:),...
%     'linewidth',lwidth,'edgealpha',malph4); hold on;
% plot(linspace(-3,4,2),[0,0],'--k','LineWidth',2); hold on;
% plot([0,0],linspace(-1,5,2),'--k','LineWidth',2); hold on;

% plot(linspace(roundout(min(TrigAir,[],'all')),roundout(max(TrigAir,[],'all')),2),[0,0],'--k','LineWidth',2); hold on;
% plot([0,0],linspace(roundout(min(TrigDia,[],'all')),roundout(max(TrigDia,[],'all')),2),'--k','LineWidth',2); hold on;

% 'linestyle','-','edgecolor',colors(2,:),...
%     'linewidth',lwidth,'edgealpha',malph2)
% ax2 = subplot(1,2,2); 
pp2 = figure;
% partS2b = scatter(TrigAir(randselection_base),TrigAbd(randselection_base),[],colors(2,:),'filled'); hold on;
% partP2b = plot(mTrigAir_base(ring_ind(1,1):ring_ind(1,2),1),mTrigAbd_base(ring_ind(1,1):ring_ind(1,2),1),'Color',colors(2,:),'LineWidth',lwidth); hold on;
% partS2b = scatter(TrigAir(randselection_base),TrigAbd(randselection_base),[],colors(2,:),'filled'); hold on;
% partP2b = patchline(mTrigAir_base(ring_ind(1,1):ring_ind(1,2),1),mTrigAbd_base(ring_ind(1,1):ring_ind(1,2),1),...
%     'linestyle','-','edgecolor',colors(2,:),...
%     'linewidth',lwidth,'edgealpha',malph4); hold on;
% plot(linspace(-3,4,2),[0,0],'--k','LineWidth',2); hold on;
% plot([0,0],linspace(-5,12,2),'--k','LineWidth',2); hold on;

% plot(linspace(roundout(min(TrigAir,[],'all')),roundout(max(TrigAir,[],'all')),2),[0,0],'--k','LineWidth',2); hold on;
% plot([0,0],linspace(roundout(min(TrigAbd,[],'all')),roundout(max(TrigAbd,[],'all')),2),'--k','LineWidth',2); hold on;

%preallocation
mMahDist_resp = nan(length(timelag),postinj_nBins);
baseringdef = sqrt((diff(mTrigDia_base(:,1)).^2)+(diff(mTrigAir_base(:,1)).^2)+(diff(mTrigAbd_base(:,1)).^2));
for j = 1:postinj_nBins
    binstart = find(strt_ind>(((j-1)*bindur)*60*1000)+postinj_period,1,'first');   
    binend = find(strt_ind<(j*bindur*60*1000)+postinj_period,1,'last');
    TAbdresp = TrigAbd(:,binstart:binend); TAirresp = TrigAir(:,binstart:binend); TDiaresp = TrigDia(:,binstart:binend);
    hfig = figure;
    dpk = diff(pkind2(binstart:binend));
    histogram(dpk);
    [zzz,~] = ginput(2); dpk(dpk<zzz(1)) = []; dpk(dpk>zzz(2)) = [];
    m_resprate(j+1,1) = mean(dpk./1000);
    close(hfig)
    % err_resprate(j+1,1) = std(dpk)./1000./sqrt(length(pkind2(binstart:binend)));
    ring_ind(j+1,1) = half_triglen-round(0.4*m_resprate(j+1,1)*1000,0); 
    ring_ind(j+1,2) = ring_ind(j+1,1)+round((m_resprate(j+1,1)*1000),0);
    ring_sample(j+1) = ring_ind(j+1,2)-ring_ind(j+1,1)+1;
    mTrigDia_resp(:,j) = mean(TDiaresp,2);
    mTrigAir_resp(:,j) = mean(TAirresp,2);
    mTrigAbd_resp(:,j) = mean(TAbdresp,2);
    
    % errTrigDia_resp(:,j) = std(TDiaresp,[],2)./sqrt(binend-binstart+1);
    % errTrigAir_resp(:,j) = std(TAirresp,[],2)./sqrt(binend-binstart+1);
    % errTrigAbd_resp(:,j) = std(TAbdresp,[],2)./sqrt(binend-binstart+1);

    errTrigDia_resp(:,j) = std(TDiaresp,[],2);
    errTrigAir_resp(:,j) = std(TAirresp,[],2);
    errTrigAbd_resp(:,j) = std(TAbdresp,[],2);
    
    figure(f1);
    random_mat_resp = ones(length(timelag),binend-binstart);
    random_mat_resp(ring_ind(j+1,1):ring_ind(j+1,2),:) = rand(ring_sample(j+1),binend-binstart); 
    randselection_resp = false(length(timelag),binend-binstart);
    randselection_resp(random_mat_resp<(d_sample/((binend-binstart)*ring_sample(j+1)))) = true;
    sampled_resp = sum(randselection_resp,"all");
    scs.(strcat('sc',num2str(3+j))) = scatter3(reshape(TAbdresp(randselection_resp),sampled_resp,1),reshape(TAirresp(randselection_resp),sampled_resp,1),reshape(TDiaresp(randselection_resp),...
    sampled_resp,1),[],colors(j+2,:),'filled'); hold on;
    p2 = plot3(mTrigAbd_resp(ring_ind(j+1,1):ring_ind(j+1,2),j),mTrigAir_resp(ring_ind(j+1,1):ring_ind(j+1,2),j),mTrigDia_resp(ring_ind(j+1,1):ring_ind(j+1,2),j),'Color',colors(j+2,:),'LineWidth',2.5);
    
    figure(pp1)
    partS1r.(strcat('n',num2str(j))) = scatter(TAirresp(randselection_resp),TDiaresp(randselection_resp),[],colors(j+2,:),'filled'); hold on;
    % partP1r = plot(mTrigAir_resp(ring_ind(j+1,1):ring_ind(j+1,2),j),mTrigDia_resp(ring_ind(j+1,1):ring_ind(j+1,2),j),'Color',colors(j+2,:),'LineWidth',lwidth); hold on;
    partP1r = patchline(mTrigAir_resp(ring_ind(j+1,1):ring_ind(j+1,2),j),mTrigDia_resp(ring_ind(j+1,1):ring_ind(j+1,2),j),...
        'linestyle','-','edgecolor',colors(j+2,:),...
    'linewidth',lwidth,'edgealpha',malph4); hold on;
    
    figure(pp2)
    partS2r.(strcat('n',num2str(j))) = scatter(TAirresp(randselection_resp),TAbdresp(randselection_resp),[],colors(j+2,:),'filled'); hold on;
    % partP2r = plot(mTrigAir_resp(ring_ind(j+1,1):ring_ind(j+1,2),j),mTrigAbd_resp(ring_ind(j+1,1):ring_ind(j+1,2),j),'Color',colors(j+2,:),'LineWidth',lwidth); hold on;
    partP2r = patchline(mTrigAir_resp(ring_ind(j+1,1):ring_ind(j+1,2),j),mTrigAbd_resp(ring_ind(j+1,1):ring_ind(j+1,2),j),...
        'linestyle','-','edgecolor',colors(j+2,:),...
    'linewidth',lwidth,'edgealpha',malph4); hold on;

    mTrigDiaresp_mbase(:,j) = mean(TrigDia_minusbase(:,binstart:binend),2);
    mTrigAirresp_mbase(:,j) = mean(TrigAir_minusbase(:,binstart:binend),2);
    mTrigAbdresp_mbase(:,j) = mean(TrigAbd_minusbase(:,binstart:binend),2);

    % errTrigDiaresp_mbase(:,j) = std(TrigDia_minusbase(:,binstart:binend),[],2)./sqrt(binend-binstart+1);
    % errTrigAirresp_mbase(:,j) = std(TrigAir_minusbase(:,binstart:binend),[],2)./sqrt(binend-binstart+1);
    % errTrigAbdresp_mbase(:,j) = std(TrigAbd_minusbase(:,binstart:binend),[],2)./sqrt(binend-binstart+1); 
    
    errTrigDiaresp_mbase(:,j) = std(TrigDia_minusbase(:,binstart:binend),[],2);
    errTrigAirresp_mbase(:,j) = std(TrigAir_minusbase(:,binstart:binend),[],2);
    errTrigAbdresp_mbase(:,j) = std(TrigAbd_minusbase(:,binstart:binend),[],2);

    mEucDist_resp(:,j) = mean(EucDistcycle(:,binstart:binend),2);
    % errEucDist_resp(:,j) = std(EucDistcycle(:,binstart:binend,[],2)./sqrt(binend-binstart+1);
    errEucDist_resp(:,j) = std(EucDistcycle(:,binstart:binend),[],2);

    % mMahDist_resp(:,j) = mean(MahDist(:,binstart:binend),2);
    for k = 1:length(timelag) %This version calcs Mahal for the mean resp (the ring) compared to the baseline cloud at each timelag (gives identical results to above)
        mMahDist_resp(k,j) = mahal([mTrigDia_resp(k,j),mTrigAir_resp(k,j),mTrigAbd_resp(k,j)],...
            horzcat(TrigDia(k,1:preinj_period)',TrigAir(k,1:preinj_period)',...
        TrigAbd(k,1:preinj_period)'));
    end
   
%Ring length plots
    ringdef(:,j) = sqrt((diff(mTrigDia_resp(:,j)).^2)+(diff(mTrigAir_resp(:,j)).^2)+(diff(mTrigAbd_resp(:,j)).^2));
    rrcumlen(:,j) = cumsum(ringdef(:,j))-cumsum(baseringdef); %relative lengths (how much deformity relative to baseline)
    rrlength(:,j) = ringdef(:,j)-baseringdef;
end
%Baseline
figure(pp1);
partS1b = scatter(TrigAir(randselection_base),TrigDia(randselection_base),[],colors(2,:),'filled'); hold on;
partP1b = patchline(mTrigAir_base(ring_ind(1,1):ring_ind(1,2),1),mTrigDia_base(ring_ind(1,1):ring_ind(1,2),1),...
    'linestyle','-','edgecolor',colors(2,:),...
    'linewidth',lwidth,'edgealpha',malph4); hold on;
plot(linspace(-3,4,2),[0,0],'--k','LineWidth',2); hold on;
plot([0,0],linspace(-1,5,2),'--k','LineWidth',2); hold on;

figure(pp2)
partS2b = scatter(TrigAir(randselection_base),TrigAbd(randselection_base),[],colors(2,:),'filled'); hold on;
partP2b = patchline(mTrigAir_base(ring_ind(1,1):ring_ind(1,2),1),mTrigAbd_base(ring_ind(1,1):ring_ind(1,2),1),...
    'linestyle','-','edgecolor',colors(2,:),...
    'linewidth',lwidth,'edgealpha',malph4); hold on;
plot(linspace(-3,4,2),[0,0],'--k','LineWidth',2); hold on;
plot([0,0],linspace(-5,12,2),'--k','LineWidth',2); hold on;

malph1 = 0.13; msize = 36;
% malph2 = 0.015;
s3.MarkerFaceAlpha = malph1;
partS1b.MarkerFaceAlpha = malph2;
partS2b.MarkerFaceAlpha = malph2;

s3.SizeData = msize;
% partS1b.SizeData = 1;
% partS2b.SizeData = 1;
% s3.MarkerEdgeAlpha = malph1;
% partS1b.MarkerEdgeAlpha = malph2;
% partS2b.MarkerEdgeAlpha = malph2;
for i = 1:postinj_nBins
    scs.(strcat('sc',num2str(3+i))).MarkerFaceAlpha = malph1;
    partS1r.(strcat('n',num2str(i))).MarkerFaceAlpha = malph2;
    partS2r.(strcat('n',num2str(i))).MarkerFaceAlpha = malph2;

    scs.(strcat('sc',num2str(3+i))).SizeData = msize;
    % partS1r.(strcat('n',num2str(i))).SizeData = 1;
    % partS2r.(strcat('n',num2str(i))).SizeData = 1;
    % scs.(strcat('sc',num2str(3+i))).MarkerEdgeAlpha = malph1;
    % partS1r.(strcat('n',num2str(i))).MarkerEdgeAlpha = malph2;
    % partS2r.(strcat('n',num2str(i))).MarkerEdgeAlpha = malph2;
end
figure(f1);
ax3d.Projection = 'orthographic';
cm = colormap(f1,colors(2:end-1,:));
c = colorbar;
cbarlabel = {'Baseline','2min Post','4min Post','6min Post','8min Post','10min Post','12min Post','14min Post', '16min Post','18min Post','20min Post'}';
% c.TickLabels = mat2cell(linspace(0,postinj_nBins,postinj_nBins+1),1,postinj_nBins+1);
c.TickLabels = cbarlabel;
xlabel('\int Abdominal (SD)'); ylabel('Airflow (SD)'); zlabel('\int Diaphragm (SD)');
xlim([-5 12]); ylim([-3 4]); zlim([-1 5])
set(gca,'ydir','reverse');
% view([-150,30]);
view([-40,20]);

% figure(pp1); set(gca,'xdir','reverse'); 
figure(pp1);
cm2 = colormap(pp1,colors(2:end-1,:));
c2 = colorbar; c2.TickLabels = cbarlabel;
xlabel('Airflow (STD)'); ylabel('\int Diaphragm (STD)')
xlim([-3 4]); ylim([-1 5]);
% figure(pp2); set(gca,'xdir','reverse'); 
figure(pp2);
cm2 = colormap(pp2,colors(2:end-1,:));
c2 = colorbar; c2.TickLabels = cbarlabel;
xlabel('Airflow (STD)'); ylabel('\int Abdominal (STD)')
xlim([-3 4]); ylim([-5 12]);
%% Distance figure explanation
% d2 = figure;
% axD2 = axes();
% malph3 = 0.065;
% base3 = scatter3(reshape(TrigAbd(randselection_base),sampled_base,1),reshape(TrigAir(randselection_base),sampled_base,1),reshape(TrigDia(randselection_base),...
%     sampled_base,1),[],'k','filled'); hold on;
% basep3 = patchline(mTrigAbd_base(ring_ind(1,1):ring_ind(1,2),1),mTrigAir_base(ring_ind(1,1):ring_ind(1,2),1),mTrigDia_base(ring_ind(1,1):ring_ind(1,2),1),...
%     'linestyle','-','edgecolor','k',...
%     'linewidth',lwidth,'edgealpha',0.4); hold on;
% % colororder(f1,colors)
% base3.MarkerFaceAlpha = malph3;
% 
% binstart = find(strt_ind>(((1-1)*bindur)*60*1000)+postinj_period,1,'first');   
% binend = find(strt_ind<(1*bindur*60*1000)+postinj_period,1,'last');
% 
% TAbdresp = TrigAbd(:,binstart:binend); TAirresp = TrigAir(:,binstart:binend); TDiaresp = TrigDia(:,binstart:binend);
% random_mat_resp = ones(length(timelag),binend-binstart);
%     random_mat_resp(ring_ind(2,1):ring_ind(2,2),:) = rand(ring_sample(2),binend-binstart); 
%     randselection_resp = false(length(timelag),binend-binstart);
%     randselection_resp(random_mat_resp<(d_sample/((binend-binstart)*ring_sample(2)))) = true;
% hold on;
% % s = rng;
% % rng(s);
% resp3 = scatter3(TAbdresp(randselection_resp),TAirresp(randselection_resp),TDiaresp(randselection_resp),...
%         [],[0.8500 0.3250 0.0980],'filled');
% respp3 = patchline(mTrigAbd_resp(ring_ind(2,1):ring_ind(2,2),1),mTrigAir_resp(ring_ind(2,1):ring_ind(2,2),1),mTrigDia_resp(ring_ind(2,1):ring_ind(2,2),1),...
%     'linestyle','-','edgecolor',[0.8500 0.3250 0.0980],...
%     'linewidth',lwidth,'edgealpha',0.4);
% resp3.MarkerFaceAlpha = malph3;
% hold on
% 
% pkpointB = find(max(mTrigAir_base)==mTrigAir_base);
% pkpointR = find(max(mTrigAir_resp(:,1))==mTrigAir_resp(:,1));
% qp = scatter3(mTrigAbd_resp(pkpointR,1),mTrigAir_resp(pkpointR,1),mTrigDia_resp(pkpointR,1),...
%     50,colors(3,:),'filled');
% s = rng;
% rng(s);
% ind = randi(preinj_period,1,1000);
% qp2 = scatter3(TrigAbd(pkpointB,ind),TrigAir(pkpointB,ind),TrigDia(pkpointB,ind),...
%     [],colors(1,:),'filled');
% qp2.MarkerFaceAlpha = 0.35;
% 
% plot3([mTrigAbd_resp(pkpointR,1);mTrigAbd_base(pkpointB,1)],...
%     [mTrigAir_resp(pkpointR,1);mTrigAir_base(pkpointB,1)],...
%     [mTrigDia_resp(pkpointR,1);mTrigDia_base(pkpointB,1)],...
%     'LineStyle','-','Color','m','LineWidth',3)
% 
% figure(d2);
% axD2.Projection = 'orthographic';
% % cm = colormap(d2,colors(2:end-1,:));
% % c = colorbar;
% % cbarlabel = {'Baseline','2min Post','4min Post','6min Post','8min Post','10min Post','12min Post','14min Post', '16min Post','18min Post','20min Post'}';
% % % c.TickLabels = mat2cell(linspace(0,postinj_nBins,postinj_nBins+1),1,postinj_nBins+1);
% % c.TickLabels = cbarlabel;
% xlabel('\int Abdominal (SD)'); ylabel('Airflow (SD)'); zlabel('\int Diaphragm (SD)');
% set(gca,'ydir','reverse');
% view([-40,20]);
%% Animation Creation
% %Requires above.
% 
% rotation_dt = 2;
% playback_dt = 0.075;
% nFrames = (360/rotation_dt); %spans 360 degrees of rotation
% figure(f1);
% malph1 = 1; msize = 5;
% s3.MarkerFaceAlpha = malph1;
% s3.SizeData = msize;
% for i = 1:postinj_nBins
%     scs.(strcat('sc',num2str(3+i))).MarkerFaceAlpha = malph1;
%     scs.(strcat('sc',num2str(3+i))).SizeData = msize;
% end
% ax3d.PositionConstraint = 'innerposition';
% ax3d.DataAspectRatioMode = 'manual';
% ax3d.PlotBoxAspectRatioMode = 'manual';
% ax3d.Projection = 'perspective';
% view([270,0])
% % [image, alpha] = export_fig(gcf,'transparent','nocrop','RGB');
% % [image,map] = rgb2ind(image,colors,'nodither');
% 
% rgbim = print('-RGBImage','-r300');
% frame = im2frame(rgbim);
% [image,map] = rgb2ind(frame.cdata,colors,'nodither');
% 
% % frame = getframe(f1);
% % image = nan(nFrames,nFrames,1,nFrames)
% % [image,map] = rgb2ind(frame.cdata,colors,'nodither');
% 
% for i = 2:(nFrames-1)
%     view([270-(i*rotation_dt),0])
%     % [image(:,:,1,i), alpha(:,:,1,i)] = export_fig(gcf);
%     % [image(:,:,1,i),map] = rgb2ind(image(:,:,1,i),colors,'nodither');
%     rgbim = print('-RGBImage','-r300');
%     frame = im2frame(rgbim);
%     image(:,:,1,i) = rgb2ind(frame.cdata,colors,'nodither');
% 
%     % frame = getframe(f1);
%     % image(:,:,1,i) = rgb2ind(frame.cdata,colors,'nodither');
% end
% % gifwrite(image, colors, {alpha},'C:\Users\mpros\Documents\MATLAB\CLAYTON\SilviaClay_Breathing\Animation2.gif', playback_dt,inf)
% 
% imwrite(image,map,'Supplemental_1.gif','DelayTime',playback_dt,'LoopCount',inf)
% 
% beep
% pause(0.5)
% beep
%% Calculations across time
%segments of the triggered cycle based on airflow and resprate
% pre_insp = nan(2,postinj_nBins+1);
% insp = nan(2,postinj_nBins+1);
% post_insp = nan(2,postinj_nBins+1);
a_pos = nan(2,postinj_nBins+1,3);

% a_pos(2,1,1) = find(mTrigAir_base(half_triglen:end)>0,1,'first') + half_triglen;
a_pos(2,1,1) = half_triglen+1;
% a_pos(1,1,1) = a_pos(2,1,1) - find(flip(mTrigAir_base(1:half_triglen))>0-thresh,1,'first');
cf = figure; plot(mTrigAir_base,'k'); disp('choose preinsp start'); [ps,~] = ginput(1);
if isempty(ps)
    % a_pos(1,1,1) = half_triglen;
    a_pos(1,1,1) = half_triglen - round((0.25*m_resprate(1)*1000),0); %say its supposed to be a qarter cycle backwards (based on p-114)
else
    a_pos(1,1,1) = round(ps,0);
end
close(cf)
thresh = 0.05;
a_pos(1,1,2) = a_pos(2,1,1); 
a_pos(2,1,2) = a_pos(1,1,2) + find(mTrigAir_base(a_pos(2,1,1):end)<0-thresh,1,'first');
a_pos(1,1,3) = a_pos(2,1,2); 
a_pos(2,1,3) = min([a_pos(1,1,3) + find(mTrigAir_base(a_pos(2,1,2):end)>0-thresh,1,'first'),a_pos(1,1,1)+round(m_resprate(1)*1000,0)]);

diadist = zeros(postinj_nBins,3);
airdist = zeros(postinj_nBins,3);
abddist = zeros(postinj_nBins,3);
eucdist = zeros(postinj_nBins,3);
mahaldist = zeros(postinj_nBins,3);
relringlen = zeros(postinj_nBins,3);
for i = 1:postinj_nBins
    % a_pos(2,i+1,1) = find(mTrigAir_resp(half_triglen:end,i)>0-thresh,1,'first')+half_triglen;
    a_pos(2,i+1,1) = half_triglen+1;

    cf = figure; plot(mTrigAir_resp(:,i)); disp('choose preinsp start'); [ps,~] = ginput(1);
    if isempty(ps)
        % a_pos(1,i+1,1) = half_triglen; %makes it so that the index difference is 1 which gets solved below
        a_pos(1,i+1,1) = half_triglen - round((0.25*m_resprate(i)*1000),0);
    else
        a_pos(1,i+1,1) = round(ps,0);
    end
    close(cf)
    % a_pos(1,i+1,1) = a_pos(2,i+1,1)-find(flip(mTrigAir_resp(1:half_triglen,i))>0-thresh,1,'first');
    a_pos(1,i+1,2) = a_pos(2,i+1,1); 
    a_pos(2,i+1,2) = a_pos(1,i+1,2) + find(mTrigAir_resp(a_pos(2,i+1,1):end,i)<0-thresh,1,'first');
    a_pos(1,i+1,3) = a_pos(2,i+1,2); 
    a_pos(2,i+1,3) = min([a_pos(1,i+1,3) + find(mTrigAir_resp(a_pos(2,i+1,2):end,i)>0-thresh,1,'first'),a_pos(1,i+1,1)+round(m_resprate(i+1)*1000,0)]);
    
    %calcs
    for pos = 1:3 %for each period
        % if a_pos(2,1,pos) - a_pos(1,1,pos) > 10
            n1 = a_pos(2,1,pos) - a_pos(1,1,pos)+1;
            % if a_pos(2,i+1,pos) - a_pos(1,i+1,pos) > 10
                n2 = a_pos(2,i+1,pos) - a_pos(1,i+1,pos)+1;
                    %Dia
                diatrigseg_base = mTrigDia_base(a_pos(1,1,pos):a_pos(2,1,pos),1);
                d1 = sum(diatrigseg_base);
                diatrigseg_r = mTrigDia_resp(a_pos(1,i+1,pos):a_pos(2,i+1,pos),i);
                d2 = sum(diatrigseg_r); %area of the response (all amplitudes are *1sample wide)
                n_diadist(i,pos) = (d2./(n2-n1))-d1;  
                rawdiadist(i,pos) = d2-d1;
                    %Air
                airtrigseg_base = mTrigAir_base(a_pos(1,1,pos):a_pos(2,1,pos),1);
                d1 = sum(airtrigseg_base);            
                airtrigseg_r = mTrigAir_resp(a_pos(1,i+1,pos):a_pos(2,i+1,pos),i);
                d2 = sum(airtrigseg_r);
                n_airdist(i,pos) = (d2./(n2-n1))-d1;
                rawairdist(i,pos) = d2-d1;
                    %Abd
                abdtrigseg_base = mTrigAbd_base(a_pos(1,1,pos):a_pos(2,1,pos),1);
                d1 = sum(abdtrigseg_base);
                abdtrigseg_r = mTrigAbd_resp(a_pos(1,i+1,pos):a_pos(2,i+1,pos),i);
                d2 = sum(abdtrigseg_r);
                n_abddist(i,pos) = (d2./(n2-n1))-d1;
                rawabddist(i,pos) = d2-d1;
                %Stretching
                  timevec = linspace(1,n2,n2); timevec = timevec*(n1/n2);
                  diastretch = nan(n1,1);
                  airstretch = nan(n1,1);
                  abdstretch = nan(n1,1);
                  for j = 1:n1
                    diastretch(j,1) = mean(diatrigseg_r(timevec>(j-1)&timevec<=(j)));  
                    airstretch(j,1) = mean(airtrigseg_r(timevec>(j-1)&timevec<=(j)));  
                    abdstretch(j,1) = mean(abdtrigseg_r(timevec>(j-1)&timevec<=(j)));  
                  end
                    diadiff = diastretch-diatrigseg_base;
                    airdiff = airstretch-airtrigseg_base;
                    abddiff = abdstretch-abdtrigseg_base;
                    diadist(i,pos) = sum(diadiff,'omitnan'); n2_diadist(i,pos) = diadist(i,pos)./(n1-sum(isnan(diadiff)));
                    airdist(i,pos) = sum(airdiff,'omitnan'); n2_airdist(i,pos) = airdist(i,pos)./(n1-sum(isnan(airdiff)));
                    abddist(i,pos) = sum(abddiff,'omitnan'); n2_abddist(i,pos) = abddist(i,pos)./(n1-sum(isnan(abddiff)));
                    eucdist(i,pos) = sum(sqrt((diadiff.^2)+(airdiff.^2)+(abddiff.^2)),'omitnan'); 
                    n2_eucdist(i,pos) = eucdist(i,pos)./(n1-sum(isnan(diadiff)));
                    mahaldiff = nan(n1,1);
                    for k = 1:n1
                        mahaldiff(k,1) = mahal([diastretch(k,1),airstretch(k,1),abdstretch(k,1)],...
                            horzcat(TrigDia(a_pos(1,1,pos)+k-1,1:preinj_period)',TrigAir(a_pos(1,1,pos)+k-1,1:preinj_period)',...
                            TrigAbd(a_pos(1,1,pos)+k-1,1:preinj_period)'));
                    end
                    mahaldist(i,pos) = sum(mahaldiff,'omitnan'); 
                    n2_mahaldist(i,pos) = mahaldist(i,pos)./(n1-sum(isnan(mahaldiff)));

            % rawdiadist(i,pos) = sum(mTrigDiaresp_mbase(a_pos(1,i+1,pos):a_pos(2,i+1,pos),i));
            % % n_diadist(i,pos) = diadist(i,pos)./(a_pos(2,i+1,pos)-a_pos(1,i+1,pos)); %not divided by 1000 becuase distances are areas with base of 1 sample
            % rawairdist(i,pos) = sum(mTrigAirresp_mbase(a_pos(1,i+1,pos):a_pos(2,i+1,pos),i));
            % % n_airdist(i,pos) = airdist(i,pos)./(a_pos(2,i+1,pos)-a_pos(1,i+1,pos));
            % rawabddist(i,pos) = sum(mTrigAbdresp_mbase(a_pos(1,i+1,pos):a_pos(2,i+1,pos),i));
            % % n_abddist(i,pos) = abddist(i,pos)./(a_pos(2,i+1,pos)-a_pos(1,i+1,pos));
            % raweucdist(i,pos) = sum(mEucDist_resp(a_pos(1,i+1,pos):a_pos(2,i+1,pos),i));
            % % n_eucdist(i,pos) = eucdist(i,pos)./(a_pos(2,i+1,pos)-a_pos(1,i+1,pos));
            % rawmahaldist(i,pos) = sum(mMahDist_resp(a_pos(1,i+1,pos):a_pos(2,i+1,pos),i));
            % % n_mahaldist(i,pos) = mahaldist(i,pos)./(a_pos(2,i+1,pos)-a_pos(1,i+1,pos));
            % rawrelringlen(i,pos) = sum(rrlength(a_pos(1,i+1,pos):a_pos(2,i+1,pos),i)); 
            % 
            % seg_ringlen_r(i,pos) = sum(ringdef(a_pos(1,i+1,pos):a_pos(2,i+1,pos),i))./(n2);
            % n_relringlen(i,pos) = seg_ringlen_r(i,pos)-(sum(baseringdef(a_pos(1,1,pos):a_pos(2,1,pos)))./n1);
            % 
    end %for each period
end %for each postinj bin

% figure
% subplot(3,1,1); yyaxis left; plot(m_resprate*1000); hold on; (plot(a_pos(2,:,1)-a_pos(1,:,1)+a_pos(2,:,2)-a_pos(1,:,2)+a_pos(2,:,3)-a_pos(1,:,3),'--b'));
% yyaxis right; plot(a_pos(2,:,1)-a_pos(1,:,1));
% subplot(3,1,2); yyaxis left; plot(m_resprate*1000);
% yyaxis right; plot(a_pos(2,:,2)-a_pos(1,:,2));
% subplot(3,1,3); yyaxis left; plot(m_resprate*1000);
% yyaxis right; plot(a_pos(2,:,3)-a_pos(1,:,3));
% 
% figure
% subplot(6,3,1); plot(diadist(:,1),'k'); title('STRETCH'); subplot(6,3,2); plot(diadist(:,2),'k'); subplot(6,3,3); plot(diadist(:,3),'k'); 
% subplot(6,3,4); plot(airdist(:,1),'m'); subplot(6,3,5); plot(airdist(:,2),'m'); subplot(6,3,6); plot(airdist(:,3),'m'); 
% subplot(6,3,7); plot(abddist(:,1),'r'); subplot(6,3,8); plot(abddist(:,2),'r'); subplot(6,3,9); plot(abddist(:,3),'r'); 
% subplot(6,3,10); plot(eucdist(:,1),'g'); subplot(6,3,11); plot(eucdist(:,2),'g'); subplot(6,3,12); plot(eucdist(:,3),'g'); 
% subplot(6,3,13); plot(mahaldist(:,1),'b'); subplot(6,3,14); plot(mahaldist(:,2),'b'); subplot(6,3,15); plot(mahaldist(:,3),'b'); 
% subplot(6,3,16); plot(rawrelringlen(:,1),'c'); subplot(6,3,17); plot(rawrelringlen(:,2),'c'); subplot(6,3,18); plot(rawrelringlen(:,3),'c'); 
% 
figure
subplot(3,3,1); plot(rawdiadist(:,1),'k'); title('RAW'); subplot(3,3,2); plot(rawdiadist(:,2),'k'); subplot(3,3,3); plot(rawdiadist(:,3),'k'); 
subplot(3,3,4); plot(rawairdist(:,1),'m'); subplot(3,3,5); plot(rawairdist(:,2),'m'); subplot(3,3,6); plot(rawairdist(:,3),'m'); 
subplot(3,3,7); plot(rawabddist(:,1),'r'); subplot(3,3,8); plot(rawabddist(:,2),'r'); subplot(3,3,9); plot(rawabddist(:,3),'r'); 
% subplot(6,3,10); plot(raweucdist(:,1),'g'); subplot(6,3,11); plot(raweucdist(:,2),'g'); subplot(6,3,12); plot(raweucdist(:,3),'g'); 
% subplot(6,3,13); plot(rawmahaldist(:,1),'b'); subplot(6,3,14); plot(rawmahaldist(:,2),'b'); subplot(6,3,15); plot(rawmahaldist(:,3),'b'); 
% subplot(6,3,16); plot(rawrelringlen(:,1),'c'); subplot(6,3,17); plot(rawrelringlen(:,2),'c'); subplot(6,3,18); plot(rawrelringlen(:,3),'c'); 

% figure
% subplot(6,3,1); plot(n_diadist(:,1),'k'); subplot(6,3,2); plot(n_diadist(:,2),'k'); subplot(6,3,3); plot(n_diadist(:,3),'k'); 
% subplot(6,3,4); plot(n_airdist(:,1),'m'); subplot(6,3,5); plot(n_airdist(:,2),'m'); subplot(6,3,6); plot(n_airdist(:,3),'m'); 
% subplot(6,3,7); plot(n_abddist(:,1),'r'); subplot(6,3,8); plot(n_abddist(:,2),'r'); subplot(6,3,9); plot(n_abddist(:,3),'r'); 
% subplot(6,3,10); plot(n2_eucdist(:,1),'g'); subplot(6,3,11); plot(n2_eucdist(:,2),'g'); subplot(6,3,12); plot(n2_eucdist(:,3),'g'); 
% subplot(6,3,13); plot(n2_mahaldist(:,1),'b'); subplot(6,3,14); plot(n2_mahaldist(:,2),'b'); subplot(6,3,15); plot(n2_mahaldist(:,3),'b'); 
% subplot(6,3,16); plot(n_relringlen(:,1),'c'); subplot(6,3,17); plot(n_relringlen(:,2),'c'); subplot(6,3,18); plot(n_relringlen(:,3),'c'); 

% figure
% subplot(6,3,1); plot(n2_diadist(:,1),'k'); title('NORM2'); subplot(6,3,2); plot(n2_diadist(:,2),'k'); subplot(6,3,3); plot(n2_diadist(:,3),'k'); 
% subplot(6,3,4); plot(n2_airdist(:,1),'m'); subplot(6,3,5); plot(n2_airdist(:,2),'m'); subplot(6,3,6); plot(n2_airdist(:,3),'m'); 
% subplot(6,3,7); plot(n2_abddist(:,1),'r'); subplot(6,3,8); plot(n2_abddist(:,2),'r'); subplot(6,3,9); plot(n2_abddist(:,3),'r'); 
% subplot(6,3,10); plot(n2_eucdist(:,1),'g'); subplot(6,3,11); plot(n2_eucdist(:,2),'g'); subplot(6,3,12); plot(n2_eucdist(:,3),'g'); 
% subplot(6,3,13); plot(n2_mahaldist(:,1),'b'); subplot(6,3,14); plot(n2_mahaldist(:,2),'b'); subplot(6,3,15); plot(n2_mahaldist(:,3),'b'); 
% subplot(6,3,16); plot(n_relringlen(:,1),'c'); subplot(6,3,17); plot(n_relringlen(:,2),'c'); subplot(6,3,18); plot(n_relringlen(:,3),'c'); 


% t = linspace(2,20,10); 
% figure
% subplot(3,3,1); plot(t,n2_airdist(:,1),'Color',groupcolor); title('Late Expiratory Period'); ylabel('Mean Airflow Difference'); set(gca,'Box','off'); ylim([-0.5 0.8]);
% subplot(3,3,2); plot(t,n2_airdist(:,2),'Color',groupcolor); title('Inspiratory Period'); ylabel('Mean Airflow Difference'); set(gca,'Box','off'); ylim([-0.5 0.8]);
% subplot(3,3,3); plot(t,n2_airdist(:,3),'Color',groupcolor); title('Early Expiratory Period'); ylabel('Mean Airflow Difference'); set(gca,'Box','off'); ylim([-0.5 0.8]);
% subplot(3,3,4); plot(t,n2_diadist(:,1),'Color',groupcolor); ylabel('Mean \int Diaphragm Difference'); set(gca,'Box','off'); ylim([-0.2 0.7]);
% subplot(3,3,5); plot(t,n2_diadist(:,2),'Color',groupcolor); ylabel('Mean \int Diaphragm Difference');set(gca,'Box','off'); ylim([-0.2 0.7]);
% subplot(3,3,6); plot(t,n2_diadist(:,3),'Color',groupcolor); ylabel('Mean \int Diaphragm Difference'); set(gca,'Box','off'); ylim([-0.2 0.7]);
% subplot(3,3,7); plot(t,n2_abddist(:,1),'Color',groupcolor); ylabel('Mean \int Abdominal Difference'); xlabel('Time Elapsed Since Injection (min)'); set(gca,'Box','off'); ylim([-0.1 5]);
% subplot(3,3,8); plot(t,n2_abddist(:,2),'Color',groupcolor); ylabel('Mean \int Abdominal Difference'); xlabel('Time Elapsed Since Injection (min)'); set(gca,'Box','off'); ylim([-0.1 5]);
% subplot(3,3,9); plot(t,n2_abddist(:,3),'Color',groupcolor);  ylabel('Mean \int Abdominal Difference'); xlabel('Time Elapsed Since Injection (min)'); set(gca,'Box','off'); ylim([-0.1 5]);
% 
% figure
% subplot(2,3,1); plot(t,n2_eucdist(:,1),'Color',groupcolor); title('Late Expiratory Period'); ylabel('Mean Euclidean Distance'); set(gca,'Box','off'); ylim([0 5]);
% subplot(2,3,2); plot(t,n2_eucdist(:,2),'Color',groupcolor); title('Inspiratory Period'); ylabel('Mean Euclidean Distance'); set(gca,'Box','off'); ylim([0 5]);
% subplot(2,3,3); plot(t,n2_eucdist(:,3),'Color',groupcolor); title('Early Expiratory Period'); ylabel('Mean Euclidean Distance');  set(gca,'Box','off'); ylim([0 5]);
% subplot(2,3,4); plot(t,n2_mahaldist(:,1),'Color',groupcolor); ylabel('Mean Mahalanobis Distance'); xlabel('Time Elapsed Since Injection (min)'); set(gca,'Box','off');  ylim([0 5000]);
% subplot(2,3,5); plot(t,n2_mahaldist(:,2),'Color',groupcolor); ylabel('Mean Mahalanobis Distance'); xlabel('Time Elapsed Since Injection (min)'); set(gca,'Box','off');  ylim([0 5000]);
% subplot(2,3,6); plot(t,n2_mahaldist(:,3),'Color',groupcolor); ylabel('Mean Mahalanobis Distance'); xlabel('Time Elapsed Since Injection (min)'); set(gca,'Box','off');  ylim([0 5000]);

%% Cycle averaged plots
groupcolor = 'b';
figure
% subplot(3,1,1); shadedErrorBar(timelag,mTrigAir_base,errTrigAir_base,'lineprops',groupcolor); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); xlim([-1 1]); ylabel('Airflow STD')
% subplot(3,1,2); shadedErrorBar(timelag,mTrigDia_base,errTrigDia_base,'lineprops',groupcolor); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]);  xlim([-1 1]); ylabel('\int Diaphragm (STD)')
% subplot(3,1,3); shadedErrorBar(timelag,mTrigAbd_base,errTrigAbd_base,'lineprops',groupcolor); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]);  xlim([-1 1]); ylabel('\int Abdominal (STD)')
subplot(3,1,1); shadedErrorBar(timelag,mTrigAir_base,errTrigAir_base,'lineprops',groupcolor); ylim([-3 4]); xlim([-1 1]); ylabel('Airflow STD')
subplot(3,1,2); shadedErrorBar(timelag,mTrigDia_base,errTrigDia_base,'lineprops',groupcolor); ylim([-1 5]);  xlim([-1 1]); ylabel('\int Diaphragm (STD)')
subplot(3,1,3); shadedErrorBar(timelag,mTrigAbd_base,errTrigAbd_base,'lineprops',groupcolor); ylim([-5 12]);  xlim([-1 1]); ylabel('\int Abdominal (STD)')

xlabel('Lag (s)'); 

%% period definition figure
% g1 = ((a_pos(1,2,1)-half_triglen+1)./1000);
% g2 = ((a_pos(2,2,1)-half_triglen+1)./1000);
% g3 = ((a_pos(2,2,2)-half_triglen+1)./1000);
% g4 = ((a_pos(2,2,3)-half_triglen+1)./1000);

g1 = a_pos(1,2,1);
g2 = a_pos(2,2,1);
g3 = a_pos(2,2,2);
g4 = a_pos(2,2,3);
figure
% plot(timelag,mTrigAir_resp(:,1),'Color',groupcolor,'LineWidth',2); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)');
subplot(3,1,1); plot(timelag(1:g1),...
    mTrigAir_resp(1:g1,1),'Color','k','LineWidth',2); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)');
hold on
plot(timelag(g1:g2),...
    mTrigAir_resp(g1:g2,1),'Color','b','LineWidth',2); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)');
plot(timelag(g2:g3),...
    mTrigAir_resp(g2:g3,1),'Color','r','LineWidth',2); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)');
plot(timelag(g3:g4),...
    mTrigAir_resp(g3:g4,1),'Color','c','LineWidth',2); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)');
plot(timelag(g4:end),...
    mTrigAir_resp(g4:end,1),'Color','k','LineWidth',2); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)');

hold on; plot([1 1].*((a_pos(1,2,1)-half_triglen+1)./1000),[min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,2,1)-half_triglen+1)./1000),[min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,2,2)-half_triglen+1)./1000),[min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,2,3)-half_triglen+1)./1000),[min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25],'--g','LineWidth',2);
xlim([-1 1]); set(gca,'Box','off');

subplot(3,1,2); plot(timelag(1:g1),...
    mTrigDia_resp(1:g1,1),'Color','k','LineWidth',2); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('\int Diaphragm (STD)');
hold on
plot(timelag(g1:g2),...
    mTrigDia_resp(g1:g2,1),'Color','b','LineWidth',2); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('\int Diaphragm (STD)');
plot(timelag(g2:g3),...
    mTrigDia_resp(g2:g3,1),'Color','r','LineWidth',2); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('\int Diaphragm (STD)');
plot(timelag(g3:g4),...
    mTrigDia_resp(g3:g4,1),'Color','c','LineWidth',2); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('\int Diaphragm (STD)');
plot(timelag(g4:end),...
    mTrigDia_resp(g4:end,1),'Color','k','LineWidth',2); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('\int Diaphragm (STD)');

hold on; plot([1 1].*((a_pos(1,2,1)-half_triglen+1)./1000),[min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,2,1)-half_triglen+1)./1000),[min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,2,2)-half_triglen+1)./1000),[min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,2,3)-half_triglen+1)./1000),[min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25],'--g','LineWidth',2);
xlim([-1 1]); set(gca,'Box','off');

subplot(3,1,3); plot(timelag(1:g1),...
    mTrigAbd_resp(1:g1,1),'Color','k','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 5]); ylabel('\int Abdominal (STD)');
hold on
plot(timelag(g1:g2),...
    mTrigAbd_resp(g1:g2,1),'Color','b','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 5]); ylabel('\int Abdominal (STD)');
plot(timelag(g2:g3),...
    mTrigAbd_resp(g2:g3,1),'Color','r','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 5]); ylabel('\int Abdominal (STD)');
plot(timelag(g3:g4),...
    mTrigAbd_resp(g3:g4,1),'Color','c','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 5]); ylabel('\int Abdominal (STD)');
plot(timelag(g4:end),...
    mTrigAbd_resp(g4:end,1),'Color','k','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 5]); ylabel('\int Abdominal (STD)');

% subplot(3,1,3); plot(timelag(1:g1),...
%     mTrigAbd_resp(1:g1,1),'Color','k','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('\int Abdominal (STD)');
% hold on
% plot(timelag(g1:g2),...
%     mTrigAbd_resp(g1:g2,1),'Color','b','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('\int Abdominal (STD)');
% plot(timelag(g2:g3),...
%     mTrigAbd_resp(g2:g3,1),'Color','r','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('\int Abdominal (STD)');
% plot(timelag(g3:g4),...
%     mTrigAbd_resp(g3:g4,1),'Color','c','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('\int Abdominal (STD)');
% plot(timelag(g4:end),...
%     mTrigAbd_resp(g4:end,1),'Color','k','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('\int Abdominal (STD)');

% hold on; plot([1 1].*((a_pos(1,2,1)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25],'--g','LineWidth',2);
% plot([1 1].*((a_pos(2,2,1)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25],'--g','LineWidth',2);
% plot([1 1].*((a_pos(2,2,2)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25],'--g','LineWidth',2);
% plot([1 1].*((a_pos(2,2,3)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25],'--g','LineWidth',2);

hold on; plot([1 1].*((a_pos(1,2,1)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 5],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,2,1)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 5],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,2,2)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 5],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,2,3)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 5],'--g','LineWidth',2);

xlim([-1 1]); xlabel('Lag (s)');
set(gca,'Box','off');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

g1 = a_pos(1,1,1);
g2 = a_pos(2,1,1);
g3 = a_pos(2,1,2);
g4 = a_pos(2,1,3);
figure
% plot(timelag,mTrigAir_resp(:,1),'Color',groupcolor,'LineWidth',2); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)');
subplot(3,1,1); plot(timelag(1:g1),...
    mTrigAir_base(1:g1,1),'Color','k','LineWidth',2); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)');
hold on
plot(timelag(g1:g2),...
    mTrigAir_base(g1:g2,1),'Color','b','LineWidth',2); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)');
plot(timelag(g2:g3),...
    mTrigAir_base(g2:g3,1),'Color','r','LineWidth',2); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)');
plot(timelag(g3:g4),...
    mTrigAir_base(g3:g4,1),'Color','c','LineWidth',2); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)');
plot(timelag(g4:end),...
    mTrigAir_base(g4:end,1),'Color','k','LineWidth',2); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)');

hold on; plot([1 1].*((a_pos(1,1,1)-half_triglen+1)./1000),[min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,1,1)-half_triglen+1)./1000),[min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,1,2)-half_triglen+1)./1000),[min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,1,3)-half_triglen+1)./1000),[min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25],'--g','LineWidth',2);
xlim([-1 1]); set(gca,'Box','off');

subplot(3,1,2); plot(timelag(1:g1),...
    mTrigDia_base(1:g1,1),'Color','k','LineWidth',2); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('\int Diaphragm (STD)');
hold on
plot(timelag(g1:g2),...
    mTrigDia_base(g1:g2,1),'Color','b','LineWidth',2); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('\int Diaphragm (STD)');
plot(timelag(g2:g3),...
    mTrigDia_base(g2:g3,1),'Color','r','LineWidth',2); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('\int Diaphragm (STD)');
plot(timelag(g3:g4),...
    mTrigDia_base(g3:g4,1),'Color','c','LineWidth',2); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('\int Diaphragm (STD)');
plot(timelag(g4:end),...
    mTrigDia_base(g4:end,1),'Color','k','LineWidth',2); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('\int Diaphragm (STD)');

hold on; plot([1 1].*((a_pos(1,1,1)-half_triglen+1)./1000),[min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,1,1)-half_triglen+1)./1000),[min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,1,2)-half_triglen+1)./1000),[min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,1,3)-half_triglen+1)./1000),[min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25],'--g','LineWidth',2);
xlim([-1 1]); set(gca,'Box','off');

subplot(3,1,3); plot(timelag(1:g1),...
    mTrigAbd_base(1:g1,1),'Color','k','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 5]); ylabel('\int Abdominal (STD)');
hold on
plot(timelag(g1:g2),...
    mTrigAbd_base(g1:g2,1),'Color','b','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 5]); ylabel('\int Abdominal (STD)');
plot(timelag(g2:g3),...
    mTrigAbd_base(g2:g3,1),'Color','r','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 5]); ylabel('\int Abdominal (STD)');
plot(timelag(g3:g4),...
    mTrigAbd_base(g3:g4,1),'Color','c','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 5]); ylabel('\int Abdominal (STD)');
plot(timelag(g4:end),...
    mTrigAbd_base(g4:end,1),'Color','k','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 5]); ylabel('\int Abdominal (STD)');

% subplot(3,1,3); plot(timelag(1:g1),...
%     mTrigAbd_resp(1:g1,1),'Color','k','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('\int Abdominal (STD)');
% hold on
% plot(timelag(g1:g2),...
%     mTrigAbd_resp(g1:g2,1),'Color','b','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('\int Abdominal (STD)');
% plot(timelag(g2:g3),...
%     mTrigAbd_resp(g2:g3,1),'Color','r','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('\int Abdominal (STD)');
% plot(timelag(g3:g4),...
%     mTrigAbd_resp(g3:g4,1),'Color','c','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('\int Abdominal (STD)');
% plot(timelag(g4:end),...
%     mTrigAbd_resp(g4:end,1),'Color','k','LineWidth',2); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('\int Abdominal (STD)');

% hold on; plot([1 1].*((a_pos(1,2,1)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25],'--g','LineWidth',2);
% plot([1 1].*((a_pos(2,2,1)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25],'--g','LineWidth',2);
% plot([1 1].*((a_pos(2,2,2)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25],'--g','LineWidth',2);
% plot([1 1].*((a_pos(2,2,3)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25],'--g','LineWidth',2);

hold on; plot([1 1].*((a_pos(1,1,1)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 5],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,1,1)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 5],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,1,2)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 5],'--g','LineWidth',2);
plot([1 1].*((a_pos(2,1,3)-half_triglen+1)./1000),[min(TrigAbd,[],'all')-0.25 5],'--g','LineWidth',2);

xlim([-1 1]); xlabel('Lag (s)');
set(gca,'Box','off');
%% 3d explanation
g1 = a_pos(1,2,1);
g2 = a_pos(2,2,1);
g3 = a_pos(2,2,2);
g4 = a_pos(2,2,3);
gB1 = a_pos(1,1,1);
gB2 = a_pos(2,1,1);
gB3 = a_pos(2,1,2);
gB4 = a_pos(2,1,3);
figure
plot3(mTrigAbd_resp(1300:g1,1),mTrigAir_resp(1300:g1,1),mTrigDia_resp(1300:g1,1),'Color','k','LineWidth',lwidth);
hold on;
plot3(mTrigAbd_resp(g1:g2,1),mTrigAir_resp(g1:g2,1),mTrigDia_resp(g1:g2,1),'Color','b','LineWidth',lwidth);
plot3(mTrigAbd_resp(g2:g3,1),mTrigAir_resp(g2:g3,1),mTrigDia_resp(g2:g3,1),'Color','r','LineWidth',lwidth);
plot3(mTrigAbd_resp(g3:g4,1),mTrigAir_resp(g3:g4,1),mTrigDia_resp(g3:g4,1),'Color','c','LineWidth',lwidth);
plot3(mTrigAbd_resp(g4:3000,1),mTrigAir_resp(g4:3000,1),mTrigDia_resp(g4:3000,1),'Color','k','LineWidth',lwidth);

plot3(mTrigAbd_base(1500:gB1,1),mTrigAir_base(1500:gB1,1),mTrigDia_base(1500:gB1,1),'LineStyle','--','Color','k','LineWidth',lwidth);
hold on;
plot3(mTrigAbd_base(gB1:gB2,1),mTrigAir_base(gB1:gB2,1),mTrigDia_base(gB1:gB2,1),'LineStyle','--','Color','b','LineWidth',lwidth);
plot3(mTrigAbd_base(gB2:gB3,1),mTrigAir_base(gB2:gB3,1),mTrigDia_base(gB2:gB3,1),'LineStyle','--','Color','r','LineWidth',lwidth);
plot3(mTrigAbd_base(gB3:gB4,1),mTrigAir_base(gB3:gB4,1),mTrigDia_base(gB3:gB4,1),'LineStyle','--','Color','c','LineWidth',lwidth);
plot3(mTrigAbd_base(gB4:3000,1),mTrigAir_base(gB4:3000,1),mTrigDia_base(gB4:3000,1),'LineStyle','--','Color','k','LineWidth',lwidth);


plot3(linspace(-5,12,2),[0 0],[0 0],'--k','LineWidth',2); hold on;
plot3([0 0],linspace(-3,4,2),[0 0],'--k','LineWidth',2); hold on;
plot3([0,0],[0,0],linspace(-1,5,2),'--k','LineWidth',2); hold on;
xlim([-5 12]); ylim([-3 4]); zlim([-1 5]);
set(gca,'ydir','reverse');  
% view([-150,30]);
view([-40,20]);
xlabel('\int Abdominal (SD)'); ylabel('Airflow (SD)'); zlabel('\int Diaphragm (SD)');
set(gca,'Box','off');
figure
plot(mTrigAir_resp(1300:g1,1),mTrigDia_resp(1300:g1,1),'Color','k','Linewidth',lwidth); hold on;
plot(mTrigAir_resp(g1:g2,1),mTrigDia_resp(g1:g2,1),'Color','b','Linewidth',lwidth); 
plot(mTrigAir_resp(g2:g3,1),mTrigDia_resp(g2:g3,1),'Color','r','Linewidth',lwidth);
plot(mTrigAir_resp(g3:g4,1),mTrigDia_resp(g3:g4,1),'Color','c','Linewidth',lwidth);
plot(mTrigAir_resp(g4:3000,1),mTrigDia_resp(g4:3000,1),'Color','k','Linewidth',lwidth);

plot(mTrigAir_base(1500:gB1,1),mTrigDia_base(1500:gB1,1),'LineStyle','--','Color','k','Linewidth',lwidth); hold on;
plot(mTrigAir_base(gB1:gB2,1),mTrigDia_base(gB1:gB2,1),'LineStyle','--','Color','b','Linewidth',lwidth); 
plot(mTrigAir_base(gB2:gB3,1),mTrigDia_base(gB2:gB3,1),'LineStyle','--','Color','r','Linewidth',lwidth);
plot(mTrigAir_base(gB3:gB4,1),mTrigDia_base(gB3:gB4,1),'LineStyle','--','Color','c','Linewidth',lwidth);
plot(mTrigAir_base(gB4:3000,1),mTrigDia_base(gB4:3000,1),'LineStyle','--','Color','k','Linewidth',lwidth);

hold on;
plot(linspace(-3,4,2),[0,0],'--k','LineWidth',2); hold on;
plot([0,0],linspace(-1,5,2),'--k','LineWidth',2); 
xlim([-3 4 ]); ylim([-1 5]); xlabel('Airflow (SD)'); ylabel('\int Diaphragm (SD)');
set(gca,'Box','off');

figure
plot(mTrigAir_resp(1300:g1,1),mTrigAbd_resp(1300:g1,1),'Color','k','Linewidth',lwidth); hold on;
plot(mTrigAir_resp(g1:g2,1),mTrigAbd_resp(g1:g2,1),'Color','b','Linewidth',lwidth);
plot(mTrigAir_resp(g2:g3,1),mTrigAbd_resp(g2:g3,1),'Color','r','Linewidth',lwidth);
plot(mTrigAir_resp(g3:g4,1),mTrigAbd_resp(g3:g4,1),'Color','c','Linewidth',lwidth);
plot(mTrigAir_resp(g4:3000,1),mTrigAbd_resp(g4:3000,1),'Color','k','Linewidth',lwidth);

plot(mTrigAir_base(1500:gB1,1),mTrigAbd_base(1500:gB1,1),'LineStyle','--','Color','k','Linewidth',lwidth); hold on;
plot(mTrigAir_base(gB1:gB2,1),mTrigAbd_base(gB1:gB2,1),'LineStyle','--','Color','b','Linewidth',lwidth);
plot(mTrigAir_base(gB2:gB3,1),mTrigAbd_base(gB2:gB3,1),'LineStyle','--','Color','r','Linewidth',lwidth);
plot(mTrigAir_base(gB3:gB4,1),mTrigAbd_base(gB3:gB4,1),'LineStyle','--','Color','c','Linewidth',lwidth);
plot(mTrigAir_base(gB4:3000,1),mTrigAbd_base(gB4:3000,1),'LineStyle','--','Color','k','Linewidth',lwidth);

plot(linspace(-3,4,2),[0,0],'--k','LineWidth',2); hold on;
plot([0,0],linspace(-5,12,2),'--k','LineWidth',2); hold on;
xlim([-3 4]); ylim([-5 12]); xlabel('Airflow (SD)'); ylabel('\int Abdominal (SD)');
set(gca,'Box','off');

% fig1 = 1; fig2 = 10; %1 is the first response
% figure
% ax1 = subplot(3,2,1); shadedErrorBar(timelag,mTrigAir_resp(:,fig1),errTrigAir_resp(:,fig1),'lineprops',groupcolor); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)'); xlim([-1 1]); set(gca,'Box','off');
% hold on; 
% % scatter((a_pos(:,fig1,1)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig1,2)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig1,3)-half_triglen+1)./1000,[0 0],'g');
% title(strcat(num2str((fig1-1)*bindur),'-',num2str((fig1)*bindur),'min Post-Injection'));
% ax2 = subplot(3,2,2); shadedErrorBar(timelag,mTrigAir_resp(:,fig2),errTrigAir_resp(:,fig2),'lineprops',groupcolor); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); xlim([-1 1]); set(gca,'Box','off');
% % hold on; shadedErrorBar(timelag,mTrigDia_resp(:,6),errTrigDia_resp(:,6),'lineprops','k');
% title(strcat(num2str((fig2-1)*bindur),'-',num2str((fig2)*bindur),'min Post-Injection'));
% hold on; 
% % scatter((a_pos(:,fig2,1)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig2,2)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig2,3)-half_triglen+1)./1000,[0 0],'g');
% ax3 = subplot(3,2,3); shadedErrorBar(timelag,mTrigDia_resp(:,fig1),errTrigDia_resp(:,fig1),'lineprops',groupcolor); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('\int Diaphragm (STD)'); xlim([-1 1]); set(gca,'Box','off');
% hold on; 
% % scatter((a_pos(:,fig1,1)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig1,2)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig1,3)-half_triglen+1)./1000,[0 0],'g');
% ax4 = subplot(3,2,4); shadedErrorBar(timelag,mTrigDia_resp(:,fig2),errTrigDia_resp(:,fig2),'lineprops',groupcolor); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); xlim([-1 1]); set(gca,'Box','off');
% % hold on; shadedErrorBar(timelag,mTrigAir_resp(:,6),errTrigAir_resp(:,6),'lineprops','m');
% hold on; 
% % scatter((a_pos(:,fig2,1)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig2,2)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig2,3)-half_triglen+1)./1000,[0 0],'g');
% ax5 = subplot(3,2,5); shadedErrorBar(timelag,mTrigAbd_resp(:,fig1),errTrigAbd_resp(:,fig1),'lineprops',groupcolor); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('\int Abdominal (STD)'); xlim([-1 1]); set(gca,'Box','off');
% hold  on; xlabel('Lag (s)');
% % scatter((a_pos(:,fig1,1)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig1,2)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig1,3)-half_triglen+1)./1000,[0 0],'g');
% ax6 = subplot(3,2,6); shadedErrorBar(timelag,mTrigAbd_resp(:,fig2),errTrigAbd_resp(:,fig2),'lineprops',groupcolor); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); xlim([-1 1]); set(gca,'Box','off');
% hold on; xlabel('Lag (s)');
% % scatter((a_pos(:,fig2,1)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig2,2)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig2,3)-half_triglen+1)./1000,[0 0],'g');
% % hold on; shadedErrorBar(timelag,mTrigAbd_resp(:,6),errTrigAbd_resp(:,6),'lineprops','r');
% linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x'); 

%% Triggered cycles
fig1 = 1; fig2 = 10; %1 is the first response
figure
ax1 = subplot(3,2,1); shadedErrorBar(timelag,mTrigAir_resp(:,fig1),errTrigAir_resp(:,fig1),'lineprops',groupcolor); ylim([-3 4]); ylabel('Airflow (STD)'); xlim([-1 1]); set(gca,'Box','off');
hold on; 
% scatter((a_pos(:,fig1,1)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig1,2)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig1,3)-half_triglen+1)./1000,[0 0],'g');
title(strcat(num2str((fig1-1)*bindur),'-',num2str((fig1)*bindur),'min Post-Injection'));
ax2 = subplot(3,2,2); shadedErrorBar(timelag,mTrigAir_resp(:,fig2),errTrigAir_resp(:,fig2),'lineprops',groupcolor); ylim([-3 4]); xlim([-1 1]); set(gca,'Box','off');
% hold on; shadedErrorBar(timelag,mTrigDia_resp(:,6),errTrigDia_resp(:,6),'lineprops','k');
title(strcat(num2str((fig2-1)*bindur),'-',num2str((fig2)*bindur),'min Post-Injection'));
hold on; 
% scatter((a_pos(:,fig2,1)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig2,2)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig2,3)-half_triglen+1)./1000,[0 0],'g');
ax3 = subplot(3,2,3); shadedErrorBar(timelag,mTrigDia_resp(:,fig1),errTrigDia_resp(:,fig1),'lineprops',groupcolor); ylim([-1 5]); ylabel('\int Diaphragm (STD)'); xlim([-1 1]); set(gca,'Box','off');
hold on; 
% scatter((a_pos(:,fig1,1)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig1,2)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig1,3)-half_triglen+1)./1000,[0 0],'g');
ax4 = subplot(3,2,4); shadedErrorBar(timelag,mTrigDia_resp(:,fig2),errTrigDia_resp(:,fig2),'lineprops',groupcolor); ylim([-1 5]); xlim([-1 1]); set(gca,'Box','off');
% hold on; shadedErrorBar(timelag,mTrigAir_resp(:,6),errTrigAir_resp(:,6),'lineprops','m');
hold on; 
% scatter((a_pos(:,fig2,1)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig2,2)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig2,3)-half_triglen+1)./1000,[0 0],'g');
ax5 = subplot(3,2,5); shadedErrorBar(timelag,mTrigAbd_resp(:,fig1),errTrigAbd_resp(:,fig1),'lineprops',groupcolor); ylim([-5 12]); ylabel('\int Abdominal (STD)'); xlim([-1 1]); set(gca,'Box','off');
hold on; xlabel('Lag (s)');
% scatter((a_pos(:,fig1,1)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig1,2)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig1,3)-half_triglen+1)./1000,[0 0],'g');
ax6 = subplot(3,2,6); shadedErrorBar(timelag,mTrigAbd_resp(:,fig2),errTrigAbd_resp(:,fig2),'lineprops',groupcolor); ylim([-5 12]); xlim([-1 1]); set(gca,'Box','off');
hold on; xlabel('Lag (s)');
% scatter((a_pos(:,fig2,1)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig2,2)-half_triglen+1)./1000,[0 0],'g'); scatter((a_pos(:,fig2,3)-half_triglen+1)./1000,[0 0],'g');
% hold on; shadedErrorBar(timelag,mTrigAbd_resp(:,6),errTrigAbd_resp(:,6),'lineprops','r');
linkaxes([ax1 ax2 ax3 ax4 ax5 ax6],'x'); 


 
%Overlays;
% binstart1 = find(strt_ind>(((fig1-1)*bindur)*60*1000)+postinj_period,1,'first');   
% binend1 = find(strt_ind<(fig1*bindur*60*1000)+postinj_period,1,'last');
% binstart2 = find(strt_ind>(((fig2-1)*bindur)*60*1000)+postinj_period,1,'first');   
% binend2 = find(strt_ind<(fig2*bindur*60*1000)+postinj_period,1,'last');
% figure
% subplot(3,2,1); plot(timelag,TrigDia(:,binstart1:binend1),'k'); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('STD')
% subplot(3,2,2); plot(timelag,TrigDia(:,binstart2:binend2),'k'); ylim([min(TrigDia,[],'all')-0.25 max(TrigDia,[],'all')+0.25]); ylabel('STD')
% subplot(3,2,3); plot(timelag,TrigAir(:,binstart1:binend1),'m'); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('STD')
% subplot(3,2,4); plot(timelag,TrigAir(:,binstart2:binend2),'m'); ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('STD')
% subplot(3,2,5); plot(timelag,TrigAbd(:,binstart1:binend1),'r'); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('STD')
% subplot(3,2,6); plot(timelag,TrigAbd(:,binstart2:binend2),'r'); ylim([min(TrigAbd,[],'all')-0.25 max(TrigAbd,[],'all')+0.25]); ylabel('STD')

%% Differences Explanatory Figures
c_strt = half_triglen-round((0.4*m_resprate(1,1)*1000),0);
c_fin = c_strt+round((m_resprate(fig1+1,1)*1000),0);
c_len = c_fin-c_strt+1;
in_strtR = a_pos(1,fig1+1,2); in_finR = a_pos(2,fig1+1,2); %zooming into the inspiratory period response

in_strt = a_pos(1,1,2); in_fin = a_pos(2,1,2); %zooming into the inspiratory period during baseline
% in_len = in_fin-in_strt+1;
% in_len2 = in_fin2-in_strt2+1;

lte_strt_R = a_pos(1,fig1+1,1); lte_fin_R = a_pos(2,fig1+1,1); 
lte_strt_B = a_pos(1,1,1); lte_fin_B = a_pos(2,1,1); 
pstI_strt_R = a_pos(1,fig1+1,3); pstI_fin_R = a_pos(2,fig1+1,3); 
pstI_strt_B = a_pos(1,1,3); pstI_fin_B = a_pos(2,1,3); 

x = linspace(0,c_fin-c_strt,c_len)'./1000;
figure
plot(x,mTrigAir_resp(c_strt:c_fin,fig1),'Color','k','LineWidth',2); 
hold on; plot(x,mTrigAir_base(c_strt:c_fin,1),'LineStyle','--','Color','k','LineWidth',2);
ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (SD)'); xlabel('Lag (s)');

% fill([x(in_strt-c_strt:in_fin-c_strt)' flip(x(in_strt2-c_strt:in_fin2-c_strt))'],...
%     [mTrigAir_base(in_strt:in_fin,1)' flip(mTrigAir_resp(in_strt2:in_fin2,fig1))'],'g','LineStyle','none');
% set(gca,'Box','off');

%baseline in green, response in red
fill([x(lte_strt_B-c_strt:lte_fin_B-c_strt)' flip(x(lte_strt_B-c_strt:lte_fin_B-c_strt))'],...
    [mTrigAir_base(lte_strt_B:lte_fin_B,1)' zeros(1,length(x(lte_strt_B-c_strt:lte_fin_B-c_strt)))],...
    'g','LineStyle','none','FaceAlpha',0.5);
hold on
fill([x(lte_strt_R-c_strt:lte_fin_R-c_strt)' flip(x(lte_strt_R-c_strt:lte_fin_R-c_strt))'],...
    [mTrigAir_resp(lte_strt_R:lte_fin_R,1)' zeros(1,length(x(lte_strt_R-c_strt:lte_fin_R-c_strt)))],...
    'b','LineStyle','none','FaceAlpha',0.5);

%baseline in green, response in red
fill([x(in_strt-c_strt:in_fin-c_strt)' flip(x(in_strt-c_strt:in_fin-c_strt))'],...
    [mTrigAir_base(in_strt:in_fin,1)' zeros(1,length(x(in_strt-c_strt:in_fin-c_strt)))],...
    'g','LineStyle','none','FaceAlpha',0.5);
hold on
fill([x(in_strtR-c_strt:in_finR-c_strt)' flip(x(in_strtR-c_strt:in_finR-c_strt))'],...
    [mTrigAir_resp(in_strtR:in_finR,1)' zeros(1,length(x(in_strtR-c_strt:in_finR-c_strt)))],...
    'r','LineStyle','none','FaceAlpha',0.5);

%baseline in green, response in red
fill([x(pstI_strt_B-c_strt:pstI_fin_B-c_strt)' flip(x(pstI_strt_B-c_strt:pstI_fin_B-c_strt))'],...
    [mTrigAir_base(pstI_strt_B:pstI_fin_B,1)' zeros(1,length(x(pstI_strt_B-c_strt:pstI_fin_B-c_strt)))],...
    'g','LineStyle','none','FaceAlpha',0.5);
hold on
fill([x(pstI_strt_R-c_strt:pstI_fin_R-c_strt)' flip(x(pstI_strt_R-c_strt:pstI_fin_R-c_strt))'],...
    [mTrigAir_resp(pstI_strt_R:pstI_fin_R,1)' zeros(1,length(x(pstI_strt_R-c_strt:pstI_fin_R-c_strt)))],...
    'c','LineStyle','none','FaceAlpha',0.5);
set(gca,'Box','off');
% xlim([-1 1]);

% figure
% % subplot(1,2,1); plot(x(in_strt-c_strt:in_fin-c_strt),mTrigAir_base(in_strt:in_fin,1),'LineStyle','--','Color',groupcolor,'LineWidth',2); hold on;
% % ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)'); xlabel('Lag (s)');
% % plot(x(in_strt2-c_strt:in_fin2-c_strt),mTrigAir_resp(in_strt2:in_fin2,fig1),'Color',groupcolor,'LineWidth',2);
% % % fill([x(in_strt:in_fin)' flip(x(in_strt:in_fin))'],...
% % %     [mTrigAir_base(in_strt:in_fin,1)' flip(mTrigAir_resp(in_strt:in_fin,fig1))'],'g');
% % fill([x(in_strt-c_strt:in_fin-c_strt)' flip(x(in_strt2-c_strt:in_fin2-c_strt))'],...
% %     [mTrigAir_base(in_strt:in_fin,1)' flip(mTrigAir_resp(in_strt2:in_fin2,fig1))'],'g','LineStyle','none');
% % set(gca,'Box','off');
% 
% subplot(1,2,1);
% plot(x(in_strt-c_strt:in_fin-c_strt),mTrigAir_base(in_strt:in_fin,1),'LineStyle','--','Color',groupcolor,'LineWidth',2); hold on;
% ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)'); xlabel('Lag (s)');
% plot(x(in_strt-c_strt:in_fin-c_strt),mTrigAir_resp(in_strt:in_fin,fig1),'Color',groupcolor,'LineWidth',2);
% % fill([x(in_strt:in_fin)' flip(x(in_strt:in_fin))'],...
% %     [mTrigAir_base(in_strt:in_fin,1)' flip(mTrigAir_resp(in_strt:in_fin,fig1))'],'g');
% fill([x(in_strt-c_strt:in_fin-c_strt)' flip(x(in_strt-c_strt:in_fin-c_strt))'],...
%     [mTrigAir_base(in_strt:in_fin,1)' flip(mTrigAir_resp(in_strt:in_fin,fig1))'],'g','LineStyle','none');
% set(gca,'Box','off');
% 
% % timevec*(n1/n2);
% timevec = linspace(0,in_len2-1,in_len2)'./1000; timevec = timevec*((in_len)/(in_len2));
% timevec = timevec +((in_strt-c_strt)./1000);
% subplot(1,2,2); plot(x(in_strt-c_strt:in_fin-c_strt),mTrigAir_base(in_strt:in_fin,1),'LineStyle','--','Color',groupcolor,'LineWidth',2); hold on;
% ylim([min(TrigAir,[],'all')-0.25 max(TrigAir,[],'all')+0.25]); ylabel('Airflow (STD)'); xlabel('Lag (s)');
% plot(timevec,mTrigAir_resp(in_strt2:in_fin2,fig1),'Color',groupcolor,'LineWidth',2);
% % fill([x(in_strt2:in_fin2)' flip(x(in_strt2:in_fin2))'],...
% %     [mTrigAir_base(in_strt2:in_fin2,1)' flip(mTrigAir_resp(in_strt2:in_fin2,fig1))'],'g');
% fill([x(in_strt-c_strt:in_fin-c_strt)' flip(timevec)'],...
%     [mTrigAir_base(in_strt:in_fin,1)' flip(mTrigAir_resp(in_strt2:in_fin2,fig1))'],'g','LineStyle','none');
% set(gca,'Box','off');

%Differences 
% figure
% subplot(3,2,1); shadedErrorBar(timelag,mTrigDiaresp_mbase(:,fig1),errTrigDiaresp_mbase(:,fig1),'lineprops','k'); ylim([min(TrigDia_minusbase,[],'all') max(TrigDia_minusbase,[],'all')]); ylabel('STD')
% title(strcat(num2str((fig1-1)*bindur),'-',num2str((fig1)*bindur),'min Post-Injection'));
% subplot(3,2,2); shadedErrorBar(timelag,mTrigDiaresp_mbase(:,fig2),errTrigDiaresp_mbase(:,fig2),'lineprops','k'); ylim([min(TrigDia_minusbase,[],'all') max(TrigDia_minusbase,[],'all')]);
% % hold on; shadedErrorBar(timelag,mTrigDiaresp_mbase(:,6),errTrigDiaresp_mbase(:,6),'lineprops','k');
% title(strcat(num2str((fig2-1)*bindur),'-',num2str((fig2)*bindur),'min Post-Injection'));
% 
% subplot(3,2,3); shadedErrorBar(timelag,mTrigAirresp_mbase(:,fig1),errTrigAirresp_mbase(:,fig1),'lineprops','m'); ylim([min(TrigAir_minusbase,[],'all') max(TrigAir_minusbase,[],'all')]); ylabel('STD')
% subplot(3,2,4); shadedErrorBar(timelag,mTrigAirresp_mbase(:,fig2),errTrigAirresp_mbase(:,fig2),'lineprops','m'); ylim([min(TrigAir_minusbase,[],'all') max(TrigAir_minusbase,[],'all')]);
% % hold on; shadedErrorBar(timelag,mTrigAirresp_mbase(:,6),errTrigAirresp_mbase(:,6),'lineprops','m');
% 
% subplot(3,2,5); shadedErrorBar(timelag,mTrigAbdresp_mbase(:,fig1),errTrigAbdresp_mbase(:,fig1),'lineprops','r'); ylim([min(TrigAbd_minusbase,[],'all') max(TrigAbd_minusbase,[],'all')]); ylabel('STD')
% subplot(3,2,6); shadedErrorBar(timelag,mTrigAbdresp_mbase(:,fig2),errTrigAbdresp_mbase(:,fig2),'lineprops','r'); ylim([min(TrigAbd_minusbase,[],'all') max(TrigAbd_minusbase,[],'all')]);
% % hold on; shadedErrorBar(timelag,mTrigAbdresp_mbase(:,6),errTrigAbdresp_mbase(:,6),'lineprops','r');

% %Euclidean Distance
% figure
% subplot(3,3,[1 4]); shadedErrorBar(timelag,mEucDist_base,errEucDist_base,'lineprops','k'); ylim([min(EucDistcycle,[],'all') max(EucDistcycle,[],'all')]); ylabel('Distance');
% subplot(3,3,[2 5]); shadedErrorBar(timelag,mEucDist_resp(:,fig1),errEucDist_resp(:,fig1),'lineprops','k'); ylim([min(EucDistcycle,[],'all') max(EucDistcycle,[],'all')]); 
% title(strcat(num2str((fig1-1)*bindur),'-',num2str((fig1)*bindur),'min Post-Injection'));
% subplot(3,3,[3 6]); shadedErrorBar(timelag,mEucDist_resp(:,fig2),errEucDist_resp(:,fig2),'lineprops','k'); ylim([min(EucDistcycle,[],'all') max(EucDistcycle,[],'all')]);
% % hold on; shadedErrorBar(timelag,mTrigDiaresp_mbase(:,6),errTrigDiaresp_mbase(:,6),'lineprops','k');
% title(strcat(num2str((fig2-1)*bindur),'-',num2str((fig2)*bindur),'min Post-Injection'));
% 
% subplot(3,3,7); shadedErrorBar(timelag,mTrigAir_resp(:,fig1),errTrigAir_resp(:,fig1),'lineprops','m'); ylim([min(TrigAir,[],'all') max(TrigAir,[],'all')]); ylabel('STD')
% subplot(3,3,8); shadedErrorBar(timelag,mTrigAir_resp(:,fig1),errTrigAir_resp(:,fig1),'lineprops','m'); ylim([min(TrigAir,[],'all') max(TrigAir,[],'all')]); 
% subplot(3,3,9); shadedErrorBar(timelag,mTrigAir_resp(:,fig1),errTrigAir_resp(:,fig1),'lineprops','m'); ylim([min(TrigAir,[],'all') max(TrigAir,[],'all')]);
% 
% %Mahal Dist
% % figure
% % subplot(3,3,[1 4]); plot(timelag,MahDist(:,1:preinj_period),'k'); hold on; plot(timelag,ones(length(timelag),1).*chi2inv(0.95,2),'--g');
% % ylim([min(MahDist,[],'all') max(MahDist,[],'all')]); ylabel('Distance')
% % subplot(3,3,[2 5]); plot(timelag,MahDist(:,binstart1:binend1),'k'); hold on; plot(timelag,ones(length(timelag),1).*chi2inv(0.95,2),'--g');
% % ylim([min(MahDist,[],'all') max(MahDist,[],'all')]);
% % title(strcat(num2str((fig1-1)*bindur),'-',num2str((fig1)*bindur),'min Post-Injection'));
% % subplot(3,3,[3 6]); plot(timelag,MahDist(:,binstart2:binend2),'k'); hold on; plot(timelag,ones(length(timelag),1).*chi2inv(0.95,2),'--g');
% % ylim([min(MahDist,[],'all') max(MahDist,[],'all')]);
% % title(strcat(num2str((fig2-1)*bindur),'-',num2str((fig2)*bindur),'min Post-Injection'));
% % 
% % subplot(3,3,7); shadedErrorBar(timelag,mTrigAir_resp(:,fig1),errTrigAir_resp(:,fig1),'lineprops','m'); ylim([min(TrigAir,[],'all') max(TrigAir,[],'all')]); ylabel('STD')
% % subplot(3,3,8); shadedErrorBar(timelag,mTrigAir_resp(:,fig1),errTrigAir_resp(:,fig1),'lineprops','m'); ylim([min(TrigAir,[],'all') max(TrigAir,[],'all')]); 
% % subplot(3,3,9); shadedErrorBar(timelag,mTrigAir_resp(:,fig1),errTrigAir_resp(:,fig1),'lineprops','m'); ylim([min(TrigAir,[],'all') max(TrigAir,[],'all')]);
% 
% figure
% subplot(3,3,[1 4]); plot(timelag,mMahDist_base,'k'); hold on; plot(timelag,ones(length(timelag),1).*chi2inv(0.95,2),'--g');
% ylim([min(mMahDist_resp+2,[],'all') max(mMahDist_resp+2,[],'all')]); ylabel('Distance')
% subplot(3,3,[2 5]); plot(timelag,mMahDist_resp(:,fig1),'k'); hold on; plot(timelag,ones(length(timelag),1).*chi2inv(0.95,2),'--g');
% ylim([min(mMahDist_resp+2,[],'all') max(mMahDist_resp+2,[],'all')]);
% title(strcat(num2str((fig1-1)*bindur),'-',num2str((fig1)*bindur),'min Post-Injection'));
% subplot(3,3,[3 6]); plot(timelag,mMahDist_resp(:,fig2),'k'); hold on; plot(timelag,ones(length(timelag),1).*chi2inv(0.95,2),'--g');
% ylim([min(mMahDist_resp+2,[],'all') max(mMahDist_resp+2,[],'all')]);
% title(strcat(num2str((fig2-1)*bindur),'-',num2str((fig2)*bindur),'min Post-Injection'));
% 
% subplot(3,3,7); shadedErrorBar(timelag,mTrigAir_resp(:,fig1),errTrigAir_resp(:,fig1),'lineprops','m'); ylim([min(TrigAir,[],'all') max(TrigAir,[],'all')]); ylabel('STD')
% subplot(3,3,8); shadedErrorBar(timelag,mTrigAir_resp(:,fig1),errTrigAir_resp(:,fig1),'lineprops','m'); ylim([min(TrigAir,[],'all') max(TrigAir,[],'all')]); 
% subplot(3,3,9); shadedErrorBar(timelag,mTrigAir_resp(:,fig1),errTrigAir_resp(:,fig1),'lineprops','m'); ylim([min(TrigAir,[],'all') max(TrigAir,[],'all')]);

%% Save
file = uigetfile;
load(file);
% fname = INFO.Name;

fnum = input('Filenumber=');
fname = strcat('file',int2str(fnum));
% annettesummary.(fname).fnumber = fnum;
annettesummary.(fname).params = struct('Triglength',half_triglen*2,'nBinsPostInj',postinj_nBins);
annettesummary.periods{1,fnum} = a_pos;

% % annettesummary.diadist{1,end+1} = diadist; annettesummary.n_diadist{1,end+1} = n_diadist;
% % annettesummary.airdist{1,end+1} = airdist; annettesummary.n_airdist{1,end+1} = n_airdist;
% % annettesummary.abdist{1,end+1} = abddist; annettesummary.n_abdist{1,end+1} = n_abddist;
% % annettesummary.eucdist{1,end+1} = eucdist; annettesummary.n_eucdist{1,end+1} = n_eucdist;
% % annettesummary.mahaldist{1,end+1} = mahaldist; annettesummary.n_mahaldist{1,end+1} = n_mahaldist;
% % annettesummary.rrlen{1,end+1} = relringlen; annettesummary.n_rrlen{1,end+1} = n_relringlen;
% % annettesummary.respperiod{1,end+1} = m_resprate;

annettesummary.diadist{1,fnum} = rawdiadist; annettesummary.n_diadist{1,fnum} = n2_diadist;
% annettesummary.rawdisdist{1,fnum} = 
annettesummary.airdist{1,fnum} = rawairdist; annettesummary.n_airdist{1,fnum} = n2_airdist;
annettesummary.abdist{1,fnum} = rawabddist; annettesummary.n_abdist{1,fnum} = n2_abddist;

annettesummary.eucdist{1,fnum} = eucdist; annettesummary.n_eucdist{1,fnum} = n2_eucdist;
annettesummary.mahaldist{1,fnum} = mahaldist; annettesummary.n_mahaldist{1,fnum} = n2_mahaldist;
% annettesummary.rrlen{1,fnum} = rawrelringlen; annettesummary.n_rrlen{1,fnum} = n_relringlen;
annettesummary.rrlen{1,fnum} = nan(10,3); annettesummary.n_rrlen{1,fnum} = nan(10,3);

annettesummary.respperiod{1,fnum} = m_resprate;

annettesummary.mTrigDia_base{1,fnum} = mTrigDia_base; annettesummary.eTrigDia_base{1,fnum} = errTrigDia_base;
annettesummary.mTrigAir_base{1,fnum} = mTrigAir_base; annettesummary.eTrigAir_base{1,fnum} = errTrigAir_base;
annettesummary.mTrigAbd_base{1,fnum} = mTrigAbd_base; annettesummary.eTrigAbd_base{1,fnum} = errTrigAbd_base;
annettesummary.mTrigDia_resp{1,fnum} = mTrigDia_resp; annettesummary.eTrigDia_resp{1,fnum} = errTrigDia_resp;
annettesummary.mTrigAir_resp{1,fnum} = mTrigAir_resp; annettesummary.eTrigAir_resp{1,fnum} = errTrigAir_resp;
annettesummary.mTrigAbd_resp{1,fnum} = mTrigAbd_resp; annettesummary.eTrigAbd_resp{1,fnum} = errTrigAbd_resp;
save('Annette2023.mat','annettesummary','-v7.3')
