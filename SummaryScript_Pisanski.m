%Computes averages and stats across the whole dataset
%Requires the data structure from the other mains script for this project.

%Written by Mitchell Prostebby, 2023
nRec = 32;
for i = 1:33
% for i = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 28 29 30 31 32 33]
    rsp(:,i) = annettesummary.respperiod{1,i};
    % rsp(:,i) = rsp(:,i)./rsp(1,i); %normalized to baseline
    for j = 1:3
        ddd(:,i,j) = annettesummary.diadist{1,i}(:,j);
        arr(:,i,j) = annettesummary.airdist{1,i}(:,j);
        abb(:,i,j) = annettesummary.abdist{1,i}(:,j);
        euu(:,i,j) = annettesummary.eucdist{1,i}(:,j);
        mmm(:,i,j) = annettesummary.mahaldist{1,i}(:,j);
        rrr(:,i,j) = annettesummary.rrlen{1,i}(:,j);

        nddd(:,i,j) = annettesummary.n_diadist{1,i}(:,j);
        narr(:,i,j) = annettesummary.n_airdist{1,i}(:,j);
        nabb(:,i,j) = annettesummary.n_abdist{1,i}(:,j);
        neuu(:,i,j) = annettesummary.n_eucdist{1,i}(:,j);
        nmmm(:,i,j) = annettesummary.n_mahaldist{1,i}(:,j);
        nrrr(:,i,j) = annettesummary.n_rrlen{1,i}(:,j);

        periods(i,:,j) = annettesummary.periods{1,i}(2,:,j)-annettesummary.periods{1,i}(1,:,j);
        p_rat(:,i,j) = (periods(i,:,j)./1000)./(rsp(:,i)');
    end
        
    t_resp(i,:,1) = (annettesummary.periods{1,i}(2,:,2)-annettesummary.periods{1,i}(1,:,2))./1000;
    t_resp(i,:,2) = (rsp(:,i)')-t_resp(i,:,1); 
end
sum_per = sum(periods,3); %has animals across the rows, and elapsed time across the columns 
diff_per = rsp'-(sum_per./1000);

%Removing experiments from experiments with errors
nG = 6;
g = [1 6;7 11;12 16;17 21;22 28;29 33];
rmve_ind = [11 26];
ddd(:,rmve_ind,:) = NaN;
arr(:,rmve_ind,:) = NaN;
abb(:,rmve_ind,:) = NaN;
euu(:,rmve_ind,:) = NaN;
mmm(:,rmve_ind,:) = NaN;
rrr(:,rmve_ind,:) = NaN;
nddd(:,rmve_ind,:) = NaN;
narr(:,rmve_ind,:) = NaN;
nabb(:,rmve_ind,:) = NaN;
neuu(:,rmve_ind,:) = NaN;
nmmm(:,rmve_ind,:) = NaN;
nrrr(:,rmve_ind,:) = NaN;
periods(rmve_ind,:,:) = NaN;
p_rat(:,rmve_ind,:) = NaN;
rsp(:,rmve_ind) = NaN;
t_resp(rmve_ind,:,1) = NaN;
t_resp(rmve_ind,:,2) = NaN;

for k = 1:nG
    for j = 1:3
        m_ddd(:,j,k) = mean(ddd(:,g(k,1):g(k,2),j),2,'omitnan'); e_ddd(:,j,k) = std(ddd(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(ddd(1,g(k,1):g(k,2),j)),'all'));
        m_arr(:,j,k) = mean(arr(:,g(k,1):g(k,2),j),2,'omitnan'); e_arr(:,j,k) = std(arr(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(ddd(1,g(k,1):g(k,2),j)),'all'));
        m_abb(:,j,k) = mean(abb(:,g(k,1):g(k,2),j),2,'omitnan'); e_abb(:,j,k) = std(abb(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(ddd(1,g(k,1):g(k,2),j)),'all'));
        m_euu(:,j,k) = mean(euu(:,g(k,1):g(k,2),j),2,'omitnan'); e_euu(:,j,k) = std(euu(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(ddd(1,g(k,1):g(k,2),j)),'all'));
        m_mmm(:,j,k) = mean(mmm(:,g(k,1):g(k,2),j),2,'omitnan'); e_mmm(:,j,k) = std(mmm(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(ddd(1,g(k,1):g(k,2),j)),'all'));
        m_rrr(:,j,k) = mean(rrr(:,g(k,1):g(k,2),j),2,'omitnan'); e_rrr(:,j,k) = std(rrr(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(ddd(1,g(k,1):g(k,2),j)),'all'));
    
        m_nddd(:,j,k) = mean(nddd(:,g(k,1):g(k,2),j),2,'omitnan'); e_nddd(:,j,k) = std(nddd(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(nddd(1,g(k,1):g(k,2),j)),'all'));
        m_narr(:,j,k) = mean(narr(:,g(k,1):g(k,2),j),2,'omitnan'); e_narr(:,j,k) = std(narr(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(nddd(1,g(k,1):g(k,2),j)),'all'));
        m_nabb(:,j,k) = mean(nabb(:,g(k,1):g(k,2),j),2,'omitnan'); e_nabb(:,j,k) = std(nabb(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(nddd(1,g(k,1):g(k,2),j)),'all'));
        m_neuu(:,j,k) = mean(neuu(:,g(k,1):g(k,2),j),2,'omitnan'); e_neuu(:,j,k) = std(neuu(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(nddd(1,g(k,1):g(k,2),j)),'all'));
        m_nmmm(:,j,k) = mean(nmmm(:,g(k,1):g(k,2),j),2,'omitnan'); e_nmmm(:,j,k) = std(nmmm(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(nddd(1,g(k,1):g(k,2),j)),'all'));
        m_nrrr(:,j,k) = mean(nrrr(:,g(k,1):g(k,2),j),2,'omitnan'); e_nrrr(:,j,k) = std(nrrr(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(nddd(1,g(k,1):g(k,2),j)),'all'));
    
        m_per(k,:,j) = mean(periods(g(k,1):g(k,2),:,j),1,'omitnan'); e_per(k,:,j) = std(periods(g(k,1):g(k,2),:,j),[],1,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(periods(g(k,1):g(k,2),1,j)),'all'));
        m_prat(:,j,k) = mean(p_rat(:,g(k,1):g(k,2),j),2,'omitnan'); e_prat(:,j,k) = std(p_rat(:,g(k,1):g(k,2),j),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(periods(g(k,1):g(k,2),1,j)),'all'));
    
    end

    m_ti(:,k) = mean(t_resp(g(k,1):g(k,2),:,1),1,'omitnan')'; e_ti(:,k) = std(t_resp(g(k,1):g(k,2),:,1),[],1,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(t_resp(g(k,1):g(k,2),1,1)),'all'));
    m_te(:,k) = mean(t_resp(g(k,1):g(k,2),:,2),1,'omitnan')'; e_te(:,k) = std(t_resp(g(k,1):g(k,2),:,2),[],1,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(t_resp(g(k,1):g(k,2),1,2)),'all'));
    
    m_rsp(:,k) = mean(rsp(:,g(k,1):g(k,2)),2,'omitnan'); e_rsp(:,k) = std(rsp(:,g(k,1):g(k,2)),[],2,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(rsp(1,g(k,1):g(k,2))),'all'));
    m_sumper(k,:) = mean(sum_per(g(k,1):g(k,2),:),1,'omitnan'); e_sumper(k,:) = std(sum_per(g(k,1):g(k,2),:),[],1,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(sum_per(g(k,1):g(k,2),1)),'all'));
    m_diffper(k,:) = mean(diff_per(g(k,1):g(k,2),:),1,'omitnan'); e_diffper(k,:) = std(diff_per(g(k,1):g(k,2),:),[],1,'omitnan')./sqrt(g(k,2)-g(k,1)+1-sum(isnan(diff_per(g(k,1):g(k,2),1)),'all'));
end

%% Plotting
time = linspace(2,20,10)';
time2 = linspace(0,20,11)';
colors = ['m';'r';'g';'y';'b';'k'];
ps = 0.15;

        dkwB = NaN(10,3);
        rkwB = NaN(10,3);
        bkwB = NaN(10,3);
        ekwB = NaN(10,3);
        mkwB = NaN(10,3);
        rlkwB = NaN(10,3);
for j = 1:3
    for m = 1:10 %each timepoint
        d_testmat = []; 
        gmat = [];
        gmat2 = [];
        r_testmat = [];
        b_testmat = [];
        e_testmat = [];
        m_testmat = [];
        % rl_testmat = [];

        %kw test
        for k = 1:nG-1 %wihtout control
            lenG = g(k,2)-g(k,1)+1; 
            % d_gmat(:,j) = strcat(d_gmat(:,j),repmat(strcat("g",num2str(k)),lenG,1));
            % gmat(end+1:end+lenG,1) = repmat(k,lenG,1);
            % 
        end
        % gmat = [1 1 1 1 1 1 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5]';
        gmat = [1 1 1 1 1 1 2 2 2 2 3 3 3 3 3 4 4 4 4 4 5 5 5 5 5 5 6 6 6 6 6]';
        for k = 1:nG %with control
            lenG = g(k,2)-g(k,1)+1; 
            % d_gmat(:,j) = strcat(d_gmat(:,j),repmat(strcat("g",num2str(k)),lenG,1));
            % gmat2(end+1:end+lenG,1) = repmat(k,lenG,1);
        end
        
        % inds = [1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 28];
        inds = [1 2 3 4 5 6 7 8 9 10 12 13 14 15 16 17 18 19 20 21 22 23 24 25 27 28 29 30 31 32 33];
        d_testmat(1:length(inds),1) = permute(ddd(m,inds,j),[2 1 3]); 
        [dKWp(m,j),dKWstats.(strcat('dia_',num2str(m),'_',num2str(j))),~] = kruskalwallis(d_testmat,gmat,'off');
        dDUNN{m,j} = dunn(d_testmat(:,1)',gmat');

        r_testmat(1:length(inds),1) = permute(arr(m,inds,j),[2 1 3]); 
        [rKWp(m,j),rKWstats.(strcat('air_',num2str(m),'_',num2str(j))),~] = kruskalwallis(r_testmat,gmat,'off');
        rDUNN{m,j} = dunn(r_testmat(:,1)',gmat');

        b_testmat(1:length(inds),1) = permute(abb(m,inds,j),[2 1 3]); 
        [bKWp(m,j),bKWstats.(strcat('abd_',num2str(m),'_',num2str(j))),~] = kruskalwallis(b_testmat,gmat,'off');
        bDUNN{m,j} = dunn(b_testmat(:,1)',gmat');

        e_testmat(1:length(inds),1) = permute(neuu(m,inds,j),[2 1 3]); 
        [eKWp(m,j),eKWstats.(strcat('euc_',num2str(m),'_',num2str(j))),~] = kruskalwallis(e_testmat,gmat,'off');
        eDUNN{m,j} = dunn(e_testmat(:,1)',gmat');

        m_testmat(1:length(inds),1) = permute(nmmm(m,inds,j),[2 1 3]); 
        [mKWp(m,j),mKWstats.(strcat('mah_',num2str(m),'_',num2str(j))),~] = kruskalwallis(m_testmat,gmat,'off');
        mDUNN{m,j} = dunn(m_testmat(:,1)',gmat');

        % rl_testmat(1:length(inds),1) = permute(rrr(m,inds,j),[2 1 3]);
        % [rlKWp(m,j),~,rlKWstats.(strcat('rrl_',num2str(m),'_',num2str(j)))] = kruskalwallis(rl_testmat,gmat,'off');
        % rlDUNN{m,j} = dunn(rl_testmat(:,1)',gmat');
        % c = multcompare(stats);
    end
    dkwB(dKWp(:,j)<0.05,j) = 1;
    rkwB(rKWp(:,j)<0.05,j) = 1;
    bkwB(bKWp(:,j)<0.05,j) = 1;
    ekwB(eKWp(:,j)<0.05,j) = 1;
    mkwB(mKWp(:,j)<0.05,j) = 1;
end

figure
for k = 1:nG
subplot(3,3,1); shadedErrorBar(time,m_arr(:,1,k),e_arr(:,1,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Airflow Area (norm.)'); hold on;
% ylim([-0.5 0.8]); 
ylim([-150 200]); 
title('Late Expiratory Period');
subplot(3,3,2); shadedErrorBar(time,m_arr(:,2,k),e_arr(:,2,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Airflow Area (norm.)'); hold on;
% ylim([-0.5 0.8]); 
ylim([-150 200]); 
title('Inspiratory Period');
subplot(3,3,3); shadedErrorBar(time,m_arr(:,3,k),e_arr(:,3,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Airflow Area (norm.)'); hold on;
% ylim([-0.5 0.8]); 
ylim([-150 200]); 
title('Post-Inspiratory Period');
        
subplot(3,3,4); shadedErrorBar(time,m_ddd(:,1,k),e_ddd(:,1,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Diaphragm Area (norm.)'); hold on;
% ylim([-0.4 0.8]); 
ylim([-100 150]); 
subplot(3,3,5); shadedErrorBar(time,m_ddd(:,2,k),e_ddd(:,2,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Diaphragm Area (norm.)'); hold on;
% ylim([-0.4 0.8]);
ylim([-100 150]); 
subplot(3,3,6); shadedErrorBar(time,m_ddd(:,3,k),e_ddd(:,3,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Diaphragm Area (norm.)'); hold on;
% ylim([-0.4 0.8]);
ylim([-100 150]); 

subplot(3,3,7); shadedErrorBar(time,m_abb(:,1,k),e_abb(:,1,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Abdominal Area (norm.)');hold on;
% ylim([-0.2 3.6]);
ylim([-50 1500]);
subplot(3,3,8); shadedErrorBar(time,m_abb(:,2,k),e_abb(:,2,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Abdominal Area (norm.)');hold on;
% ylim([-0.2 3.6]);
ylim([-50 1500]);
subplot(3,3,9); shadedErrorBar(time,m_abb(:,3,k),e_abb(:,3,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Abdominal Area (norm.)');hold on;
% ylim([-0.2 3.6]);
ylim([-50 1500]);
end
subplot(3,3,1); scatter(time,rkwB(:,1).*175,'*k'); 
subplot(3,3,2); scatter(time,rkwB(:,2).*175,'*k'); 
subplot(3,3,3); scatter(time,rkwB(:,3).*175,'*k');

subplot(3,3,4); scatter(time,dkwB(:,1).*130,'*k');
subplot(3,3,5); scatter(time,dkwB(:,2).*130,'*k');
subplot(3,3,6); scatter(time,dkwB(:,3).*130,'*k');

subplot(3,3,7); scatter(time,bkwB(:,1).*1300,'*k');
subplot(3,3,8); scatter(time,bkwB(:,2).*1300,'*k'); 
subplot(3,3,9); scatter(time,bkwB(:,3).*1300,'*k'); 

figure
for k = 1:nG
    subplot(2,3,1); shadedErrorBar(time,m_neuu(:,1,k),e_neuu(:,1,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Mean Euclidean Distance (norm.)'); hold on;
    ylim([0 3.6]); title('Late Expiratory Period');
    subplot(2,3,2); shadedErrorBar(time,m_neuu(:,2,k),e_neuu(:,2,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Mean Euclidean Distance (norm.)'); hold on;
    ylim([0 3.6]); title('Inspiratory Period');
    subplot(2,3,3); shadedErrorBar(time,m_neuu(:,3,k),e_neuu(:,3,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Mean Euclidean Distance (norm.)'); hold on;
    ylim([0 3.6]); title('Post-Inspiratory Period');

    subplot(2,3,4); shadedErrorBar(time,m_nmmm(:,1,k),e_nmmm(:,1,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Mean Mahalanobis Distance'); hold on;
    ylim([0 4100]);
    subplot(2,3,5); shadedErrorBar(time,m_nmmm(:,2,k),e_nmmm(:,2,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Mean Mahalanobis Distance'); hold on;
    ylim([0 4100]);
    subplot(2,3,6); shadedErrorBar(time,m_nmmm(:,3,k),e_nmmm(:,3,k),'lineprops',colors(k,:),'patchSaturation',ps); ylabel('Mean Mahalanobis Distance'); hold on;
    ylim([0 4100]);
end
subplot(2,3,1); scatter(time,ekwB(:,1).*3.5,'*k');
subplot(2,3,2); scatter(time,ekwB(:,2).*3.5,'*k'); 
subplot(2,3,3); scatter(time,ekwB(:,3).*3.5,'*k'); 

subplot(2,3,4); scatter(time,mkwB(:,1).*4000,'*k');
subplot(2,3,5); scatter(time,mkwB(:,2).*4000,'*k'); 
subplot(2,3,6); scatter(time,mkwB(:,3).*4000,'*k');

%Box plots
for j = 1:3
    for i = 1:nRec
        maxdiffind(i,j) = find(nmmm(:,i,j)==max(nmmm(:,i,j)));
        box(i,j) = nmmm(maxdiffind(i,j),i,j);
    end
    mx_testmat(i,1) = box(i,j);
    [MXp(j),~,MXstats.(strcat('mx_',num2str(j)))] = kruskalwallis(mx_testmat,gmat2,'off');
    mxDUNN{j} = dunn(mx_testmat(:,1)',gmat2');
end

newgmap = [6,1;2,2;5,3;3,4;1,5;4,6];
% newgmap = [1,5;2,2;3,4;4,6;5,3;6,1];%Old group designation to new group
glabel = {'Control','-0.25','0.35','0.5','0.6','0.9'};
gmat3 = []; newbox = []; box2 = [];

for j = 1:3
    gmat3 = [];
    for n = 1:nG
        
        grp = repmat(newgmap(n,1),(g(n,2)-g(n,1)+1),1);
        gmat3 = vertcat(gmat3,grp); 
        % gmat3(g(k,1):g(k,2),1) = ones(g(k,2)-g(k,1)+1,1).*newgmap(k,2);
        newbox = vertcat(newbox,box(g(newgmap(n,1),1):g(newgmap(n,1),2),j));
        % box2(gmat2==k,j) = box(g(k,1):g(k,2),j); 
    end
    box2(:,j) = newbox; newbox = []; 
end
figure
subplot(1,3,1); boxplot(box2(:,1),gmat2); ylabel('Mean Mahalanobis Distance'); set(gca,'Box','off');
set(gca,'XTickLabel',glabel); set(gca,'Box','off'); ylim([0 7000]); 
subplot(1,3,2); boxplot(box2(:,2),gmat2); ylabel('Mean Mahalanobis Distance'); set(gca,'Box','off');
set(gca,'XTickLabel',glabel); set(gca,'Box','off'); ylim([0 1000]); 
subplot(1,3,3); boxplot(box2(:,3),gmat2); ylabel('Mean Mahalanobis Distance'); set(gca,'Box','off');
set(gca,'XTickLabel',glabel); set(gca,'Box','off'); ylim([0 1500]); 

for j = 1:3
        [BXp(j),~,BXstats.(strcat('bx_',num2str(j)))] = kruskalwallis(box2(:,j),gmat3,'off');
        bxDUNN{j} = dunn(box2(:,j)',gmat3');
end
% for j ==2 (is normal)
[BXXp,~,BXXstats.(strcat('bx_',num2str(2)))] = anova1(box2(:,2)',gmat3');
% [c,m] = multcompare(BXXstats.(strcat('bx_',num2str(2))),'CriticalValueType','bonferroni'); 
[c,m] = multcompare(BXXstats.(strcat('bx_',num2str(2)))); 

for m = 1:11
    rs_testmat(:,1) = rsp(m,inds);
    [RSp(m),~,RSstats.(strcat('rs_',num2str(m)))] = kruskalwallis(rs_testmat,gmat,'off');
    rsDUNN{m} = dunn(rs_testmat(:,1)',gmat');

    ti_testmat(:,1) = t_resp(inds,m,1);
    [TIp(m),TIstats.(strcat('ti_',num2str(m))),~] = kruskalwallis(ti_testmat,gmat,'off');
    tiDUNN{m} = dunn(ti_testmat(:,1)',gmat');

    te_testmat(:,1) = t_resp(inds,m,2);
    [TEp(m),TEstats.(strcat('te_',num2str(m))),~] = kruskalwallis(te_testmat,gmat,'off');
    teDUNN{m} = dunn(te_testmat(:,1)',gmat');
end
tikwB = nan(11,1); tikwB(TIp<0.05) = 1;
tekwB = nan(11,1); tekwB(TEp<0.05) = 1;
rskwB = nan(11,1); rskwB(RSp<0.05) = 1;

figure
% time2 = linspace(0,20,11)';
for k = 1:nG
    subplot(1,3,1); shadedErrorBar(time2,m_rsp(:,k),e_rsp(:,k),'lineprops',colors(k,:),'patchSaturation',ps); hold on
    ylabel('Respiratory Period (s)');
    subplot(1,3,2); shadedErrorBar(time2,m_ti(:,k),e_ti(:,k),'lineprops',colors(k,:),'patchSaturation',ps); hold on
    ylabel('Inspiratory Period (s)');
    subplot(1,3,3); shadedErrorBar(time2,m_te(:,k),e_te(:,k),'lineprops',colors(k,:),'patchSaturation',ps); hold on
    ylabel('Expiratory Period (s)'); 
end

subplot(1,3,1); scatter(time2,rskwB.*1.33,'*k');
subplot(1,3,2); scatter(time2,tikwB.*1.02,'*k');
subplot(1,3,3); scatter(time2,tekwB.*1.55,'*k');