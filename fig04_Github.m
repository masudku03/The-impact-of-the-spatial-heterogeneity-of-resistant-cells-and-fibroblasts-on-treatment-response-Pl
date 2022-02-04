clear; clc; tMax=2000; T =(1:tMax)'; TCTC = zeros(30,1); TCTR = zeros(30,1); TCTU = zeros(30,1); TATC = zeros(30,1); TATR = zeros(30,1); TATU = zeros(30,1);
%%%%from 210312 clamp

KC1 = table2array(readtable('C:\SimulationData\210312\CT_K_1_m0.0_rho1.020.0189_clumpInput_010_0.csv', 'HeaderLines',4)); TCTC(1) = find(KC1(1:tMax,4)>(1.2*5000),1);
for ii = 1:29
    filename = sprintf('C:\\SimulationData\\210312\\CT_K_1_m0.0_rho1.020.0189_clumpInput_010_%d.csv',ii)
    tempKC = table2array(readtable(filename, 'HeaderLines',4)); KC1 = KC1 + tempKC; TCTC(ii) = find(tempKC(1:tMax,4)>(1.2*5000),1);
end
KC1 = KC1./30;
figure(10)
subplot(6,6,4:6)
plot(T,KC1(1:tMax,4),'-k','LineWidth',2); hold on; area(T,KC1(1:tMax,4),'LineStyle','none','FaceColor','#0072BD'); hold on;% 'Color',[0, 0.4470, 0.7410]
area(T,KC1(1:tMax,3),'LineStyle','none','FaceColor','#D95319'); subtitle('Clumped'); ylabel('Population'); xlabel('Time(Days)');
plot([0 mean(TCTC)],[6000 6000],'c-','LineWidth',1); 
plot([mean(TCTC) mean(TCTC)],[0 6000],'c-','LineWidth',1); 
title('(b)'); ax = gca; ax.TitleHorizontalAlignment = 'left';

%%Random
KC1 = table2array(readtable('C:\SimulationData\210312\CT_K_1_m0.0_rho1.020.0189_randomInput_010_0.csv', 'HeaderLines',4));  TCTR(1) = find(KC1(1:tMax,4)>(1.2*5000),1);
for ii = 1:29
    filename = sprintf('C:\\SimulationData\\210312\\CT_K_1_m0.0_rho1.020.0189_randomInput_010_%d.csv',ii)
    tempKC = table2array(readtable(filename, 'HeaderLines',4)); KC1 = KC1 + tempKC; TCTR(ii) = find(tempKC(1:tMax,4)>(1.2*5000),1);
end

KC1 = KC1./30;
figure(10)
subplot(6,6,10:12)
plot(1:tMax,KC1(1:tMax,4),'-k','LineWidth',2); hold on; area(1:tMax,KC1(1:tMax,4),'LineStyle','none','FaceColor','#0072BD'); hold on;%'Color',[0.8500, 0.3250, 0.0980] 
area(1:tMax,KC1(1:tMax,3),'LineStyle','none','FaceColor','#D95319'); subtitle('Random'); ylabel('Population'); xlabel('Time(Days)');
plot([0 mean(TCTR)],[6000 6000],'c-','LineWidth',1); 
plot([mean(TCTR) mean(TCTR)],[0 6000],'c-','LineWidth',1); 

%%Uniform
KC1 = table2array(readtable('C:\SimulationData\210312\CT_K_1_m0.0_rho1.020.0189_unifInput_0.csv', 'HeaderLines',4)); TCTU(1) = find(KC1(1:tMax,4)>(1.2*5000),1);
for ii = 1:29
    filename = sprintf('C:\\SimulationData\\210312\\CT_K_1_m0.0_rho1.020.0189_unifInput_%d.csv',ii)
    tempKC = table2array(readtable(filename, 'HeaderLines',4)); KC1 = KC1 + tempKC; TCTR(ii) = find(tempKC(1:tMax,4)>(1.2*5000),1);
end
KC1 = KC1./30;
figure(10)
subplot(6,6,16:18)
plot(1:tMax,KC1(1:tMax,4),'-k','LineWidth',2); hold on;area(1:tMax,KC1(1:tMax,4),'LineStyle','none','FaceColor','#0072BD'); hold on;%'Color',[0.9290, 0.6940, 0.1250] 
area(1:tMax,KC1(1:tMax,3),'LineStyle','none','FaceColor','#D95319'); legend('$\overline{N}(t)$','$\overline{S}(t)$','$\overline{R}(t)$','Interpreter','latex');%,'Color',[0.9290, 0.6940, 0.1250]
plot([0 mean(TCTU)],[6000 6000],'c-','LineWidth',1); 
plot([mean(TCTU) mean(TCTU)],[0 6000],'c-','LineWidth',1); 
subtitle('Uniform'); ylabel('Population'); xlabel('Time(Days)'); 

%% neighborhood
clear;clc;  %%% 210312
figure(10)
subplot(6,6,1)
mymap = [1 1 1; 0 0.4470 0.7410; 0.8500 0.3250 0.0980]; colormap(mymap);
imagesc(csvread('C:\SimulationData\210312\CT_K1_m000_Clump_f010_r29\CT_K_1_m0.0_rho1.020.0189_clumpInput_010_time_1_29.csv'),[0.2,1.5]);xticks([]); yticks([]); ylabel('Clumped'); title('(a)'); ax = gca; ax.TitleHorizontalAlignment = 'left';
subplot(6,6,2)
imagesc(csvread('C:\SimulationData\210312\CT_K1_m000_Clump_f010_r29\CT_K_1_m0.0_rho1.020.0189_clumpInput_010_time_120_29.csv'),[0.2,1.5]);xticks([]); yticks([]); 
subplot(6,6,3)
imagesc(csvread('C:\SimulationData\210312\CT_K1_m000_Clump_f010_r29\CT_K_1_m0.0_rho1.020.0189_clumpInput_010_time_1999_29.csv'),[0.2,1.5]);xticks([]); yticks([]);

meanRc = zeros(30,3); meanSc = zeros(30,3); meanEc = zeros(30,3);

for ii = 0:29
    filename = sprintf('C:\\SimulationData\\210312\\CT_K1_m000_Clump_f010_r%d\\CT_K_1_m0.0_rho1.020.0189_clumpInput_010_time_1_%d.csv',ii,ii)
    DC{ii+1} = csvread(filename);
end
j = 1;
for i = 1:30
[meanRc(i,j), meanSc(i,j), meanEc(i,j)] = nbdHood(DC{i}); 
end

for ii = 0:29
    filename = sprintf('C:\\SimulationData\\210312\\CT_K1_m000_Clump_f010_r%d\\CT_K_1_m0.0_rho1.020.0189_clumpInput_010_time_60_%d.csv',ii,ii)
    DC{ii+1} = csvread(filename);
end
j = 2;
for i = 1:30
[meanRc(i,j), meanSc(i,j), meanEc(i,j)] = nbdHood(DC{i}); 
end

for ii = 0:29
    filename = sprintf('C:\\SimulationData\\210312\\CT_K1_m000_Clump_f010_r%d\\CT_K_1_m0.0_rho1.020.0189_clumpInput_010_time_300_%d.csv',ii,ii)
    DC{ii+1} = csvread(filename);
end
j = 3;
for i = 1:30
[meanRc(i,j), meanSc(i,j), meanEc(i,j)] = nbdHood(DC{i}); 
end


%% %% Random 210312
subplot(6,6,7)
imagesc(csvread('C:\SimulationData\210312\CT_K1_m000_Random_f010_r29\CT_K_1_m0.0_rho1.020.0189_randomInput_010_time_1_29.csv'),[0.2,1.5]);xticks([]); yticks([]); ylabel('Random');
subplot(6,6,8)
imagesc(csvread('C:\SimulationData\210312\CT_K1_m000_Random_f010_r29\CT_K_1_m0.0_rho1.020.0189_randomInput_010_time_120_29.csv'),[0.2,1.5]);xticks([]); yticks([]);
subplot(6,6,9)
imagesc(csvread('C:\SimulationData\210312\CT_K1_m000_Random_f010_r29\CT_K_1_m0.0_rho1.020.0189_randomInput_010_time_1999_29.csv'),[0.2,1.5]);xticks([]); yticks([]);

meanRrR = zeros(30,3); meanSrR = zeros(30,3); meanErR = zeros(30,3);

for ii = 0:29
    filename = sprintf('C:\\SimulationData\\210312\\CT_K1_m000_Random_f010_r%d\\CT_K_1_m0.0_rho1.020.0189_randomInput_010_time_1_%d.csv',ii,ii)
    DC{ii+1} = csvread(filename);
end
j = 1;
for i = 1:30
[meanRrR(i,j), meanSrR(i,j), meanErR(i,j)] = nbdHood(DC{i}); 
end

for ii = 0:29
    filename = sprintf('C:\\SimulationData\\210312\\CT_K1_m000_Random_f010_r%d\\CT_K_1_m0.0_rho1.020.0189_randomInput_010_time_60_%d.csv',ii,ii)
    DC{ii+1} = csvread(filename);
end
j = 2;
for i = 1:30
[meanRrR(i,j), meanSrR(i,j), meanErR(i,j)] = nbdHood(DC{i}); 
end

for ii = 0:29
    filename = sprintf('C:\\SimulationData\\210312\\CT_K1_m000_Random_f010_r%d\\CT_K_1_m0.0_rho1.020.0189_randomInput_010_time_300_%d.csv',ii,ii)
    DC{ii+1} = csvread(filename);
end
j = 3;
for i = 1:30
[meanRrR(i,j), meanSrR(i,j), meanErR(i,j)] = nbdHood(DC{i}); 
end


%% %% Uniform 210312
subplot(6,6,13)
imagesc(csvread('C:\SimulationData\210312\CT_K_1_m0.0_rho1.020.0189_unifInput_time_1_29.csv'),[0.2,1.5]);xticks([]); yticks([]); ylabel('Uniform'); xlabel('Day 1');
subplot(6,6,14)
imagesc(csvread('C:\SimulationData\210312\CT_K_1_m0.0_rho1.020.0189_unifInput_time_120_29.csv'),[0.2,1.5]);xticks([]); yticks([]); xlabel('Day 120');
subplot(6,6,15)
imagesc(csvread('C:\SimulationData\210312\CT_K_1_m0.0_rho1.020.0189_unifInput_time_1999_29.csv'),[0.2,1.5]);xticks([]); yticks([]);xlabel(' Day 2000');


meanRuU = zeros(30,3); meanSuU = zeros(30,3);  meanEuU = zeros(30,3);

for ii = 0:29
    filename = sprintf('C:\\SimulationData\\210312\\CT_K1_m000_Random_f010_r%d\\CT_K_1_m0.0_rho1.020.0189_randomInput_010_time_1_%d.csv',ii,ii)
    DC{ii+1} = csvread(filename);
end
j = 1;
for i = 1:30
[meanRuU(i,j), meanSuU(i,j), meanEuU(i,j)] = nbdHood(DC{i}); 
end

for ii = 0:29
    filename = sprintf('C:\\SimulationData\\210312\\CT_K1_m000_Random_f010_r%d\\CT_K_1_m0.0_rho1.020.0189_randomInput_010_time_60_%d.csv',ii,ii)
    DC{ii+1} = csvread(filename);
end
j = 2;
for i = 1:30
[meanRuU(i,j), meanSuU(i,j), meanEuU(i,j)] = nbdHood(DC{i}); 
end

for ii = 0:29
    filename = sprintf('C:\\SimulationData\\210312\\CT_K1_m000_Random_f010_r%d\\CT_K_1_m0.0_rho1.020.0189_randomInput_010_time_300_%d.csv',ii,ii)
    DC{ii+1} = csvread(filename);
end
j = 3;
for i = 1:30
[meanRuU(i,j), meanSuU(i,j), meanEuU(i,j)] = nbdHood(DC{i}); 
end

subplot(6,6,[22:24,28:30,34:36])
boxchart([ones(30,1);5*ones(30,1);9*ones(30,1)],[meanRc(:,1);meanRc(:,2);meanRc(:,3)]); hold on; 
boxchart([1.5*ones(30,1);5.5*ones(30,1);9.5*ones(30,1)],[meanRrR(:,1);meanRrR(:,2);meanRrR(:,3)]); hold on;
boxchart([2*ones(30,1);6*ones(30,1);10*ones(30,1)],[meanRuU(:,1);meanRuU(:,2);meanRuU(:,3)]); box on; 
title('(d)'); ax = gca; ax.TitleHorizontalAlignment = 'left'; hold on; 

xticks([1.5 5.5 9.5]); xticklabels({'1','60','300'}); ylabel('$\mathcal{N}_{R\overline{R}}^i$','Interpreter','latex'); xlabel('Time(Days)'); xlim([0 10.6]);

subplot(6,6,[19:21,25:27,31:33])
hold on;
boxchart([ones(30,1);5*ones(30,1);9*ones(30,1)],[meanEc(:,1);meanEc(:,2);meanEc(:,3)]); hold on;
boxchart([1.5*ones(30,1);5.5*ones(30,1);9.5*ones(30,1)],[meanErR(:,1);meanErR(:,2);meanErR(:,3)]); hold on;
boxchart([2*ones(30,1);6*ones(30,1);10*ones(30,1)],[meanEuU(:,1);meanEuU(:,2);meanEuU(:,3)]);
xticks([1.5 5.5 9.5]); xticklabels({'1','60','300'}); ylabel('$\mathcal{N}_{E\overline{R}}^i$','Interpreter','latex'); xlabel('Time(Days)'); xlim([0 10.6]);
box on; title('(c)'); ax = gca; ax.TitleHorizontalAlignment = 'left'; hold on; legend('i=c','i=r','i=u'); 

function [meanR, meanS, meanE] = nbdHood(DC)    
       [rR,cR] = find(DC==2);
    countR = zeros(length(rR),1); countS = zeros(length(rR),1); countE = zeros(length(rR),1);
    for ii=1:1:length(rR)
        if (rR(ii)==1 && cR(ii)~=1 && cR(ii)~=100)
            nbd = [DC(rR(ii),cR(ii)-1), DC(rR(ii),cR(ii)+1), DC(rR(ii)+1,cR(ii))];
        elseif (rR(ii)==100 && cR(ii)~=1 && cR(ii)~=100)
            nbd = [DC(rR(ii),cR(ii)-1), DC(rR(ii),cR(ii)+1), DC(rR(ii)-1,cR(ii))];
        elseif (rR(ii)~=1 && rR(ii)~=100 && cR(ii)==1)
            nbd = [DC(rR(ii),cR(ii)+1), DC(rR(ii)+1,cR(ii)), DC(rR(ii)-1,cR(ii))];
        elseif (rR(ii)~=1 && rR(ii)~=100 && cR(ii)==100)
            nbd = [DC(rR(ii),cR(ii)-1), DC(rR(ii)-1,cR(ii)), DC(rR(ii)+1,cR(ii))];
        elseif (rR(ii)==100 && cR(ii)==100)
            nbd = [DC(rR(ii),cR(ii)-1), DC(rR(ii)-1,cR(ii))];
        elseif (rR(ii)==1 && cR(ii)==1)
            nbd = [DC(rR(ii),cR(ii)+1), DC(rR(ii)+1,cR(ii))];
        elseif (rR(ii)==100 && cR(ii)==1)
            nbd = [DC(rR(ii),cR(ii)+1), DC(rR(ii)-1,cR(ii))];
        elseif (rR(ii)==1 && cR(ii)==100)
            nbd = [DC(rR(ii)+1,cR(ii)), DC(rR(ii),cR(ii)-1)];
        else
            nbd = [DC(rR(ii),cR(ii)-1), DC(rR(ii),cR(ii)+1), DC(rR(ii)-1,cR(ii)), DC(rR(ii)+1,cR(ii))];    
        end

        countR(ii) = length(nbd(nbd==2)); countS(ii) = length(nbd(nbd==1)); countE(ii) = length(nbd(nbd==0));
    end
    meanR = mean(countR); meanS = mean(countS); meanE = mean(countE);
end