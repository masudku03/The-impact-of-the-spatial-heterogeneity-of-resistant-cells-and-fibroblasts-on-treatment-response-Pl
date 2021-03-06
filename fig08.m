clear;
figure(62)
day = [60 165 250 550];
    mymap = [1 1 1; 0.4940 0.1840 0.5560; 0 0.4470 0.7410; 0.8500 0.3250 0.0980]; colormap(mymap);
% Fib_SQC.csv Fib_SQO.csv Fib_SQ2.csv Fib_OU1.csv Fib_clump.csv are the csv
% files with infromation regarding the fibroblas structure.
% the other csv file should be generated by execeuting the java file    
for ii = 1:1:numel(day)
    subplot(6,4,ii)
    filename = sprintf('C:\\SimulationData\\210312\\CT_K1_m000_Clump_f010_r29\\CT_K_1_m0.0_rho1.020.0189_clumpInput_010_time_%d_29.csv',day(ii))
    im = imagesc(csvread(filename));xticks([]); yticks([]); im.AlphaData = 0.8;
	ylabel('NoF');

    subplot(6,4,4+ii)
    filename = sprintf('C:\\SimulationData\\K1CTm000FBSQCFI200\\CT_K_1_m0.0_rho1.020.0189SQC200.0_clumpInput_010_time_%d_29.csv',day(ii))
    imagesc(0.5*csvread(sprintf('C:\\SimulationData\\Fib_SQC.csv')));xticks([]); yticks([]); hold on; imC = imagesc(csvread(filename)); imC.AlphaData = .8; 
    ylabel('FC');
    subplot(6,4,8+ii)
    filename = sprintf('C:\\SimulationData\\K1CTm000FBSQOFI200\\CT_K_1_m0.0_rho1.020.0189SQO200.0_clumpInput_010_time_%d_29.csv',day(ii))
    imagesc(0.5*csvread(sprintf('C:\\SimulationData\\Fib_SQO.csv')));xticks([]); yticks([]); hold on; imC = imagesc(csvread(filename)); imC.AlphaData = .8; 
    ylabel('FCp');
    subplot(6,4,12+ii)
    filename = sprintf('C:\\SimulationData\\K1CTm000FBSQ2FI200\\CT_K_1_m0.0_rho1.020.0189SQ2200.0_clumpInput_010_time_%d_29.csv',day(ii))
    imagesc(0.5*csvread(sprintf('C:\\SimulationData\\Fib_SQ2.csv')));xticks([]); yticks([]); hold on; im2 = imagesc(csvread(filename)); im2.AlphaData = .8; 
    ylabel('FSq');
    subplot(6,4,16+ii)
    filename = sprintf('C:\\SimulationData\\K1CTm000FBOU1FI200\\CT_K_1_m0.0_rho1.020.0189OU1200.0_clumpInput_010_time_%d_29.csv',day(ii))
    imagesc(0.5*csvread(sprintf('C:\\SimulationData\\Fib_OU1.csv')));xticks([]); yticks([]); hold on; im2 = imagesc(csvread(filename)); im2.AlphaData = .8; 
    ylabel('FSp');
    subplot(6,4,20+ii)
    filename = sprintf('C:\\SimulationData\\K1CTm000FBclumpsFI200_1000\\CT_K_1_m0.0_rho1.020.0189_RDC_200.0clumpInput_010_time_%d_29.csv',day(ii))
    imagesc(0.5*csvread(sprintf('C:\\SimulationData\\Fib_clump.csv')));xticks([]); yticks([]); hold on; im2 = imagesc(csvread(filename)); im2.AlphaData = .8; xlabel(['Day ',num2str(day(ii))]);
	ylabel('FR');
end
