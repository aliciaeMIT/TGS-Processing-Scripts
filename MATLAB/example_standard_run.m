% directory location
pname = 'example_data\eurofer_temperature_study'; %ADJUST
dir=cd(pname);
cd(dir)
%First part of file name and number of spots
str_base='Eurofer-HT_weld_2022-06-06-06.40um-';  %ADJUST
%Change filename modifying number here
spots=0:4; %ADJUST
%Baseline handling if no baseline arguments are provided, the code still works - it just doesn't do any baseline correction
baselineBool = 1; % 1 for "yes do baseline subtraction", 0 for "no don't do that" 
POSbaselineStr = pname + string('\') + 'Eurofer-HT_weld_2022-06-06-06.40um-baseline-POS-1.txt';
NEGbaselineStr = pname + string('\') + 'Eurofer-HT_weld_2022-06-06-06.40um-baseline-NEG-1.txt';
% Calibrated grating spacings
grat = 6.4611; %in um %ADJUST

% Initialize all outputs
freq=0;
freq_err=0;
speed=0;
diff=0;
diff_err=0;
tau=zeros(3,45); % The second dimension of this array needs to match the total number of times TGSPhaseAnalysis is run
tau_err = 0;
A = 0;
AErr = 0;
beta = 0;
betaErr = 0;
B = 0;
BErr = 0;
theta = 0;
thetaErr = 0;
C = 0;
CErr = 0;
formatSpec = '%.1f';
%Set up the output files with their headers, clears the file at beginning
printFile = pname + string('\') + str_base + string('_postprocessing.txt');
fid1 = fopen( printFile, 'wt' );
fprintf(fid1,'%s', 'run_name grating_value[um] SAW_freq[Hz] SAW_freq_error[Hz] A[Wm^-2] A_err[Wm^-2] alpha[m^2s^-1] alpha_err[m2s-1] beta[s^0.5] beta_err[s^0.5] B[Wm^-2] B_err[Wm^-2] theta theta_err tau[s] tau_err[s] C[Wm^-2] C_err[Wm^-2]');
fclose(fid1);
temprange = [24,50,100,150,200,250,300,350,400,450,500,550,600,650,700];


%Run through each of the files, sweeping across gratings and peak-skips
for i=1:length(temprange)
    for j=0:2
        close all
        pos_str=strcat(pname,string('\'),str_base,num2str(temprange(i)),'C_spot',num2str(j),'-POS-1.txt');
        neg_str=strcat(pname,string('\'),str_base,num2str(temprange(i)),'C_spot',num2str(j),'-NEG-1.txt');
        %POSbaselineStr = strcat(str_base,num2str(i),'-POS-1.txt'); %ADJUST
        %NEGbaselineStr = strcat(str_base,num2str(i),'-POS-1.txt'); %ADJUST
        
        disp(['current run is: ', pos_str, '   with grating: ', num2str(grat),' um']);
        disp(['current baseline is: ', POSbaselineStr]);
        [freq,freq_err,speed,diff,diff_err,tau(:,3*(i-1)+j+1),tauErr, A, AErr, beta, betaErr, B, BErr, theta, thetaErr, C, CErr]=TGSPhaseAnalysis(pos_str,neg_str,grat,2,0, baselineBool, POSbaselineStr, NEGbaselineStr);
        
        fid1 = fopen( printFile, 'a' );
        fprintf(fid1, string('\n%s %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g'), pos_str, grat, freq, freq_err, A, AErr, diff, diff_err, beta, betaErr, B, BErr, theta, thetaErr, tau(3,3*(i-1)+j+1), tauErr, C, CErr);
        fclose(fid1);
    end
end
close all
disp('F L A W L E S S');
