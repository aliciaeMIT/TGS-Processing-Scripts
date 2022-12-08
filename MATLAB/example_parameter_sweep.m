%------------------------------------------------------------------------------------------------------------------------------------------------------
%Script to sweep over some input values of TGSPhaseAnalysis for the purposes of parameter sensitivity analysis
%This will dump 2 files for each spot you have taken - one sweeping across the grating parameter and one across the initial peak skip (fit start point)
%In the files will be values and errors for each of the seven fitted parameters in the TGS equation as well as the SAW frequency.
%------------------------------------------------------------------------------------------------------------------------------------------------------
% directory location
pname = 'example_data\iron_parameter_study'; %ADJUST
dir=cd(pname);
cd(dir)
%First part of file name and number of spots
str_base='IronIrradiated_Fe3plus-2022-01-04-06.40um-unirradiated_';  %ADJUST
%Change filename modifying number here
spots=5:20; %ADJUST
%Baseline handling - you always have to supply some form of baseline anyway for this to work, even if it is not used!
baselineBool = 1; % 1 for "yes do baseline subtraction", 0 for "no don't do that" 
POSbaselineStr = pname + string('\IronIrradiated_Fe3plus-2022-01-04-06.40um-baseline-POS-1.txt'); %ADJUST
NEGbaselineStr = pname + string('\IronIrradiated_Fe3plus-2022-01-04-06.40um-baseline-NEG-1.txt'); %ADJUST

% Calibrated grating spacing
grat=6.4190; %in um %ADJUST

%END INPUTS (broadly)----------------------------------------------------------------------------------------------------------------------------------

% Initialize all outputs
freq=0;
freq_err=0;
speed=0;
diff=0;
diff_err=0;
tau=zeros(3,length(spots));
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
disp('Begin parameter sweeping');

%Run through each of the files, sweeping across gratings and peak-skips
for jj=1:length(spots)
    close all
    pos_str=strcat(pname,string('\'),str_base,num2str(spots(jj)),'-POS-1.txt');
    neg_str=strcat(pname,string('\'),str_base,num2str(spots(jj)),'-NEG-1.txt');
    
    %Set up the output files with their headers, clears the file at beginning
    %Grating sensitivity  
    printFileGratingSensitivity = pname + string('\') + str_base + num2str(spots(jj)) + string('_grating_sensitivity.txt');
    fidGratingSensitivity = fopen( printFileGratingSensitivity, 'wt' );
    fprintf(fidGratingSensitivity,'%s', 'grating_mod grating_value[um] SAW_freq[Hz] SAW_freq_error[Hz] A[Wm^-2] A_err[Wm^-2] alpha[m^2s^-1] alpha_err[m2s-1] beta[s^0.5] beta_err[s^0.5] B[Wm^-2] B_err[Wm^-2] theta theta_err tau[s] tau_err[s] C[Wm^-2] C_err[Wm^-2]');
    fclose(fidGratingSensitivity);
    
    %Loop for the grating spacing scan, default to recommended peak skip of 2.
    %Choose your range of scan here, e.g. i=-5:5 for a sweep between 95% and 105% of grating spacing.
    for i=-5:5 
        gratingMod = 1+i/100; %can edit this from 100 to 1000 if you want a finer sweep will take more time obvs.
        disp(['current run is: ', pos_str, '   with grating modifier: ', num2str(gratingMod)]);
        [freq,freq_err,speed,diff,diff_err,tau(:,jj),tauErr, A, AErr, beta, betaErr, B, BErr, theta, thetaErr, C, CErr]=TGSPhaseAnalysis(pos_str,neg_str,grat*gratingMod,2,0,1,POSbaselineStr,NEGbaselineStr);
        fidGratingSensitivity = fopen( printFileGratingSensitivity, 'a' );
        fprintf(fidGratingSensitivity, string('\n%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g'), gratingMod, gratingMod*grat, freq, freq_err, A, AErr, diff, diff_err, beta, betaErr, B, BErr, theta, thetaErr, tau(3,jj), tauErr, C, CErr);
        fclose(fidGratingSensitivity);
        close all
    end
    
    %Peak skip sensitivity
    printFilePeakSkipSensitivity = pname + string('\') + str_base + num2str(spots(jj)) + string('_peak_skip_sensitivity.txt');
    fidPeakSkipSensitivity = fopen( printFilePeakSkipSensitivity, 'wt' );
    fprintf(fidPeakSkipSensitivity,'%s', 'half_peak_skip_number SAW_freq[Hz] SAW_freq_error[Hz] A[Wm^-2] A_err[Wm^-2] alpha[m^2s^-1] alpha_err[m2s-1] beta[s^0.5] beta_err[s^0.5] B[Wm^-2] B_err[Wm^-2] theta theta_err tau[s] tau_err[s] C[Wm^-2] C_err[Wm^-2]');
    fclose(fidPeakSkipSensitivity);
    
    %Loop for the peak skipping scan.
    for i=1:4 %The TGS_phase_fft script is hard-coded to freak out above 4
        disp(['current run is: ', pos_str, '   with peak-skip parameter: ', num2str(i)]);
        [freq,freq_err,speed,diff,diff_err,tau(:,jj),tauErr, A, AErr, beta, betaErr, B, BErr, theta, thetaErr, C, CErr]=TGSPhaseAnalysis(pos_str,neg_str,grat,i,0, 1, POSbaselineStr, NEGbaselineStr);
        fidPeakSkipSensitivity = fopen( printFilePeakSkipSensitivity, 'a' );
        fprintf(fidPeakSkipSensitivity, string('\n%.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g'), i, freq, freq_err, A, AErr, diff, diff_err, beta, betaErr, B, BErr, theta, thetaErr, tau(3,jj), tauErr, C, CErr);
        fclose(fidPeakSkipSensitivity);
        close all
    end
end

disp('F L A W L E S S');