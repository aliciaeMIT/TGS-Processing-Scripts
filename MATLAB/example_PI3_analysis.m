
% directory location
pname = 'C:\Users\apcwy\OneDrive - Massachusetts Institute of Technology\MIT_postdoc\TGS\Data\Kyocera\Kyocera_W_W\Tungsten_data\2024-02-05'; %ADJUST

%First part of file name and number of spots
str_base='Tungsten_Kyocera-2024-02-05-03.40um-spot03-';  %ADJUST
%Change filename modifying number here
numberOfRuns=305; %ADJUST

% Calibrated grating spacings
grat = 3.5217; %in um %ADJUST

% Initialize all outputs
freq=0;
freq_err=0;
speed=0;
diff=0;
diff_err=0;
tau=zeros(3,numberOfRuns); 
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
file_date_time = ' ';
formatSpec = '%.1f';
%Set up the output files with their headers, clears the file at beginning
printFile = pname + string('\') + str_base + string('_postprocessing.txt');
fid1 = fopen( printFile, 'wt' );
fprintf(fid1,'%s', 'run_name grating_value[um] date_time SAW_freq[Hz] SAW_freq_error[Hz] A[Wm^-2] A_err[Wm^-2] alpha[m^2s^-1] alpha_err[m2s-1] beta[s^0.5] beta_err[s^0.5] B[Wm^-2] B_err[Wm^-2] theta theta_err tau[s] tau_err[s] C[Wm^-2] C_err[Wm^-2]');
fclose(fid1);

%Run through each of the files, sweeping across gratings and peak-skips
for j=1:numberOfRuns
    close all
    pos_str=strcat(pname,string('\'),str_base,'POS-',num2str(j),'.txt');
    neg_str=strcat(pname,string('\'),str_base,'NEG-',num2str(j),'.txt');
    disp(['current run is: ', pos_str, '   with grating: ', num2str(grat),' um']);
    [freq,freq_err,speed,diff,diff_err,tau(:,3*(j)),tauErr, A, AErr, beta, betaErr, B, BErr, theta, thetaErr, C, CErr, file_date_time]=TGSPhaseAnalysis(pos_str,neg_str,grat,2,0);
    
    %print these to the postprocessing file
    fid1 = fopen( printFile, 'a' );
    fprintf(fid1, string('\n%s %s %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g'), pos_str, file_date_time, grat, freq, freq_err, A, AErr, diff, diff_err, beta, betaErr, B, BErr, theta, thetaErr, tau(3,3*(j)), tauErr, C, CErr);
    fclose(fid1);
end

close all
disp('F L A W L E S S');
