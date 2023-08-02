% directory location
pname = 'C:\Users\apcwy\OneDrive - Massachusetts Institute of Technology\MIT_postdoc\TGS\Data\PI3\Y0.7Cu0.7Mo\sample_1.1'; %ADJUST
%First part of file name
str_base='Y0.7Cu0.7Mo-sample_1.1-06.40um-spot1-';  %ADJUST
% Calibrated grating spacings
grat = 6.6213; %in um %ADJUST

% baseline handling
blbool = 1;
baselinePOS = pname + string('\') + 'Y0.7Cu0.7Mo-sample_1.1-06.40um-baseline-POS-1.txt';
baselineNEG = pname + string('\') + 'Y0.7Cu0.7Mo-sample_1.1-06.40um-baseline-NEG-1.txt';

% Initialize all outputs
freq=0;
freq_err=0;
speed=0;
diff=0;
diff_err=0;
tau=zeros(3,158); % The second dimension of this array needs to match the total number of times TGSPhaseAnalysis is run
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

%Run through each of the files, generating a new filename each time, then
%write all the outputs to a file
for i=1:158
    close all
    pos_str = strcat(pname,string('\'),str_base,"POS-",num2str(i),'.txt');
    neg_str = strcat(pname,string('\'),str_base,"NEG-",num2str(i),'.txt');

    disp(['current run is: ', pos_str, '  with grating: ', num2str(grat),' um']);
    disp(['current baseline is: ', pos_str]);
    [freq,freq_err,speed,diff,diff_err,tau(:,3*(i)),tauErr, A, AErr, beta, betaErr, B, BErr, theta, thetaErr, C, CErr]=TGSPhaseAnalysis(pos_str,neg_str,grat,2,0);

    fid1 = fopen( printFile, 'a' );
    naming_str=strcat(str_base,num2str(i));
    fprintf(fid1, string('\n%s %.8g %0.5e %0.5e %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g'), naming_str, grat, freq, freq_err, A, AErr, diff, diff_err, beta, betaErr, B, BErr, theta, thetaErr, tau(3,(i)), tauErr, C, CErr);
    fclose(fid1);
end
close all
disp('F L A W L E S S');
