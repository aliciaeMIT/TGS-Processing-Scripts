% directory location
pname = 'C:\Users\apcwy\OneDrive - Massachusetts Institute of Technology\MIT_postdoc\TGS\Data\PI3\UMoNbZr\2023-07-25'; %ADJUST
%First part of file name
str_base='UMoNbZr-2023-07-25-04.20um-time';  %ADJUST
% Calibrated grating spacings
grat = 4.3323; %in um %ADJUST

% baseline handling
blbool = 1;
baselinePOS = pname + string('\') + 'UMoNbZr-2023-07-25-04.20um-baseline-POS-1.txt';
baselineNEG = pname + string('\') + 'UMoNbZr-2023-07-25-04.20um-baseline-NEG-1.txt';

printFile = pname + string('\') + str_base + string('_postprocessing.txt');
fid1 = fopen( printFile, 'a' );

posstr = strcat(pname,string('\'),str_base,num2str(i),'-POS-1.txt');
negstr = strcat(pname,string('\'),str_base,num2str(i),'-POS-1.txt');

[freq_final,freq_error,speed,diffusivity,diffusivity_err,tau, tauErr, paramrA, AErr, ParamBeta, BetaErr, paramB, BErr, paramTheta, thetaErr, paramC, CErr] = TGSPhaseAnalysis(posstr,negstr,grat,2,0);
disp('Frequency = ');
disp(freq_final);
disp('Frequency error = ');
disp(freq_error);
disp('Thermal diffusivity = ');
disp(diffusivity);
disp('Thermal diffusivity error = ');
disp(diffusivity_err);

fid1 = fopen( printFile, 'a' );
    naming_str=strcat(str_base,num2str(i));
    fprintf(fid1, string('\n%s %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g'), naming_str, grat, freq, freq_err, A, AErr, diff, diff_err, beta, betaErr, B, BErr, theta, thetaErr, tau(3,(i+1)), tauErr, C, CErr);
    fclose(fid1);
%close all % comment this out if you want to keep the plot windows at the end