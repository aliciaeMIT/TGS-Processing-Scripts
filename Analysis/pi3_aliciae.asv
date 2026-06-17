clear all;
close all;
clc

% for March 2 2026 Ni into Ni 550C

% directory location
base_dir="/Users/aliciae/MIT Dropbox/Alicia Elliott/Irradiations/2026/03/02-March-2026/pi3/raw";


%calibration_filename_base= "Users/aliciae/MIT Dropbox/Alicia Elliott/Irradiations/2026/03/02-March-2026/pi3/raw/Tungsten_Calibration/2026-03-02/Tungsten_Calibration-2026-03-02-03.40um-spot00";
%calibration_filename_base= base_dir + "/Tungsten_Calibration/2026-03-02/Tungsten_Calibration-2026-03-02-03.40um-spot00";
baseline_studyname="Baseline_Ni_into_Ni_6MeV_550C";

studyname = "Ni_into_Ni_6MeV_550C";
rundate = "2026-03-02";
grating_set = "03.40um";
spotname = "spot00";


%baseline_basename = base_dir + "/" + baseline_studyname + "/";

%basename = "/Users/aliciae/MIT Dropbox/Alicia Elliott/Irradiations/2026/03";
%basename="Users/aliciae/MIT Dropbox/Alicia Elliott/Irradiations/2026/03/02-March-2026/pi3/raw";
basename = base_dir + "/" + studyname;

pname = basename + "/" + rundate; %ADJUST
pp_dir = basename + "/post";
mkdir(pp_dir);

%First part of file name and number of spots
str_base=studyname + "-" + rundate + "-" + grating_set + "-" + spotname + "-";
%Change filename modifying number here
numberOfRuns=12; %ADJUST
% 283

% Calibrated grating spacings
grat = 3.5124992; %3.5217; %in um %ADJUST
%3.5124992 um

% baseline subtract?

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

printFile = pp_dir + "/" + "PP_" + str_base + "matlab.csv";
fid1 = fopen( printFile, 'wt' );
fprintf(fid1,'%s', 'run_name,date_time,grating_value[um],SAW_speed[m/s],SAW_freq[Hz],SAW_freq_error[Hz],A[Wm^-2],A_err[Wm^-2],alpha[m^2s^-1],alpha_err[m2s-1],beta[s^0.5],beta_err[s^0.5],B[Wm^-2],B_err[Wm^-2],theta,theta_err,tau[s],tau_err[s],C[Wm^-2],C_err[Wm^-2]');
disp(printFile);
fclose(fid1);

%Run through each of the files, sweeping across gratings and peak-skips
for j=1:numberOfRuns
    close all
    pos_str=strcat(pname,"/",str_base,'POS-',num2str(j),'.txt');
    neg_str=strcat(pname,"/",str_base,'NEG-',num2str(j),'.txt');
    %[pos_path, pos_fname, pos_ext] = fileparts(pos_str)
    disp(['current run is: ', num2str(j), ' of: ', num2str(numberOfRuns)]);
    [freq,freq_err,speed,diff,diff_err,tau(:,3*(j)),tauErr, A, AErr, beta, betaErr, B, BErr, theta, thetaErr, C, CErr, file_date_time]=TGSPhaseAnalysis(pos_str,neg_str,grat,2,0);
    
    %print these to the postprocessing file
    fid1 = fopen( printFile, 'a' );
    %fprintf(fid1, string('\n%s %s %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g %.8g'), pos_str, file_date_time, grat, speed, freq, freq_err, A, AErr, diff, diff_err, beta, betaErr, B, BErr, theta, thetaErr, tau(3,3*(j)), tauErr, C, CErr);
    %fprintf(fid1, string('\n%s,%s,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g'), pos_str, file_date_time, grat, speed, freq, freq_err, A, AErr, diff, diff_err, beta, betaErr, B, BErr, theta, thetaErr, tau(3,3*(j)), tauErr, C, CErr);
    fprintf(fid1, string('\n%s,%s,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g'), pos_str, file_date_time, grat, speed, freq, freq_err, A, AErr, diff, diff_err, beta, betaErr, B, BErr, theta, thetaErr, tau(3,3*(j)), tauErr, C, CErr);
    fclose(fid1);
end

close all
disp('F L A W L E S S');
