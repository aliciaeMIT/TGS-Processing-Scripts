clear all;
close all;
clc;

%calibration_filename_base= '/Users/aliciae/MIT Dropbox/Alicia Elliott/Irradiations/2026/03/02-March-2026/pi3/raw/Tungsten_Calibration/2026-03-02/Tungsten_Calibration-2026-03-02-03.40um-spot00-';

%calibration_filename_base= "Users/aliciae/MIT Dropbox/Alicia Elliott/Irradiations/2026/03/02-March-2026/pi3/raw/Tungsten_Calibration/2026-03-02/Tungsten_Calibration-2026-03-02-03.40um-spot00";
calibration_filename_base= "/Users/aliciae/MIT Dropbox/Alicia Elliott/Irradiations/2026/03/02-March-2026/pi3/raw/Tungsten_Calibration/2026-03-02/Tungsten_Calibration-2026-03-02-03.40um-spot00";
%base_dir="/Users/aliciae/MIT Dropbox/Alicia Elliott/Irradiations/2026/03/03-March-2026/data/pi3";
%calibration_filename_base= base_dir + "/Tungsten_Calibration/2026-03-02/Tungsten_Calibration-2026-03-02-03.40um-spot00";

initial_grating_guess=3.4;
%Figure out true grating size from calibration run
%initialise calibration outputs
freq_cal = 0;
freq_err_cal = 0;
speed_cal = 0;
diff_cal = 0;
diff_err_cal = 0;
tau_cal=zeros(3,1);
tau_err_cal = 0;
A_cal = 0;
AErr_cal = 0;
beta_cal = 0;
betaErr_cal = 0;
B_cal = 0;
BErr_cal = 0;
theta_cal = 0;
thetaErr_cal = 0;
C_cal = 0;
CErr_cal = 0;
file_date_time_cal = ' ';

% Calibrated grating spacing
pos_str_cal=strcat(calibration_filename_base,'-POS-1.txt');
neg_str_cal=strcat(calibration_filename_base,'-NEG-1.txt');
[freq_cal,freq_err_cal,speed_cal,diff_cal,diff_err_cal,tau_cal(:,1),tauErr_cal,A_cal,AErr_cal,beta_cal,betaErr_cal,B_cal,BErr_cal,theta_cal,thetaErr_cal,C_cal,CErr_cal]=TGSPhaseAnalysis(pos_str_cal,neg_str_cal,initial_grating_guess,2,0,0, pos_str_cal, neg_str_cal);
grat=1E6*2665.9/freq_cal; %in um, speed of sound in tungsten is 2665.9 m/s
disp(['calculated grating size: ',num2str(grat,8),' um']);