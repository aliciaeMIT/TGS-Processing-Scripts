%%%%%%%%%%%%%%%%%%%
%Either run from directory containing all raw TGS data files for this exposure, or change directory here
dir=cd();
%%%%%%%%%%%%%%%%%%%

grat=<GRATING>; %in um
overlay1=<TIMESTAMP>
overlay2=<SAMPLESTAMP>

posstr=<FILENAME-POS>;
negstr=<FILENAME-NEG>;
TGSPhaseAnalysis(posstr,negstr,grat,2,0,2e-7,overlay1,overlay2);
