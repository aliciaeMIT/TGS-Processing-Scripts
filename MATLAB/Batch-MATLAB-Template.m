%%%%%%%%%%%%%%%%%%%
%Either run from directory containing all raw TGS data files for this exposure, or change directory here
dir=cd();
%%%%%%%%%%%%%%%%%%%

grat=<GRATING>; 	% pulled from metadata in first file, in microns
overlay1=<TIMESTAMP>	% pulled from metadata in first file
overlay2=<SAMPLESTAMP>	% pulled from metadata in first file

posstr=<FILENAME-POS>;
negstr=<FILENAME-NEG>;

EndTime=<END-TIME>	% automatically extracted from first file to analyze

TimeIndex=<TIME-INDEX>  % Loop through these in a script to determine the best one

TGSPhaseAnalysis(posstr,negstr,grat,2,0,EndTime,TimeIndex,overlay1,overlay2);
