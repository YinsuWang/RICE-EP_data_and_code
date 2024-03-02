function [HDD,CDD] = DegreeDaysGelegenis(HT,CT,Tmean,Tmax,Tmin,NY)
HDD = NY*(HT - Tmin).^2./(Tmean - Tmin)/4;
CDD = NY*(CT - Tmax).^2./(Tmax - Tmean)/4;