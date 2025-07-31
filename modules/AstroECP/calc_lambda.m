function [lambda_out] = calc_lambda(V_in)
%CALC_LAMBDA calculate the wavelength of a particle, with relativistic correction
%
% V_in = voltage in eV
%
% lambda_out = output wavelength in m
%
% e.g. 
% [lambda_out] = calc_lambda(V_in)
% where:
% V_in=20E3; %voltage of the microscrope in eV

%constants
h=6.626E-34; %Planck constant in Js
me=9.1E-31; %mass of electron in kg
e_charge=1.602E-19; %charge of an electron in C
c=2.998E8; %speed of light in m/s

%equation with relativistic correction
lambda_out = ( h / sqrt(2 * me * e_charge * V_in) ) / (sqrt(1 + e_charge * V_in /2/me/c^2)); % in m

end

