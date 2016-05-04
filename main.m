clc; close all; clear;

load Aircraft;
nz = 1000;
dz = b/2/nz;
z = 0:dz:b/2;

% Get wx, wy


% Airfoil section properties
nacaDigits = 2415;
chordLength = 1.5;
[WingSection] = GetSectionProperties(nacaDigits, chordLength);

