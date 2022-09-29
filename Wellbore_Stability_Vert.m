%%
%%% Arian Velayati, PhD
%%%% This script is used to find min mud pressure required to avoid shear
%%%% failure of the well (borehole breakout). MC failure criterion is used
%%%% to find the min Pmud in the vertical wellbore


close; clear; clc;

%% inputs

Pp = 3080; % Pore pressure (psi)
SHmax = 6300-Pp; Shmin = 4300-Pp; % Effective stresses (psi)
UCS = 3500; % Uniaxial CS (psi)
phi = atand(0.6); % Friction angle (deg) 
TVD = 7000; % True vertical depth (ft)

%% Calculations --> Wellbore breakout in vertical wells

% Min Mw to avoid borehole shear breakout

q = (1+sind(phi))/(1-sind(phi)); % Anisotropy factor

Pwshear = Pp + (3*SHmax - Shmin - UCS)/ (1+q)  % Min mud weight pressure (psi)
Mwmin  = (Pwshear/TVD)*(8.337/0.44) % min mud weight (ppg) to avoid breakouts

%% Breakout angle (wbo)

% Input mud weight 
Pw_w = 10; % ppg
Pw = 0.052 * Pw_w * TVD; % Mud column hydrostatic pressure (psi)

tetaB = acosd((SHmax + Shmin - UCS - (1+q)*(Pw-Pp))/(2*(SHmax - Shmin)))/2;
wbo = 180 - 2*tetaB

%% Determining SHmax from wellbore breakout angle

% Input (Other inputs from previous sections)
wbo1 = 70; % Wellbore breakout angle from the image logs or calculations.

SHmax1 = Pp + (UCS + (1 + q)*(Pw - Pp) - Shmin * (1 + 2*cosd(180-wbo1)))/(1-2*cosd(180-wbo1))
SHmax_eff = SHmax1 - Pp


%% Calculations --> Wellbore tensile fracture (Breakdown pressure, Pb)Pb couold be smaller or larger than S3
% When Pb > S3 hydraulic fracture could be created 
% When Pb < S3 and PW > Pb: short tensile fracs around the borehole are
% produced that do not propagate far from the borehole  

% Inputs

Ts = 800; % Tensile strength of the rock (psi)
E = 50*1000*145; % Young's modulus (psi)
alpha = 5*10^-6; % Linear thermal expansion coefficient (ppm/C)
DelT = 0; % Temperature change (C)
v = 0.25; % Poisson ratio

% Cals

Sdelt = (E*alpha*DelT) / (1-v);

Pb = Pp - SHmax + 3 * Shmin + Ts + Sdelt % psi
Pb_ppg = Pb/TVD * (8.337/0.44)
