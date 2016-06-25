#!/usr/bin/octave -q
%    Discrete Transistor RF Design Script
%    Copyright (C) 2016  Martin Berglund
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%
% This Octave/Matlab script is intended for designing RF amplifiers with 
% discrete transistors. The goal is to as closely as possible follow the 
% variable names from David M. Pozars Microwave Engineering book where 
% most of the formulas also originates from.
%
% The script uses scattering parameters to calculate the matching network
% required. The script does not handle biasing but tries to give some 
% suggested matching networks.
%
% To preform a calculation:
%
% 1. Somehow import the scattering parameters into S11,S12,S21,S22, and F.
%    Using csv import or just enter it by hand.
% 
% 2. Choose What parameter should be the primary parameter to design for
%    e.g. Gain, bandwidth or Noise. The choose the secondary parameter 
%    (see the table below).
%   
%    | des_prio | Priority | Fixed | note
%    |---------------------------------------------------------------
%    | 1        | Gain     |       | Maximum gain 
%    | 2        | Noise    |       | Minimum noise figure
%    | 3        | Bandwidth| Gain  | Maximum bandwidth at fixed gain
%    | 4        | Noise    | Gain  | Minimum noise at fixed gain
%
% 3. Enter Additional parameters that might be needed for the calculation.
%    e.g. gNF, Gset etc.
%
% 4. Run the scrip to look at the result.
% 
% 5. Adjust parameters accordingly.
% 
% 6. Iterate step #4 and #5 until the desired design is acquired.
%
%%

close all;
clear all;

%% Import data %%
filename = 'BFR92A_5V10mA.s2p';
startRow = 6;

fileID = fopen(filename);
A = dlmread(filename,'', startRow,0);
fclose(fileID);

Z0 = 50;

F_=A(:,1)*1e6;
S11_=A(:,2).*exp(i*pi/180*A(:,3));
S21_=A(:,4).*exp(i*pi/180*A(:,5));
S12_=A(:,6).*exp(i*pi/180*A(:,7));
S22_=A(:,8).*exp(i*pi/180*A(:,9));

idx = 10;
F0 = F_ (idx);
S11= S11_(idx);
S12= S12_(idx);
S21= S21_(idx);
S22= S22_(idx);

%F0 = 1.9e9;
%S11= 0.869*exp(i*pi/180* -159);
%S12= 0.031*exp(i*pi/180* -9);
%S21= 4.250*exp(i*pi/180* 61);
%S22= 0.507*exp(i*pi/180* -117);

clear filename startRow A fileID idx;

%% User Data
des_prio = 1;	% What to optimize fore, See above.
set_gain = 10; 	% [dB]
gNF 	 = NaN; % Optimal gS for minimized noise figure.

%% Calculate some useful variables %%
%Reflection coefficient looking into port 1 of the transistor with Zs=Z0
gIn = S11 + S12 * S21 / (1 - S22); 
%Reflection coefficient looking into port 2 of the transistor with Zl=Z0
gOut = S22 + S12 * S21 / (1 - S11);
% determinant of the scattering matrix
delta = S11 * S22 - S12 * S21;

% Load and source stability circles
ScL = conj(S22 - delta * conj(S11))/(abs(S22)^2 -abs(delta)^2);
ScS = conj(S11 - delta * conj(S22))/(abs(S11)^2 -abs(delta)^2);
SrL = abs(S12 * S21 / ( abs(S22)^2 - abs(delta)^2 ));
SrS = abs(S12 * S21 / ( abs(S11)^2 - abs(delta)^2 ));

% Rollet's stability condition
K = (1 - abs(S11)^2 - abs(S22)^2 + abs(delta)^2) / (2 * abs( S12 * S21));
% µ-test of stability
u = (1 - abs(S11)^2) / (abs(S22-delta*conj(S11)) + abs(S12*S21));

% Max gain
GSUmax = 10*log10(1 / (1 - abs(S11)^2)); % Max source matching gain, unilateral transistor.
GLUmax = 10*log10(1 / (1 - abs(S22)^2)); % Max load matching gain, unilateral transistor.
G0 = 10*log10(abs(S21)^2);		 % Transistor gain.
GTUmax = GSUmax + G0 + GLUmax;		 % Max Unilateral transducer gain.

B1 = 1 +abs(S11)^2 - abs(S22)^2 + abs(delta)^2;
B2 = 1 +abs(S22)^2 - abs(S11)^2 + abs(delta)^2;
C1 = S11 - delta*conj(S22);
C2 = S22 - delta*conj(S11);

% Conjugated impedances (maximizes gain at F0
gSmax = (B1-B1/abs(B1)*sqrt(B1^2-4*abs(C1)^2))/2/C1;
gLmax = (B2-B2/abs(B2)*sqrt(B2^2-4*abs(C2)^2))/2/C2;

%% Decide input and output matching impedances %%
switch des_prio
	case {1} % Maximum gain 
		gS = gSmax;
		gL = gLmax;
	case {2} % Minimized noise
		gS = gNF;
		gL = conj(gOut);
	case {3} % Max Bandwidth, fixed Gain, unilateral transistor is assumed.
		Gmatch = set_gain - G0; % the sum GS and GL
    
		% Normalized gain factors
		gl_ = 10.^(linspace(GLUmax, Gmatch-GSUmax, 200) / 10) / GLUmax;
		gs_ = 10.^((Gmatch - linspace(GLUmax, Gmatch-GSUmax, 200)) / 10) / GSUmax;
    
		% Calculate gain circles
		Cl_ = (gl_ * conj(S22)) ./ (1 - (1 - gl_)*abs(S22)^2);
		Rl_ = (sqrt(1 -gl_)*(1 - abs(S22)^2)) ./ (1 - (1 - gl_)*abs(S22)^2);
		Cs_ = (gs_ * conj(S11)) ./ (1 - (1 - gs_)*abs(S11)^2);
		Rs_ = (sqrt(1 -gs_)*(1 - abs(S11)^2)) ./ (1 - (1 - gs_)*abs(S11)^2);
    
		% Maximize badwidth
		[~, i_maxbw] = min( abs(abs(Cl_)-Rl_).^2 + abs(abs(Cs_)).^2-Rs_);
    
		% Pick out the impedances
		gL = (abs(Cl_(i_maxbw))-Rl_) *exp(i*angle(Cl_(i_maxbw)));
		gS = (abs(Cs_(i_maxbw))-Rs_) *exp(i*angle(Cs_(i_maxbw)));
	otherwise % Simple conjugate matching, assuming abs(S12)=0.
		gS = conj(gIn);
		gL = conj(gOut);
endswitch



% Create Matching
if size(F_) > 0
  [matchS, zS_] = create_match((1+gS)/(1-gS), Z0, F0, F_);
  [matchL, zL_] = create_match((1+gL)/(1-gL), Z0, F0, F_);
else
  [matchS, zS_] = create_match((1+gS)/(1-gS), Z0, F0);
  [matchL, zL_] = create_match((1+gL)/(1-gL), Z0, F0); 
endif

fprintf('Source circuit: %sOhm (as seen from the transistor)\n1. %s\n2. %s\n\n',num2eng(Z0*(1+gS)/(1-gS),3),matchS{1},matchS{2});
fprintf('Load circuit: %sOhm (as seen from the transistor)\n1. %s\n2. %s\n\n',num2eng(Z0*(1+gL)/(1-gL),3),matchL{1},matchL{2});
%% Calculate characterisics %%
% Gain
GL = 10*log10((1 - abs(gL)^2)/abs(1 - S22*gL)^2); % Output matching circuit gain 
GS = 10*log10((1 - abs(gS)^2)/abs(1 - S11*gS)^2);
G  =  G0 + GL + GS;

% Unilateral figure of merit
U = abs(S12)*abs(S12)*abs(S12)*abs(S12)/(1-abs(S11)^2)/(1-abs(S22)^2);
UErrU = 10*log10(1/(1-U)^2);
UErrL = 10*log10(1/(1+U)^2);

%% Present frequency sweep %%
if size(F_) > 0
  % Presentation of badwidth
  gS_ = squeeze((zS_- 1 )./(zS_ + 1));
  gL_ = squeeze((zL_- 1 )./(zL_ + 1));

  G0_ = 10*log10(abs(S21_).^2);
  GS_ = 10*log10((1 - abs(gS_).^2)./abs(1 - [S11_ S11_].*gS_).^2);
  GL_ = 10*log10((1 - abs(gL_).^2)./abs(1 - [S22_ S22_].*gL_).^2);
  
  G11_  = G0_ + GS_(:,1) + GL_(:,1);
  G12_  = G0_ + GS_(:,1) + GL_(:,2);
  G21_  = G0_ + GS_(:,2) + GL_(:,1);
  G22_  = G0_ + GS_(:,2) + GL_(:,2);
  
  
  figure
  x=[min(F_) max(F_)];
  plot(F_,G11_,'r','linewidth',4);
  hold on;
  plot(F_,G12_,'b','linewidth',4);
  plot(F_,G21_,'g','linewidth',4);
  plot(F_,G22_,'m','linewidth',4);
  legend('1>T>1','1>T>2','2>T>1','2>T>2')
  plot(x,[G+1 G+1],'k--',x,[G-1 G-1],'k--')
  text(F0,GTUmax,'----- G_{max} -----','horizontalalignment','center')
  text(x(1),G,'\pm1dB')
  plot([F0 F0],[min(G11_)-10 max(G11_)+5],'k--')
  text(F0, min(G11_),['F_{0}=' num2eng(F0) 'Hz'])
  
  title('Gain using four different matching networks');
  xlabel('Frequency [Hz]');
  hold off
  
  % Stability presentation
  
  delta_ = S11_.*S22_ - S12_.*S21_;
  u_ = (1 - abs(S11_).^2)./(abs(S22_-delta*conj(S11_)) + abs(S12_.*S21_));
  % Stability circles
  ScL_ = conj(S22_ - delta_.*conj(S11_))./(abs(S22_).^2-abs(delta_).^2);
  SrL_ = abs(S12_.*S21_)./(abs(S22_).^2 -abs(delta_).^2);
  ScS_ = conj(S11 - delta_.*conj(S22_))./(abs(S11_).^2-abs(delta_).^2);
  SrS_ = abs(S12_.*S21_)./(abs(S11_).^2 - abs(delta_).^2);
  stableS = (abs(ScS_ - gS_)-SrS_ > 0);
  stableL = (abs(ScL_ - gL_)-SrL_ > 0);
  
  figure
  plot(F_,u_,'linewidth',4);
  hold on;
  plot(F_(~stableL(:,1)),u_(~stableL(:,1)),'*r','markersize',10,'marker','square');
  plot(F_(~stableL(:,2)),u_(~stableL(:,2)),'*r','markersize',10,'marker','+');
  plot(F_(~stableS(:,1)),u_(~stableS(:,1)),'*r','markersize',10,'marker','d');
  plot(F_(~stableS(:,2)),u_(~stableS(:,2)),'*r','markersize',10,'marker','o');
  warning('off');
  legend('\mu test', 'Unstable Load (1)' , 'Unstable Load (2)' , 'Unstable Source (1)' , 'Unstable Source (2)');
  warning('on');
  plot(x,[1 1],'k--');
  text(x(1),1.01,'Unconditionally Stable');
  
  plot([F0 F0],[min(u_) max(u_)],'k--')
  text(F0, min(u_),['F_{0}=' num2eng(F0) 'Hz'])
  if find(stableS == 0 || stableL == 0)
    fprintf('Th');
  endif

  
endif

%% Plot a Smith Chart
figure
theta = 0:pi/50:2*pi;
x=sin(theta);y=cos(theta);
plot(x,y,'k'); axis equal; hold on;
plot(0.5*x+0.5,0.5*y,'k--');
plot([1,-1],[0,0],'k');
plot(sin(theta/4+pi)+1,cos(theta/4+pi)+1,'k--')
plot(sin(theta/4+pi)+1,-cos(theta/4+pi)-1,'k--')
plot(0,0,'k.')
xlim([-1 1])
ylim([-1 1])
title(['Smith chart Z_{0}=' num2str(Z0) '\Omega']);

%Stability circles
plot(SrL*sin(theta)+real(ScL),SrL*cos(theta)+imag(ScL),'b')
text(real(ScL), imag(ScL), 'Load Stability','Color','b')
plot(SrS*sin(theta)+real(ScS),SrS*cos(theta)+imag(ScS),'r')
text(real(ScS), imag(ScS), 'Source Stability','Color','r')

% 
plot(real(gL),imag(gL),'b.')
text(real(gL)+0.05,imag(gL),'\Gamma_{L}')
plot(real(gS),imag(gS),'r.')
text(real(gS)+0.05,imag(gS),'\Gamma_{S}')