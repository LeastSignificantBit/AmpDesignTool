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

close all;
clear all;
warning ('off', 'Octave:broadcast'); % We use some broadcasting

%% Import data %%
filename = 'BFR92A_5V10mA.s2p';
startRow = 6;

fileID = fopen(filename);
A = dlmread(filename,'', startRow,0);
fclose(fileID);

Z0 = 50;
gNF 	 = 0.1 + 0.2i; % Optimal gS for minimized noise figure. (set to NA if not available)
NFmin  = 2.1;         % Minimum noise figure of transistor, attained at gS = gNF

F_=A(:,1)*1e6;
S11_=A(:,2).*exp(i*pi/180*A(:,3));
S21_=A(:,4).*exp(i*pi/180*A(:,5));
S12_=A(:,6).*exp(i*pi/180*A(:,7));
S22_=A(:,8).*exp(i*pi/180*A(:,9));

idx = 11;
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
transistor_name = 'BFR92A';
transistor_VCE  = '5 V';
transistor_IC   = '10 mA';

fprintf ('## Designing an RF amplifier using the %s transistor ##\n',transistor_name);
fprintf ('The following parameters are used:\n');
if size(F_)>0
  fprintf('(%d S-parameters between %sHz - %sHz)\n', ...
          max(size(F_)), num2eng(min(F_)), num2eng(max(F_)));
endif        

fprintf ('VCE = %s,\t\tIC = %s,\nZ0 = %d Ohm,\t\tF0 = %sHz,\n',transistor_VCE,...
         transistor_IC, Z0, num2eng(F0));
if ~isna(gNF)
  fprintf('Gamma_opt = %s,\tNFmin = %.2f dB\n',num2str(gNF,3), NFmin);
endif
fprintf ('S11 = %s,\tS12 = %s,\nS21 = %s,\tS22 = %s.\n\n', num2str(S11,3),...
          num2str(S12,3),num2str(S21,3),num2str(S22,3));

clear filename startRow A fileID idx;

%% User Data
des_prio = 2;	% What to optimize for, See above.
set_gain = 10; 	% dB
NFmax    = 4;   % dB

%% What to display %%
plot_NF_C = 0;
plot_Gl_C = 0;
plot_Gs_C = 0;
if ~isna(gNF)
	plot_gNF =1;
endif

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
GSUmax = 10*log10(1 / (1 - abs(S11)^2));  % Max source matching gain, unilateral transistor.
GLUmax = 10*log10(1 / (1 - abs(S22)^2));  % Max load matching gain, unilateral transistor.
G0 = 10*log10(abs(S21)^2);		            % Transistor gain.
GTUmax = GSUmax + G0 + GLUmax;		        % Max Unilateral transducer gain.

% Conjugated impedances (maximizes gain at F0)
B1 = 1 +abs(S11)^2 - abs(S22)^2 + abs(delta)^2;
B2 = 1 +abs(S22)^2 - abs(S11)^2 + abs(delta)^2;
C1 = S11 - delta*conj(S22);
C2 = S22 - delta*conj(S11);
gSmax = (B1-B1/abs(B1)*sqrt(B1^2-4*abs(C1)^2))/2/C1;
gLmax = (B2-B2/abs(B2)*sqrt(B2^2-4*abs(C2)^2))/2/C2;

if ~isna(set_gain) && ~isnan(set_gain)  % Gain circles, unilateral transistor is assumed.
		Gmatch = set_gain - G0; % the sum GS and GL
    
		% Normalized gain factors
		gl_ = (10.^(linspace(GLUmax, Gmatch-GSUmax, 100) / 10)) ...
          / (10^(GLUmax/10 +0.00001));
		gs_ = (10.^((Gmatch - linspace(GLUmax, Gmatch-GSUmax, 100)) / 10)) ...
          /  (10^(GSUmax/10+0.00001));
    
		% Calculate gain circles
		Cl_ = (gl_ * conj(S22)) ./ (1 - (1 - gl_)*abs(S22)^2);
		Cs_ = (gs_ * conj(S11)) ./ (1 - (1 - gs_)*abs(S11)^2);
		Rl_ = (sqrt(1 -gl_)*(1 - abs(S22)^2)) ./ (1 - (1 - gl_)*abs(S22)^2);
		Rs_ = (sqrt(1 -gs_)*(1 - abs(S11)^2)) ./ (1 - (1 - gs_)*abs(S11)^2);
endif


if ~isna(gNF) % Circles of constant noise figure.
  rN = real((1+gNF)/(1-gNF)); % normalized noise resistance of the transistor at optimal noisefigure
  NF_ = logspace(log10(NFmin+0.0001),1,200); % vector of noise figures
  N_  = (NF_ - NFmin)/4/rN*abs(1+gNF)^2; %noise figure parameter
  Cf_ = gNF./(N_ +1);
  Rf_ = sqrt(N_.*(N_ + 1 - abs(gNF)^2))./(N_ + 1);
endif

%% Decide input and output matching impedances %%
switch des_prio
	case {1} % Maximum gain 
    fprintf ('Optimizing for maximum gain.\n\n');
		gS = gSmax;
		gL = gLmax;
	case {2} % Minimized noise, fixed gain
    fprintf ('Optimizing for low noise figure at a gain of %.1f dB.\n\n', set_gain);
		if find(Rs_ > abs(gNF - Cs_))
			gS = gNF;
			GS = 10*log10((1 - abs(gS)^2)/abs(1 - S11*gS)^2);
			gl = (10^((Gmatch - GS)/10))/(10^(GLUmax/10));
			i_maxbw = 1; 
			Cl_(i_maxbw) = (gl*conj(S22))/(1 - (1 - gl)*abs(S22)^2);
			Rl_(i_maxbw) = (sqrt(1 - gl)*(1 - abs(S22)^2))/(1 - (1 - gl)*abs(S22)^2);
			gL = (abs(Cl_(i_maxbw))-Rl_(i_maxbw)) *exp(i*angle(Cl_(i_maxbw)));
			plot_Gl_C = 1;
			plot_gNF = 0;
		else 
    			Coverlap = Rf_' > abs(abs(conj(Cf_')-Cs_)-Rs_); 
			i_minf = find(max(Coverlap'), 1, 'first');
			i_maxbw= find(Coverlap(i_minf,:), 1, 'first');
			gS = Cf_(i_minf) + Rf_(i_minf)*exp(i*angle(Cs_(i_maxbw) - Cf_(i_minf)));
			gL = (abs(Cl_(i_maxbw))-Rl_(i_maxbw)) *exp(i*angle(Cl_(i_maxbw)));
			plot_NF_C = 1;
			plot_Gl_C = 1;
			plot_Gs_C = 1;
		endif
	case {3} % Max Bandwidth, fixed Gain, unilateral transistor is assumed.
    fprintf ('Optimizing for bandwidth at a gain of %.1f dB.\n\n', set_gain);    
		% Maximize badwidth
		[~, i_maxbw] = min( abs(abs(Cl_)-Rl_).^2 + abs(abs(Cs_)-Rs_).^2);
    
		% Pick out the impedances
		gL = (abs(Cl_(i_maxbw))-Rl_(i_maxbw)) *exp(i*angle(Cl_(i_maxbw)));
		gS = (abs(Cs_(i_maxbw))-Rs_(i_maxbw)) *exp(i*angle(Cs_(i_maxbw)));
		plot_Gl_C = 1;
		plot_Gs_C = 1;
  case {4}
    fprintf ('Optimizing for gain with a NF of %.2f dB.\n\n', NFmax);
    N = (NFmax - NFmin)/4/rN*abs(1+gNF);
    Cf = gNF./(N + 1);
    Rf = sqrt(N.*(N + 1 - abs(gNF)^2))./(N + 1);
    if Rf > abs(Cf - gSmax)
	gS = gSmax;
    else % Bruteforse ...
        theta = linspace (0,2*pi,2000);
	ggNF_ = Rf*sin(theta) + i*Rf*cos(thetaa) + Cf;
	[~, i_ggNF] = max( 10*log10((1 - abs(ggNF_)^2)/abs(1 - S11*ggNF_)^2));
	gS = ggNF_(i_ggNF); 
    endif
    gL = gLmax;
  case {5}
    fprintf ('Optimizing for bandwidth with a NF of %.2f dB.\n\n', NFmax);
  case {6}
    fprintf ('Optimizing for gain with a gL of %sOhm.\n\n', num2eng(Z0*(1+gL)/(1-gL),3));
  case {7}
    fprintf ('Optimizing for noise figure with a gL of %sOhm.\n\n', num2eng(Z0*(1+gL)/(1-gL),3));
  case {8}
    fprintf ('Optimizing for bandwidth with a gL of %sOhm.\n\n', num2eng(Z0*(1+gL)/(1-gL),3));
  case {9}
    fprintf ('Optimizing for gain with a gS of %sOhm.\n\n', num2eng(Z0*(1+gS)/(1-gS),3));
  case {10}
    fprintf ('Optimizing for gain with a gS of %sOhm.\n\n', num2eng(Z0*(1+gS)/(1-gS),3));
	otherwise % Simple conjugate matching, assuming abs(S12)=0.
		gS = gSmax;
		gL = gLmax;
endswitch



% Create Matching
if size(F_) > 0
  [matchS, zS_] = create_match((1+gS)/(1-gS), Z0, F0, F_);
  [matchL, zL_] = create_match((1+gL)/(1-gL), Z0, F0, F_);
else
  [matchS, zS_] = create_match((1+gS)/(1-gS), Z0, F0);
  [matchL, zL_] = create_match((1+gL)/(1-gL), Z0, F0); 
endif

fprintf('Source circuit: %sOhm (as seen from the transistor)\n1. %s\n2. %s\n\n', ...
        num2eng(Z0*(1+gS)/(1-gS),3),matchS{1},matchS{2});
fprintf('Load circuit: %sOhm (as seen from the transistor)\n1. %s\n2. %s\n\n', ...
        num2eng(Z0*(1+gL)/(1-gL),3),matchL{1},matchL{2});
%% Calculate characteristics %%
% Gain
GL = 10*log10((1 - abs(gL)^2)/abs(1 - S22*gL)^2); % Output matching circuit gain 
GS = 10*log10((1 - abs(gS)^2)/abs(1 - S11*gS)^2); % Input matching circuit gain 
G  =  G0 + GL + GS;

fprintf('Transducer gain = %.2f dB (@ %sHz)\n',G,num2eng(F0));
fprintf('Source:     %.2f dB\n',GS);
fprintf('Transistor: %.2f dB\n',G0);
fprintf('Load:       %.2f dB\n',GL);

% Unilateral figure of merit
U = abs(S12)*abs(S12)*abs(S12)*abs(S12)/(1-abs(S11)^2)/(1-abs(S22)^2);
UErrU = 10*log10(1/(1-U)^2);
UErrL = 10*log10(1/(1+U)^2);
fprintf('Gain error introduced from the unilateral approximation:\n');
fprintf('%.3f < GT - GTU < %.3f dB\n\n',UErrL,UErrU);

%% Present frequency sweep %%
if size(F_) > 0
  % Presentation of bandwidth
  gS_ = squeeze((zS_- 1 )./(zS_ + 1));
  gL_ = squeeze((zL_- 1 )./(zL_ + 1));

  G0_ = 10*log10(abs(S21_).^2);
  GS_ = 10*log10((1 - abs(gS_).^2)./abs(1 - [S11_ S11_].*gS_).^2);
  GL_ = 10*log10((1 - abs(gL_).^2)./abs(1 - [S22_ S22_].*gL_).^2);
  % Max gain
  GSUmax_ = 10*log10(1 ./ (1 - abs(S11_).^2));
  GLUmax_ = 10*log10(1 ./ (1 - abs(S22_).^2));  
  GTUmax_ = GSUmax_ + G0_ + GLUmax_;
  
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
  text(F_(1),GTUmax_(1),'G_{TU_{max}}')
  plot(F_,GTUmax_,'k--')
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
  %warning('off');
  legend('\mu test', 'Unstable Load (1)' , 'Unstable Load (2)' , ...
         'Unstable Source (1)' , 'Unstable Source (2)');
  %warning('on');
  plot(x,[1 1],'k--');
  text(x(1),1.01,'Unconditionally Stable');
  
  plot([F0 F0],[min(u_) max(u_)],'k--')
  text(F0, min(u_),['F_{0}=' num2eng(F0) 'Hz'])
  if find(stableS == 0 || stableL == 0)
    fprintf('Stability:\n******* Warning *******\n\n    Instabilities ');
    fprintf('found\n\n***********************\n');
  else
    fprintf('Stability:\nNo instabilities found in the frequency range ');
    fprintf('%sHz - %sHz\n',num2eng(min(F_)),num2eng(max(F_)));
  endif
else % No frequency sweep
  fprintf('Stability:\n');
  
endif

%Stability
fprintf('µ-test value: %.4f\n',u);
fprintf('K = %.4f\n|Δ| = %.4f\n\n',K,abs(delta));


%% Plot a Smith Chart %%
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

if plot_NF_C 
  % plot NF circle
  plot(Rf_(i_minf)*sin(theta)+real(Cf_(i_minf)), Rf_(i_minf)*cos(theta)+imag(Cf_(i_minf)),'g--');
endif
if plot_gNF
  plot(real(gNF),imag(gNF),'c.') 
  text(real(gNF)+0.02,imag(gNF),'\Gamma_{NF_{opt}}')
endif
%Gain circles
if plot_Gl_C
  plot(Rl_(i_maxbw)*sin(theta)+real(Cl_(i_maxbw)), Rl_(i_maxbw)*cos(theta)+imag(Cl_(i_maxbw)),'b--');
endif
if plot_Gs_C
  plot(Rs_(i_maxbw)*sin(theta)+real(Cs_(i_maxbw)), Rs_(i_maxbw)*cos(theta)+imag(Cs_(i_maxbw)),'r--');
endif
% 
plot(real(gL),imag(gL),'b.')
text(real(gL)+0.02,imag(gL),'\Gamma_{L}')
plot(real(gS),imag(gS),'r.')
text(real(gS)+0.02,imag(gS),'\Gamma_{S}')
