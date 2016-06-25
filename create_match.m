%    Discrete Transistor RF Design Script (create_match)
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
% A small scrip that matches two impedances using lumped elements
% and also returns the apparent impedance over a frequency span F_.
%
% This script was purpose build for use together with the Discrete Transistor
% RF Design Script but can be used to find matching networks for other circuits.
% 
% Most of the theory is taken from chapter 4 in RF Circuit Design, Second Edition 
% by Christopher Bowick.
%
% I also would like to complement this with a stub matching circuit wit the same
% functionality.
%
% Explanation of inputs
% gIn   - The sought scattering parameter of the matching network looking into it.
% Z0    - Normalizing impedance for the scattering parameters.
% F0    - Design frequency of the matching network.
% F_    - A vector with all frequencies for witch gIn_ is calculated.
% gOut  - The scattering parameter of the output from the network
% gOut_ - If gOut is complex a frequency dependent gOut  Needs to be supplied 
%         in order to return an accurate gIn_.
%
% Explanation of inputs
% build - Two cell strings that describe the two network solutions.
% gIn_  - The scattering parameter for each frequency F_ and each network.
%%

function [build, gIn_] = create_match(gIn, Z0, F0, F_=F0, gOut=1, gOut_=gOut)
  
  if imag(gOut) == 0
    Z_(:,:,1) = ones(size(F_))*Z0*gOut;
    Z_(:,:,2) = ones(size(F_))*Z0*gOut;
    else
      Z_(:,:,1) = ones(size(F_))*Z0.*gOut_;
      Z_(:,:,2) = ones(size(F_))*Z0.*gOut_;
  endif
  
  omega_ = 2*pi*F_;
  omega  = 2*pi*F0;

  if real(gIn)<real(gOut) % Load is parallel with one element
  
    Rs = real(gIn)*Z0;
    
    Rp = ((imag(gOut)/real(gOut))^2 + 1)*real(gOut)*Z0;
    Q= sqrt(Rp/Rs - 1);
    
    if Rp == real(gOut)*Z0
      Xp = Rp/Q;
    else
      Xp_par = Rp*real(gOut)/imag(gOut);
      Xp = Rp*Xp_par./(Q*Xp_par+[-1 1]*Rp);
    endif    
     
    Xs = Q*Rs+[-1 1]*imag(gIn)*Z0;
    
    Xp = [1 -1].*Xp;
    Xs = [-1 1].*Xs;
    build{1} = build{2} = [num2eng(Z0.*gOut,4) 'Ohm'];;
    for ix = [1 2]
         
      if Xp(ix) < 0     % Capacitor
        C  = -1/(omega*Xp(ix));
        Zshunt_=1./(i*C*omega_);
        Z_(:,:,ix) = Zshunt_.*Z_(:,:,ix)./(Zshunt_ + Z_(:,:,ix));
        build{ix} = [ ' shunt ' num2eng(C,3) 'F |' build{ix}];
      elseif Xp(ix) > 0 % Inductor
        L  = Xp(ix)/omega;
        Zshunt_ = i*omega_*L;
        Z_(:,:,ix) = Zshunt_.*Z_(:,:,ix)./(Zshunt_ + Z_(:,:,ix));
        build{ix} = [' shunt ' num2eng(L,3) 'H |' build{ix}];
      endif
       
      if Xs(ix) < 0     % Capacitor
        C  = -1/(omega*Xs(ix));
        Z_(:,:,ix) = Z_(:,:,ix)+1./(i*C*omega_);
        build{ix} = [ 'series ' num2eng(C,3) 'F |' build{ix}];
      
      elseif Xs(ix) > 0 % Inductor
        L  = Xs(ix)/omega;
        Z_(:,:,ix) = Z_(:,:,ix)+i*omega_*L;
        build{ix} = [ 'series ' num2eng(L,3) 'H |' build{ix}];
      endif

      build{ix} = ['-->| ' build{ix}];
    endfor
    
    
  elseif real(gIn)>=real(gOut) % Load is in series with one element

    Rs = real(gOut)*Z0;
    
    Rp = ((imag(gIn)/real(gIn))^2 + 1)*real(gIn)*Z0;
    Q= sqrt(Rp/Rs - 1);
    
    if Rp == real(gIn)*Z0
      Xp = Rp/Q;
    else
      Xp_par = Rp*real(gIn)/imag(gIn);
      Xp = Rp*Xp_par./(Q*Xp_par+[1 -1]*Rp);
    endif    
     
    Xs = Q*Rs+[1 -1]*imag(gOut)*Z0;
    
    Xp = [1 -1].*Xp;
    Xs = [-1 1].*Xs;
    build{1} = build{2} = [num2eng(Z0.*gIn,4) 'Ohm'];;
    for ix = [1 2]
    
      if Xs(ix) < 0     % Capacitor
        C  = -1/(omega*Xs(ix));
        Z_(:,:,ix) = Z_(:,:,ix)+1./(i*C*omega_);
        build{ix} = [ 'series ' num2eng(C,3) 'F |' build{ix}];
      
      elseif Xs(ix) > 0 % Inductor
        L  = Xs(ix)/omega;
        Z_(:,:,ix) = Z_(:,:,ix)+i*omega_*L;
        build{ix} = [ 'series ' num2eng(L,3) 'H |' build{ix}];
      endif
      
      if Xp(ix) < 0     % Capacitor
        C  = -1/(omega*Xp(ix));
        Zshunt_=1./(i*C*omega_);
        Z_(:,:,ix) = Zshunt_.*Z_(:,:,ix)./(Zshunt_ + Z_(:,:,ix));
        build{ix} = [ ' shunt ' num2eng(C,3) 'F |' build{ix}];
      elseif Xp(ix) > 0 % Inductor
        L  = Xp(ix)/omega;
        Zshunt_ = i*omega_*L;
        Z_(:,:,ix) = Zshunt_.*Z_(:,:,ix)./(Zshunt_ + Z_(:,:,ix));
        build{ix} = [' shunt ' num2eng(L,3) 'H |' build{ix}];
      endif
       

      build{ix} = ['-->| ' build{ix}];
    endfor
    
  endif 
  gIn_ = Z_/Z0;
endfunction

function s=num2eng(d, prec=4)
  if d>=1e9
    s=[num2str(d/1e9,prec) ' G'];
  elseif d>=1e6
    s=[num2str(d/1e6,prec) ' M'];
  elseif d>=1e3
    s=[num2str(d/1e3,prec) ' k'];
  elseif d>=1
    s=[num2str(d,prec) ' '];
  elseif d>=1e-3
    s=[num2str(d/1e-3,prec) ' m'];
  elseif d>=1e-6
    s=[num2str(d/1e-6,prec) ' µ'];
  elseif d>=1e-9
    s=[num2str(d/1e-9,prec) ' n'];
  elseif d>=1e-12
    s=[num2str(d/1e-12,prec) ' p'];
  elseif d>=1e-15
    s=[num2str(d/1e-15,prec) ' f'];
  else
    s=[num2str(d/1e-18,prec) ' a'];
  endif
endfunction
