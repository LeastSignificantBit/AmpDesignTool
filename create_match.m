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
% zIn   - The sought normalized impedance of the matching network looking into it.
% Z0    - Normalizing impedance.
% F0    - Design frequency of the matching network.
% F_    - A vector with all frequencies for witch zIn_ is calculated.
% zOut  - The normalized impedance of the output from the network
% zOut_ - If zOut is complex a frequency dependent zOut  Needs to be supplied 
%         in order to return an accurate zIn_.
%
% Explanation of inputs
% build - Two cell strings that describe the two network solutions.
% zIn_  - The Impedance for each frequency F_ and each network.
%%

function [build, zIn_] = create_match(zIn, Z0, F0, F_=F0, zOut=1, zOut_=zOut)
  
  if imag(zOut) == 0
    Z_(:,:,1) = ones(size(F_))*Z0*zOut;
    Z_(:,:,2) = ones(size(F_))*Z0*zOut;
    else
      Z_(:,:,1) = ones(size(F_))*Z0.*zOut_;
      Z_(:,:,2) = ones(size(F_))*Z0.*zOut_;
  endif
  
  omega_ = 2*pi*F_;
  omega  = 2*pi*F0;

  if real(zIn)<real(zOut) % Load is parallel with one element
  
    Rs = real(zIn)*Z0;
    
    Rp = ((imag(zOut)/real(zOut))^2 + 1)*real(zOut)*Z0;
    Q= sqrt(Rp/Rs - 1);
    
    if Rp == real(zOut)*Z0
      Xp = Rp/Q;
    else
      Xp_par = Rp*real(zOut)/imag(zOut);
      Xp = Rp*Xp_par./(Q*Xp_par+[-1 1]*Rp);
    endif    
     
    Xs = Q*Rs+[-1 1]*imag(zIn)*Z0;
    
    Xp = [1 -1].*Xp;
    Xs = [-1 1].*Xs;
    build{1} = build{2} = [num2eng(Z0.*zOut,4) 'Ohm'];
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
    
    
  elseif real(zIn)>=real(zOut) % Load is in series with one element

    Rs = real(zOut)*Z0;
    
    Rp = ((imag(zIn)/real(zIn))^2 + 1)*real(zIn)*Z0;
    Q= sqrt(Rp/Rs - 1);
    
    if Rp == real(zIn)*Z0
      Xp = Rp/Q;
    else
      Xp_par = Rp*real(zIn)/imag(zIn);
      Xp = Rp*Xp_par./(Q*Xp_par+[1 -1]*Rp);
    endif    
     
    Xs = Q*Rs+[1 -1]*imag(zOut)*Z0;
    
    Xp = [1 -1].*Xp;
    Xs = [-1 1].*Xs;
    build{1} = build{2} = [num2eng(Z0.*zOut,4) 'Ohm'];
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
  zIn_ = Z_/Z0;
endfunction


