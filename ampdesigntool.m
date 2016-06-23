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

%% Import data
filename = 'BFR92A_5V10mA.s2p';
startRow = 6;

fileID = fopen(filename);
A = dlmread(filename,'', startRow,0);
fclose(fileID);

F=A(:,1);
S11=A(:,2).*exp(i*pi/180*A(:,3));
S21=A(:,4).*exp(i*pi/180*A(:,5));
S12=A(:,6).*exp(i*pi/180*A(:,7));
S22=A(:,8).*exp(i*pi/180*A(:,9));
clear filename startRow A fileID;

%% User Data
des_prio = 1;
set_gain = 10^(14/10);

%% Calculate some useful variables
