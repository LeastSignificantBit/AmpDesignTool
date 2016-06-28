# AmpDesignTool #
This Octave script is intended for designing small signal RF amplifiers 
with discrete transistors. The goal is to as closely as possible follow 
the variable names from David M. Pozars Microwave Engineering book where 
most of the formulas also originates from.

The script uses scattering parameters to calculate the matching network
required. The script does not handle biasing but tries to give some 
suggested matching networks.

#### To preform a calculation ####

1. Somehow import the scattering parameters into S11, S12, S21, S22, Z0
   and F. Using csv import or just enter it by hand.

2. Enter additional parameters that might be needed for the calculation.
   e.g. gNF, Gset, FNmin, F_, S11_ etc.

3. Choose What parameter should be the primary parameter to design for
   e.g. gain, bandwidth or noise. The choose the secondary parameter 
   (see the table below).

   des_prio | Priority | Fixed | Note
   ---------|----------|-------|---------
   1        | Gain     |       | Maximum gain 
   2        | Noise    | Gain  | Minimum noise at fixed gain
   3        | Bandwidth| Gain  | Maximum bandwidth at fixed gain
   4        | Gain     | gL    | Maximum gain at fixed load impedance
   5        | Noise    | gL    | Minimum noise at fixed load impedance
   6        | Bandwidth| gL    | Maximum bandwidth at fixed load impedance
   7        | Gain     | gS    | Maximum gain at fixed source impedance
   8        | Bandwidth| gS    | Maximum bandwidth at fixed source impedance



4. Run the scrip to look at the result.
 
5. Adjust parameters accordingly.

6. Iterate step #4 and #5 until the desired design is acquired.
