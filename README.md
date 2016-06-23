# AmpDesignTool #
This Octave script is intended for designing RF amplifiers with 
discrete transistors. The goal is to as closely as possible follow the 
variable names from David M. Pozars Microwave Engineering book where 
most of the formulas also originates from.

The script uses scattering parameters to calculate the matching network
required. The script does not handle biasing but tries to give some 
suggested matching networks.

#### To preform a calculation ####

1. Somehow import the scattering parameters into S11,S12,S21,S22, and F.
   Using csv import or just enter it by hand.

2. Choose What parameter should be the primary parameter to design for
   e.g. gain, bandwidth or noise. The choose the secondary parameter 
   (see the table below).

   des_prio | Priority | Fixed | Note
   ---------|----------|-------|---------
   1        | Gain     |       | Maximum gain 
   2        | Noise    |       | Minimum noise figure
   3        | Bandwidth| Gain  | Maximum bandwidth at fixed gain
   4        | Noise    | Gain  | Minimum noise at fixed gain

3.  Enter Additional parameters that might be needed for the calculation.
   e.g. gNF, Gset etc.

4. Run the scrip to look at the result.
 
5. Adjust parameters accordingly.

6. Iterate step #4 and #5 until the desired design is acquired.
