#!/bin/bash

var_input="%%  Input and Output"

var_extra="%-------------------------------------------------------
ICFormat                  0     % O)ASCII 1)Gadget 2)User defined
SnapshotFileBase        _ph1

%-------------------------------------------------------
% Tree related options
SpatialScale            1.0 %  x->x/SpatialScale and v->v
PartBoundary            7   % Min particles in a node to do boundary correction 
NodeSplittingCriterion  1   % 0)Alternate 1) Min Entropy 
CubicCells              0   % use 1 in spherically symmetric systems
MedianSplittingOn       0   %

%--------------------------------------------------------
% Smoothing options  AM=adaptive metric Ker=Kernel Sp=Spherical Pr=Product 
% 0) None 1)FiEstAS 2)Ker Sp Normal 3)Ker Sp AM 4)KerPr Normal 5)KerPr AM
TypeOfSmoothing      0
DesNumNgb            50   % 2-10 for Fiestas and 25-100 for Kernel
VolCorr              1    %  0) Disbale 1) Enable 

%--------------------------------------------------------
% Kernel smoothing  related options
% 0) B-Spline 1)top hat 2)Bi_weight (1-x^2)^2 3)Epanechikov 4)CIC 5)TSC 
TypeOfKernel           3 
KernelBiasCorrection   1    % 0)none 1)shift central particle 
AnisotropicKernel      0    % 0) Isotropic 1) Anisotropic 
Anisotropy             0    % fix minimum c/a minor/major axis ratio 
DesNumNgbA             128   % Neighbors for cal covar metric for Anisotropic Ker 
%--------------------------------------------------------
% other miscellaneous option
TypeListOn        0
PeriodicBoundaryOn 0
%--------------------------------------------------------"

destdir="/home/brian/Desktop/second_project_gcs/Enbid-2.0/parameterfiles/"

for f in /home/brian/Desktop/second_project_gcs/data/enbid_files/*.ascii;
do
	input=$f
	substring=$(basename $input .ascii)
	ascii=".ascii"

	#creates a parameterfile for each snapshot/Nclusters
	mkdir "$destdir$substring"

	#fills in the parameter file
	echo "$var_input" >> "$parameterfile"

	var_initcond="InitCondFile	"
	var_filename="$destdir$substring$ascii"

	echo "$var_initcond$var_filename" >> "$parameterfile"
	echo "$var_extra" >> "$parameterfile"

	#./Enbid "$parameterfile"

done
