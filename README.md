# Probabilistic Framework for Hand-Eye and Robot-World Calibration AX = YB
<br>

MATLAB implementation of Probabilistic Framework for Hand-Eye and Robot-World Calibration AX = YB


## Overview

This is a MATLAB code of probabilistic solver for hand-eye and robot-world calibration AX = YB, of which the detailed algorithm is presented in the paper entitled "Probabilistic Framework for Hand–Eye and Robot–World Calibration AX=YB" (under review). The algorithm imcorporates different individual noise distributions of measurements A_i and B_i, and also provides calibration uncertainty as an error covariance matrix.

The code is uploaded for review purpose now. The instruction and the code will be cleaned up once the paper is accepted.

## Instruction

1. See the three system noise configurations presented in the paper.
2. Calibration functions are different between noise configurations 1,2 and noise configuration 3.
	* For noise configuration2 1 and 2, call	
		```
		[X, Y] = solveAXYB_prob(A, B, X0, Y0, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf, step_R, step_p)
		```
	* For noise configuration 3, call		
		```
		[X, Y] = solveAXYB_prob_noiselessA(A, B, X_prob, Y_prob, invSig_wM, invSig_pM, step_R, step_p)
		```		
Here ``X, Y`` are the calibration results, and ``A, B`` are the measurements pairs in size of ``4 X 4 X n`` each (n is the number of meausurement pairs). ``invSig_wN, invSig_pN, invSig_wM, invSig_pM`` are the inverses of rotational and translational noise covariances of A and B, each of which is in size of ``3 X 3 X n``.

## Demos
The demo scripts are named as ``main_fig5.m``, ``main_fig6.mm``, ..., which generate the figures in the paper.

