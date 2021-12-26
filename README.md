# probabilistic_axyb
This is a MATLAB code of probabilistic solver for calibration problem AX=YB, of which the detailed algorithm is presented in the paper entitled "Probabilistic Framework for Hand–Eye and Robot–World Calibration AX=YB" (under review).  
The code is uploaded for review purpose now. It will be cleaned up once the paper is accepted.

# How to use
1. See the paper for the three system noise configuration.

2-1. For noise configuration2 1 and 2, call the following function:

[X,Y,C,J,N,M] = solveAXYB_prob(A, B, X0, Y0, invSig_wN, invSig_pN, invSig_wM, invSig_pM, noiseConf, step_R, step_p)

2-2. For noise configuration 3, call the following funciton:

[X_prob, Y_prob, J] = solveAXYB_prob_noiselessA(A, B, X_prob, Y_prob, invSig_wM, invSig_pM, step_R, step_p)

Here X,Y are the calibration results, and A, B are the measurements pairs in size of 4 X 4 X n (n is the number of meausurement pairs). invSig_wN, invSig_pN, invSig_wM, invSig_pM are the inverse of rotational and translational noise covariances of A and B, each of which is in size of 3 X 3 X n.

# demos
The demo scripts are named such as "main_fig6.m", which generate the figures in the paper.

