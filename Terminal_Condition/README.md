# Files in this folder check the terminal condition of the entrainment process:

'Periodic_Solution_JFK_I_1000lux.mat' contains the periodic solution (reference state) of the S+C3 model;

'Open_Loop_shift_notol_initial.slx' simulate the S+C3 model under reference light;

'Error_check.m' loads the Simulink file, generates initial points satisfying terminal condition and integrates state equation forward and checks state error.
