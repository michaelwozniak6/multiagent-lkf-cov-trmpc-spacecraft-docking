==========================================================================
Wozniak, Michael P
AAE568 - Final Project
Optimal Task Assignment and Robust Tube-Based Formation Control for 
Multi-Agent Spacecraft Proximity Operations & Docking
==========================================================================
--------------------------------------------------------------------------
File Structure:
--------------------------------------------------------------------------
- wozniak_aae568_finalproject_fullintegration.m
	Main implementation, fully simulates KF, COV for task 
	assignment, and tube-based MPC
	**Primary algorithm to run**

- wozniak_aae568_finalproject_mpctube_QRtuning.m
	Tuning algorithm that performs a grid search of the Q,R
	diagonal matrices, observes the number of steps where the 
	control constraints were active, and returns the gains with
	minimum control saturation (informs TRMPC implementation)

- wozniak_aae568_finalproject_nonlinearAttempt.m
	Extension to main implementation that attempts to use nonlinear
	dynamics. Warning: this algorithm does not converge and is not
	recommended. However, it showcases the full structural 
	implementation of nonlinear dynamics and demonstrates the mission
	phases where the system performance is no longer sufficient.
	This file would be a starting point for future work in 
	bridging the sim-to-real gap.
--------------------------------------------------------------------------
How to Run:
--------------------------------------------------------------------------
Case I - Default Run, Linear Dynamics (Most Recommended)
- Open wozniak_aae568_finalproject_fullintegration.m and click Run
	By default, this considers linear dynamics and converges

Case II - Tuning Q,R Gains then Running
- Open wozniak_aae568_finalproject_mpctube_QRtuning.m and click Run
	The Q, R values at each stage of the search will be printed, and
	the final Q*, R* will be printed at the end.
	Modify qGain, rGain if a different search space is desired

- Open wozniak_aae568_finalproject_fullintegration.m and navigate to
	the function prepDockingQP
	Define qBar, rBar to the Q*, R* values returned from the tuning
	algorithm
	Run wozniak_aae568_finalproject_fullintegration.m

Case III - Attempt Running with Nonlinear Dynamics
- Open wozniak_aae568_finalproject_nonlinearAttempt.m and click Run
--------------------------------------------------------------------------
Contact:
--------------------------------------------------------------------------
Michael Wozniak
woznia14@purdue.edu | michaelwozniak14@gmail.com
==========================================================================