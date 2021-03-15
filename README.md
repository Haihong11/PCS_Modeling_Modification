Insert your soft manipulator parameters in: piecewise_driver

and run it to simulate the dynamics. The differential equations are implemented in: piecewise_derivatives

The state vector is composed by one 6D vector of the section strain for each section of the manipulator. The integration quantities such as position, orientation and velocity of each cross section of the manipulator are calculated in:

piecewise_observables
