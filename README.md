# Virtual-Analog
A collection of virtual analog models in the form of MATLAB scripts.

Currently there are 4 circuit models within this repository. All of the models are implemented using a nonlinear state-space formation.

1. Fender Blackface Tone Stack - This is a full time-domain physical model of the tone stack circuit from the Fender blackface series of guitar amplifiers. It is a parametric, linear model that performs filtering of an input signal.

2. Diode Clipper - This is a full time-domain physical model of a 2 capacitor (high pass and low pass) diode clipper circuit. This script also allows the user to compare different numerical methods for discretization (trapezoid, midpoint, and BDF2).

3. Boss DS-1 Transistor Stage - This is a full time-domain physical model of the transistor stage from the Boss DS-1 circuit.

4. Marshall Guv'nor Inspired Distortion Effect - This is a guitar distortion model based on the Marshall Guv'nor circuit. It is not a full physical model, as the stages of the circuit are modeled separately. The tone stage of the circuit is a model of the Big Muff tone stage, rather than the original Marshall tone stage.
