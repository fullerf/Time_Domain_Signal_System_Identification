Time_Domain_Signal_System_Identification
========================================

Decomposes time-domain signals with an unknown strict, proper transfer function into the Takenaka-Malmquist general basis, using adaptively found system poles. The algorithm is derived from a body of literature from Tao Qian, which he calls AFD: Adaptive Fourier Decomposition.  Visit his website for a list of papers and a matlab implementation of his software: http://www.fst.umac.mo/en/staff/fsttq.html

To my knowledge, there are no implementations or algorithms published which use AFD for signals described by strict, proper transfer functions. To be clear, this means the time domain signal is causal, decays to zero at infinity, and has a transfer function which has greater polynomial order in the denominator than in the numerator.  Such signals are common in spectroscopy.  Extending AFD to this signal domain was fairly straight-forward.  This repository, written in Julia, hosts the code that implements AFD in this signal domain.

Installation:

At the moment, what I have is a small pile of functions wrapped in a module.  Clone this directory and then use the commands:
include("Path_TO_YOUR_CLONE/AFDM.jl")
using AFDM

On windows machines, in IJulia or the Forio IDE, you need to use linux style / slashes in the path, rather than the windows style \ slashes. Or you need to use double \ instead of \.  Either works.  Once you've included the funciton pile,
it gives you access to the following main functions:

1. afd(xn,Nmax,Nz): Takes a one-dimensional time domain measurement xn, which is assumed (for now) to be sampled on an even time grid.  It will then iterate Nmax times over about Nz candidate poles within the unit disk.  It returns the Nmax poles and their amplitudes.
2. reconstructTimeDomain(poles,amps,t): Takes a vector poles and their associated amplitudes and reconstructs the time domain signal at a requested time index t. t can be a float valued, which means that the reconstructed time domain signal allows you to interpolate between sampled time points.
3. czt(xn,A,W,k): The chirp z transform, evaluating the z-transform of xn along a spiral starting at A and continuing every W^n thereafter for k points.
4. UnitDiskGrid(Nz,delta): Produces about Nz points inside the unit disk, ranging from delta to 1-delta.
