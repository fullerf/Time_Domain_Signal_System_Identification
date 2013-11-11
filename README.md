Time_Domain_Signal_System_Identification
========================================

Decomposes time-domain signals with an unknown strict, proper transfer function into the Takenaka-Malmquist general basis, using adaptively found system poles. The algorithm is derived from a body of literature from Tao Qian, which he calls AFD: Adaptive Fourier Decomposition.  Visit his website for a list of papers and a matlab implementation of his software: http://www.fst.umac.mo/en/staff/fsttq.html

To my knowledge, there are no implementations or algorithms published which use AFD for signals described by strict, proper transfer functions. To be clear, this means the time domain signal is causal, decays to zero at infinity, and has a transfer function which has greater polynomial order in the denominator than in the numerator.  Such signals are common in spectroscopy.  Extending AFD to this signal domain was fairly straight-forward.  This repository, written in Julia, hosts the code that implements AFD in this signal domain.
