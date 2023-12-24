# About This Repository
This repository contains a simulator for studying **quantum rotational dynamics of linear molecules** induced by single- and/or two-color intense femtosecond (fs) laser pulses. An old version of this simulator was used to find optimal conditions of laser parameters to control molecular orientation (i.e., the head-to-tail order of linear asymmetric molecules) and make a comparison between simulated and experimental results in [[1](https://www.sciencedirect.com/science/article/abs/pii/S0009261418300095)].

# What Does This Simulator Do?
Theoretical details are explained in [[1](https://www.sciencedirect.com/science/article/abs/pii/S0009261418300095)] with some visual examples. With this simulator, quantum rotational dynamics of linear molecules induced by a combination of a near-infrared (800 nm) intense fs laser pulse and two-color (800 nm + 400 nm) intense femtosecond laser pulses can be investigated by numerically solving time-dependent [Schrödinger equation](https://en.wikipedia.org/wiki/Schrödinger_equation). Using simulated results, molecular orientation (i.e., head-to-tail order) and molecular alignment (i.e., molecular axis aligned along a specific line without considering the head-to-tail order) can be quantitatively evaluated. Results are saved in text files. 

# How to Use
1. Specify simulation parameters and save them.
2. Run main.cpp file. Use the command `bash run.sh` in Terminal.

# Keywords
- Coherent rotational Raman excitation. Non-adiabatic molecular alignment and orientation. Impulsive Raman. Intense femtosecond laser pulses.

# References
[[1](https://www.sciencedirect.com/science/article/abs/pii/S0009261418300095)] K. Sonoda, A. Iwasaki, K. Yamanouchi, and H. Hasegawa, Field-free molecular orientation of nonadiabatically aligned OCS, Chem. Phys. Lett., 693, 114-120, 2018
