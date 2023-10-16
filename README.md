# phazap

Phazap is a package to efficiently post-process gravitational wave (GW) parameter estimation data to obtain the phases and polarization state of the signal at a given detector and frequency. Details on the method are presented in [Ezquiaga, Hu, Lo (2023)](https://arxiv.org/abs/2308.06616). The key module is phase.py.

This code is used for low-latency identification of strongly lensed gravitational waves via their phase consistency by measuring their distance in the detector phase space. The relevant module including the distance statistic is tension_utils.py

Phazap builds on top of the IGWN conda enviroment https://computing.docs.ligo.org/conda/environments/igwn/ which include the standard GW packages such as LALSuite and bilby.
