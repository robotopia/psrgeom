# PSRGEOM
## A numerical simulation software package for the pulsar carousel model

This software is used in an academic paper treating the effect of pulsar emission geometry on subpulse drifting assuming the carousel model (a link to the paper will be placed here when available). Apart from instructions on how to install this software, this README only contains examples of usage that pertain directly to the paper itself.

[![DOI](https://zenodo.org/badge/158472781.svg)](https://zenodo.org/badge/latestdoi/158472781)

### Installation

This software has only been tested on Linux 4.18.12 (Arch Linux distro).

#### Prerequisites:

* [nmead](https://github.com/robotopia/nmead) - Nelder-Mead optimisation
* [emmt/newuoa](https://github.com/emmt/Algorithms) - Various numerical algorithms
* [fftw3](http://www.fftw.org/) - Fast Fourier Transform

#### Versions:

The software versions used in the paper are the following:

* psrgeom (this package): v1.4a.5 (git commit 0c6ee7c)
* nmead: v1.0.0 (git commit c554e9c)
* emmt/newuoa: (git commit a92ce79)
* fftw3: v3.3.8-1

#### Installation instructions:

```bash
$ make
$ sudo make install
```

### Usage examples

The following examples are used in the linked paper:

* Fig 4: `psr_pulsestack -a 8.5 -z 13.0 -P 1.292241446862 -4 -143.44 -N 10 -s 0.5 -S 6 -l 121 -n 1024`
* Fig 5:
```bash
psr_mcpa -a 9.0 -b 4.5 -P -1.292241446862 -p 0:360:361 -s 0.2 -L -1
psr_mcpa -a 9.0 -b 4.5 -P -1.292241446862 -p 0:360:361 -s 0.4 -L -1
psr_mcpa -a 9.0 -b 4.5 -P -1.292241446862 -p 0:360:361 -s 0.6 -L -1
psr_mcpa -a 9.0 -b 4.5 -P -1.292241446862 -p 0:360:361 -s 0.8 -L -1
psr_mcpa -a 9.0 -b 4.5 -P -1.292241446862 -p 0:360:361 -s 1.0 -L -1
```
* Fig 6 (top): `psr_pulsestack -a 9.0 -z 13.5 -P 1.292241446862 -4 -143.44 -N 10 -s 0.5 -S 6 -e B0809+74_profile.txt -l 121 -n 1024`
* Fig 6 (bottom):
```bash
psr_pulsestack -a 47 -z 48.9 -P 2.0743770781 -4 -6.65 -N 5 -s 1.0 -S 20 -e B2034+19_325MHz_profile_full_left.txt -l 100 -p 2000
psr_pulsestack -a 47 -z 48.9 -P 2.0743770781 -4 -6.65 -N 5 -s 0.8 -S 20 -e B2034+19_325MHz_profile_full_right.txt -l 101 -p 2000
```
* Fig 7 (left): `psr_pulsestack -a 8.39 -z 13.38 -P 1.292241446862 -4 -143.8320943 -N 10 -s 0.5 -S 6 -e B0809+74_profile.txt -l 121 -n 1024`
* Fig 7 (right): `psr_pulsestack -a 8.39 -z 13.38 -P 1.292241446862 -4 -143.8320943 -N 10 -s 0.25 -S 6 -e B0809+74_profile.txt -l 121 -n 1024`
* Fig 8:
```bash
psr_mcpa -a 9 -z 13.5 -P 1.292241446862 -p -110:100:360 -s 0.2 -1 -d
psr_mcpa -a 9 -z 13.5 -P 1.292241446862 -p -110:100:360 -s 0.4 -1 -d
psr_mcpa -a 9 -z 13.5 -P 1.292241446862 -p -110:100:360 -s 0.6 -1 -d
psr_mcpa -a 9 -z 13.5 -P 1.292241446862 -p -110:100:360 -s 0.8 -1 -d
psr_mcpa -a 9 -z 13.5 -P 1.292241446862 -p -110:100:360 -s 1.0 -1 -d
psr_mcpa -a 9 -z 13.5 -P 1.292241446862 -p -110:100:360 -s 0.2 -1
psr_mcpa -a 9 -z 13.5 -P 1.292241446862 -p -110:100:360 -s 0.4 -1
psr_mcpa -a 9 -z 13.5 -P 1.292241446862 -p -110:100:360 -s 0.6 -1
psr_mcpa -a 9 -z 13.5 -P 1.292241446862 -p -110:100:360 -s 0.8 -1
psr_mcpa -a 9 -z 13.5 -P 1.292241446862 -p -110:100:360 -s 1.0 -1
```

Other dependecies:

* [B0809+74_profile.txt](http://www.epta.eu.org/epndb/ascii/nsk+15/J0814+7429/B0809+74_L78237.txt)
* [B2034+19_325MHz_profile_full_left.txt and B2034+19_325MHz_profile_full_right.txt](http://adsabs.harvard.edu/abs/2017JApA...38...53R) -- from Fig 1.
