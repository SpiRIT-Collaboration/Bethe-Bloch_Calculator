# Energy deposit by heavy charged particles
There are several types of empirical formulae to describe dE/dx.  
Reffered from [PDG 2018 Reviews, chapter 33. Passage of particles through matter](http://pdg.lbl.gov/2018/reviews/rpp2018-rev-passage-particles-matter.pdf).

For a simple summary, please look at EnergyDepositStudy.pdf in this directory.

## MassEstimator.h
This file contains functions to fit dE/dx curve as a function of particle rigidity.  
Examples of root finder algorithm are described [here](https://root.cern.ch/root-finder-algorithms).


## draw.C
This macro draws Bethe-Bloch equation and Landau-Vavilov distribution.

## testRootFinder.C
This macro calculate the mass from artificially prepared parameters (dedx, momentum, ...).  
