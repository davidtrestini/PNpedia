# PNpedia

This repository was created to make lengthy post-Newtonian results:
* easier to find, by putting all known results in one place
* more user-friendly, by including associated machine-readable files (in Wolfram language format)
* more reliable, by allowing users to *endorse* a result they agree with, and raise issues when typos or disagreements have been identified

The scope of this repository is complete, fully analytical, post-Newtonian results for the two-body problem. Partial results, semi-analytical or many-body results are outside the scope of this repository.

The repository currently contains two main directories:
* "Raw data", which is compendium of publications and the machine-readable post-Newtonian results associated to them. It is itself subvided into two directories:
    * "arXiv", containing preprints on arXiv, which is used by default
    * "publications", containing peer-reviewed journal publications which either predate arXiv or have results that differ from the arXiv version 
* "Gauge-invariant quantities at future null infinity", which aims to contain machine-readable files corresponding to the state of the art of post-Newtonian predictions for gravitational waveforms. It is structured using nested directories.
In the future, I intend to include gauge-dependent post-Newtonian results, such as the equations of motions, or gauge-independent results that are not strictly speaking observables at future null infinity, such as the invariant quantities of the conservative problem.


Alongside each machine-readable file (or collection of files), there should be a README.md file. This file should contain:
* the description of the physical quantity presented
* the conventions used
* the reference(s) of the result, and potential known typos and disagreements
* a list of endorsers

Currently, the results are quite scarce, because it is intended to be a *collaborative* effort. Users are encourage to raise *pull request* to:
* upload a missing machine-readable file
* endorse an available machine-readable file
* raise an issue if they have a disagreement with an available machine-readable file
Suggestions are also welcome.
