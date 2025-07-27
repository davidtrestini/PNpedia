# PNpedia

This repository was created to make lengthy post-Newtonian results:
* easier to find, by putting all known results in one place
* more traceable, by collecting the various references in the literature in one place and flagging known discrepancies and typos
* more user-friendly, by including associated machine-readable files (in Wolfram language format)
* more reliable, by allowing users to *endorse* a result they agree with, and raise issues when typos or disagreements have been identified

Currently, the scope of this repository is complete, fully analytical, post-Newtonian results for the two-body problem in general relativity. Partial results, semi-analytical, many-body results or alternative theories of gravity are at this stage outside the scope of this repository. Post-Newtonian results from analytical self-force methods (i.e., at leading order in the mass-ratio) can be found the in the `PostNewtonianSelfForce` package of the Black Hole Perturbation Toolkit: [https://github.com/BlackHolePerturbationToolkit/PostNewtonianSelfForce](https://github.com/BlackHolePerturbationToolkit/PostNewtonianSelfForce).

The repository currently contains two main directories:
* ``Publications``, which is compendium of publications and the machine-readable post-Newtonian results associated to them. Publications are referenced preferably by their arXiv number, if it exists, otherwise by their journal reference.
* ``Core post-Newtonian quantities``, which aims to contain machine-readable files corresponding to the state of the art of post-Newtonian predictions for gravitational waveforms. It is structured using nested directories.

Alongside each machine-readable file (or collection of files), there should be a README.md file. This file should contain:
* the description of the physical quantity presented
* the notation used
* the reference(s) of the result, and potential known typos and disagreements
* a list of endorsers

Currently, the results are quite scarce, because it is intended to be a *collaborative* effort. Users are encourage to raise a *pull request* to:
* upload a missing machine-readable file
* endorse an available machine-readable file
* add or correct references or notations in the README.md files
* raise an *issue* if they have a disagreement with an available machine-readable file

Suggestions are also welcome.

If you have used this repository for your academic work, please acknowledge it by citing the repository. For this, just click on ``Cite this repository`` in the ``About`` section of the right-hand side of the page, and copy-paste the BibTeX file.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15002834.svg)](https://doi.org/10.5281/zenodo.15002834)


## Tree structure

Directories that do not contain any information yet are flagged as ``[empty]``.
Windows users might encounter problems due to long pathnames. A hack is to run the following command in the git bash:
```git config --global core.longpaths true```
