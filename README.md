# PNpedia

This repository was created to make lengthy post-Newtonian results:
* easier to find, by putting all known results in one place
* more traceable, by collecting the various references in the literature in one place and flagging known discrepencies and typos
* more user-friendly, by including associated machine-readable files (in Wolfram language format)
* more reliable, by allowing users to *endorse* a result they agree with, and raise issues when typos or disagreements have been identified

Currently, the scope of this repository is complete, fully analytical, post-Newtonian results for the two-body problem in general relativity. Partial results, semi-analytical, many-body results or alternative theories of gravity are at this stage outside the scope of this repository.

The repository currently contains two main directories:
* ``Raw data``, which is compendium of publications and the machine-readable post-Newtonian results associated to them. It is itself subvided into two directories:
    * ``arXiv``, containing preprints on arXiv, which is used by default
    * ``publications``, containing peer-reviewed journal publications which either predate arXiv or have results that differ from the arXiv version 
* ``Quantities at future null infinity``, which aims to contain machine-readable files corresponding to the state of the art of post-Newtonian predictions for gravitational waveforms. It is structured using nested directories.

In the future, I intend to include other post-Newtonian results, such as the equations of motions, Poincaré invariants, source moments, etc.


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

If you have used this repository for you academic work, please acknowledge it by citing the repository. For this, just click on ``Cite this repository`` in the ``About`` section of the right-hand side of the page, and copy-paste the BibTeX file.


