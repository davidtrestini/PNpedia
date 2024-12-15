# PNpedia

This repository was created to make lengthy post-Newtonian results:
* easier to find, by putting all known results in one place
* more traceable, by collecting the various references in the literature in one place and flagging known discrepancies and typos
* more user-friendly, by including associated machine-readable files (in Wolfram language format)
* more reliable, by allowing users to *endorse* a result they agree with, and raise issues when typos or disagreements have been identified

Currently, the scope of this repository is complete, fully analytical, post-Newtonian results for the two-body problem in general relativity. Partial results, semi-analytical, many-body results or alternative theories of gravity are at this stage outside the scope of this repository.

The repository currently contains two main directories:
* ``Raw data``, which is compendium of publications and the machine-readable post-Newtonian results associated to them. It is itself subdivided into two directories:
    * ``arXiv``, containing preprints on arXiv, which is used by default
    * ``publications``, containing peer-reviewed journal publications which either pre-date arXiv or have results that differ from the arXiv version 
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

## Tree structure

```bash
.
|-- Core post-Newtonian quantities
|   |-- Fluxes
|   |   |-- Energy flux
|   |   |   |-- Circular orbits
|   |   |   |   |-- Spinning
|   |   |   |   |   |-- Nonprecessing
|   |   |   |   |   |   |-- With tidal effects
|   |   |   |   |   |   `-- Without tidal effects
|   |   |   |   |   `-- Precessing
|   |   |   |   |       |-- Without tidal effects
|   |   |   |   |       `-- With tidal effects
|   |   |   |   `-- Nonspinning
|   |   |   |-- Elliptic orbits
|   |   |   |   |-- Spinning
|   |   |   |   |   |-- Precessing
|   |   |   |   |   |   |-- Without tidal effects
|   |   |   |   |   |   `-- With tidal effects
|   |   |   |   |   `-- Nonprecessing
|   |   |   |   |       |-- Without tidal effects
|   |   |   |   |       `-- With tidal effects
|   |   |   |   `-- Nonspinning
|   |   |   `-- Hyperbolic orbits
|   |   |       |-- Spinning
|   |   |       |   |-- Nonprecessing
|   |   |       |   `-- Precessing
|   |   |       `-- Nonspinning
|   |   |-- Center-of-mass position flux
|   |   |   |-- Circular orbits
|   |   |   |   |-- Spinning
|   |   |   |   |   |-- Nonprecessing
|   |   |   |   |   `-- Precessing
|   |   |   |   `-- Nonspinning
|   |   |   |-- Elliptic orbits
|   |   |   |   |-- Spinning
|   |   |   |   |   |-- Nonprecessing
|   |   |   |   |   `-- Precessing
|   |   |   |   `-- Nonspinning
|   |   |   `-- Hyperbolic orbits
|   |   |       |-- Spinning
|   |   |       |   |-- Nonprecessing
|   |   |       |   `-- Precessing
|   |   |       `-- Nonspinning
|   |   |-- Angular momentum flux
|   |   |   |-- Circular orbits
|   |   |   |   |-- Spinning
|   |   |   |   |   |-- Precessing
|   |   |   |   |   |   |-- With tidal effects
|   |   |   |   |   |   `-- Without tidal effects
|   |   |   |   |   `-- Nonprecessing
|   |   |   |   |       |-- With tidal effects
|   |   |   |   |       `-- Without tidal effects
|   |   |   |   `-- Nonspinning
|   |   |   |-- Elliptic orbits
|   |   |   |   |-- Spinning
|   |   |   |   |   |-- Nonprecessing
|   |   |   |   |   `-- Precessing
|   |   |   |   `-- Nonspinning
|   |   |   `-- Hyperbolic orbits
|   |   |       |-- Spinning
|   |   |       |   |-- Nonprecessing
|   |   |       |   `-- Precessing
|   |   |       `-- Nonspinning
|   |   `-- Linear momentum flux
|   |       |-- Circular orbits
|   |       |   |-- Spinning
|   |       |   |   |-- Nonprecessing
|   |       |   |   `-- Precessing
|   |       |   `-- Nonspinning
|   |       |-- Elliptic orbits
|   |       |   |-- Spinning
|   |       |   |   |-- Nonprecessing
|   |       |   |   `-- Precessing
|   |       |   `-- Nonspinning
|   |       `-- Hyperbolic orbits
|   |           |-- Spinning
|   |           |   |-- Nonprecessing
|   |           |   `-- Precessing
|   |           `-- Nonspinning
|   |-- Waveform
|   |   |-- Circular orbits
|   |   |   |-- Spinning
|   |   |   |   |-- Nonprecessing
|   |   |   |   |   |-- With tidal effects
|   |   |   |   |   `-- Without tidal effects
|   |   |   |   `-- Precessing
|   |   |   |       |-- Without tidal effects
|   |   |   |       `-- With tidal effects
|   |   |   `-- Nonspinning
|   |   |-- Elliptic orbits
|   |   |   |-- Spinning
|   |   |   |   |-- Precessing
|   |   |   |   |   |-- Without tidal effects
|   |   |   |   |   `-- With tidal effects
|   |   |   |   `-- Nonprecessing
|   |   |   |       |-- Without tidal effects
|   |   |   |       `-- With tidal effects
|   |   |   `-- Nonspinning
|   |   `-- Hyperbolic orbits
|   |       |-- Spinning
|   |       |   |-- Nonprecessing
|   |       |   `-- Precessing
|   |       `-- Nonspinning
|   `-- Conserved quantities
|       |-- Energy
|       |   |-- Circular orbits
|       |   |   |-- Spinning
|       |   |   |   |-- Nonprecessing
|       |   |   |   |   |-- With tidal effects
|       |   |   |   |   `-- Without tidal effects
|       |   |   |   `-- Precessing
|       |   |   |       |-- Without tidal effects
|       |   |   |       `-- With tidal effects
|       |   |   `-- Nonspinning
|       |   |-- Elliptic orbits
|       |   |   |-- Spinning
|       |   |   |   |-- Precessing
|       |   |   |   |   |-- Without tidal effects
|       |   |   |   |   `-- With tidal effects
|       |   |   |   `-- Nonprecessing
|       |   |   |       |-- Without tidal effects
|       |   |   |       `-- With tidal effects
|       |   |   `-- Nonspinning
|       |   `-- Hyperbolic orbits
|       |       |-- Spinning
|       |       |   |-- Nonprecessing
|       |       |   `-- Precessing
|       |       `-- Nonspinning
|       |-- Center-of-mass position
|       |   |-- Circular orbits
|       |   |   |-- Spinning
|       |   |   |   |-- Nonprecessing
|       |   |   |   `-- Precessing
|       |   |   `-- Nonspinning
|       |   |-- Elliptic orbits
|       |   |   |-- Spinning
|       |   |   |   |-- Nonprecessing
|       |   |   |   `-- Precessing
|       |   |   `-- Nonspinning
|       |   `-- Hyperbolic orbits
|       |       |-- Spinning
|       |       |   |-- Nonprecessing
|       |       |   `-- Precessing
|       |       `-- Nonspinning
|       |-- Angular momentum
|       |   |-- Circular orbits
|       |   |   |-- Spinning
|       |   |   |   |-- Nonprecessing
|       |   |   |   |   |-- Nonprecessing
|       |   |   |   |   |   |-- With tidal effects
|       |   |   |   |   |   `-- Without tidal effects
|       |   |   |   |   `-- Precessing
|       |   |   |   |       |-- Without tidal effects
|       |   |   |   |       `-- With tidal effects
|       |   |   |   `-- Precessing
|       |   |   `-- Nonspinning
|       |   |-- Elliptic orbits
|       |   |   |-- Spinning
|       |   |   |   |-- Nonprecessing
|       |   |   |   `-- Precessing
|       |   |   `-- Nonspinning
|       |   `-- Hyperbolic orbits
|       |       |-- Spinning
|       |       |   |-- Nonprecessing
|       |       |   `-- Precessing
|       |       `-- Nonspinning
|       `-- Linear momentum
|           |-- Circular orbits
|           |   |-- Spinning
|           |   |   |-- Nonprecessing
|           |   |   `-- Precessing
|           |   `-- Nonspinning
|           |-- Elliptic orbits
|           |   |-- Spinning
|           |   |   |-- Nonprecessing
|           |   |   `-- Precessing
|           |   `-- Nonspinning
|           `-- Hyperbolic orbits
|               |-- Spinning
|               |   |-- Nonprecessing
|               |   `-- Precessing
|               `-- Nonspinning
`-- Raw data
    |-- publications
    `-- arXiv
        |-- 0711.0302v2
        |-- 2304.11186v4
        |-- 2304.11185v4
        `-- 2410.16373v1
```
