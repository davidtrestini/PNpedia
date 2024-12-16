# PNpedia

This repository was created to make lengthy post-Newtonian results:
* easier to find, by putting all known results in one place
* more traceable, by collecting the various references in the literature in one place and flagging known discrepancies and typos
* more user-friendly, by including associated machine-readable files (in Wolfram language format)
* more reliable, by allowing users to *endorse* a result they agree with, and raise issues when typos or disagreements have been identified

Currently, the scope of this repository is complete, fully analytical, post-Newtonian results for the two-body problem in general relativity. Partial results, semi-analytical, many-body results or alternative theories of gravity are at this stage outside the scope of this repository.

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

## Tree structure

Directories that do not contain any information yet are flagged as ``[empty]``. The current tree structure is as follows:

```bash
.
|-- Core post-Newtonian quantities
|   |-- Conserved quantities
|   |   |-- Angular momentum [empty]
|   |   |   |-- Circular orbits [empty]
|   |   |   |   |-- Nonspinning [empty]
|   |   |   |   `-- Spinning [empty]
|   |   |   |       |-- Precessing [empty]
|   |   |   |       |   |-- With tidal effects [empty]
|   |   |   |       |   `-- Without tidal effects [empty]
|   |   |   |       `-- Nonprecessing [empty]
|   |   |   |           |-- With tidal effects [empty]
|   |   |   |           `-- Without tidal effects [empty]
|   |   |   |-- Elliptic orbits [empty]
|   |   |   |   |-- Nonspinning [empty]
|   |   |   |   `-- Spinning [empty]
|   |   |   |       |-- Nonprecessing [empty]
|   |   |   |       |   |-- With tidal effects [empty]
|   |   |   |       |   `-- Without tidal effects [empty]
|   |   |   |       `-- Precessing [empty]
|   |   |   |           |-- With tidal effects [empty]
|   |   |   |           `-- Without tidal effects [empty]
|   |   |   `-- Hyperbolic orbits [empty]
|   |   |       |-- Nonspinning [empty]
|   |   |       `-- Spinning [empty]
|   |   |           |-- Nonprecessing [empty]
|   |   |           |   |-- With tidal effects [empty]
|   |   |           |   `-- Without tidal effects [empty]
|   |   |           `-- Precessing [empty]
|   |   |               |-- With tidal effects [empty]
|   |   |               `-- Without tidal effects [empty]
|   |   |-- Center-of-mass position [empty]
|   |   |   |-- Circular orbits [empty]
|   |   |   |   |-- Nonspinning [empty]
|   |   |   |   `-- Spinning [empty]
|   |   |   |       |-- Precessing [empty]
|   |   |   |       |   |-- With tidal effects [empty]
|   |   |   |       |   `-- Without tidal effects [empty]
|   |   |   |       `-- Nonprecessing [empty]
|   |   |   |           |-- With tidal effects [empty]
|   |   |   |           `-- Without tidal effects [empty]
|   |   |   |-- Elliptic orbits [empty]
|   |   |   |   |-- Nonspinning [empty]
|   |   |   |   `-- Spinning [empty]
|   |   |   |       |-- Nonprecessing [empty]
|   |   |   |       |   |-- With tidal effects [empty]
|   |   |   |       |   `-- Without tidal effects [empty]
|   |   |   |       `-- Precessing [empty]
|   |   |   |           |-- With tidal effects [empty]
|   |   |   |           `-- Without tidal effects [empty]
|   |   |   `-- Hyperbolic orbits [empty]
|   |   |       |-- Nonspinning [empty]
|   |   |       `-- Spinning [empty]
|   |   |           |-- Nonprecessing [empty]
|   |   |           |   |-- With tidal effects [empty]
|   |   |           |   `-- Without tidal effects [empty]
|   |   |           `-- Precessing [empty]
|   |   |               |-- With tidal effects [empty]
|   |   |               `-- Without tidal effects [empty]
|   |   |-- Linear momentum [empty]
|   |   |   |-- Circular orbits [empty]
|   |   |   |   |-- Nonspinning [empty]
|   |   |   |   `-- Spinning [empty]
|   |   |   |       |-- Precessing [empty]
|   |   |   |       |   |-- With tidal effects [empty]
|   |   |   |       |   `-- Without tidal effects [empty]
|   |   |   |       `-- Nonprecessing [empty]
|   |   |   |           |-- With tidal effects [empty]
|   |   |   |           `-- Without tidal effects [empty]
|   |   |   |-- Elliptic orbits [empty]
|   |   |   |   |-- Nonspinning [empty]
|   |   |   |   `-- Spinning [empty]
|   |   |   |       |-- Nonprecessing [empty]
|   |   |   |       |   |-- With tidal effects [empty]
|   |   |   |       |   `-- Without tidal effects [empty]
|   |   |   |       `-- Precessing [empty]
|   |   |   |           |-- With tidal effects [empty]
|   |   |   |           `-- Without tidal effects [empty]
|   |   |   `-- Hyperbolic orbits [empty]
|   |   |       |-- Nonspinning [empty]
|   |   |       `-- Spinning [empty]
|   |   |           |-- Nonprecessing [empty]
|   |   |           |   |-- With tidal effects [empty]
|   |   |           |   `-- Without tidal effects [empty]
|   |   |           `-- Precessing [empty]
|   |   |               |-- With tidal effects [empty]
|   |   |               `-- Without tidal effects [empty]
|   |   `-- Energy
|   |       |-- Circular orbits
|   |       |   |-- Nonspinning
|   |       |   `-- Spinning
|   |       |       |-- Nonprecessing
|   |       |       |   |-- With tidal effects
|   |       |       |   `-- Without tidal effects
|   |       |       `-- Precessing [empty]
|   |       |           |-- Without tidal effects [empty]
|   |       |           `-- With tidal effects [empty]
|   |       |-- Elliptic orbits [empty]
|   |       |   |-- Nonspinning [empty]
|   |       |   `-- Spinning [empty]
|   |       |       |-- Nonprecessing [empty]
|   |       |       |   |-- With tidal effects [empty]
|   |       |       |   `-- Without tidal effects [empty]
|   |       |       `-- Precessing [empty]
|   |       |           |-- With tidal effects [empty]
|   |       |           `-- Without tidal effects [empty]
|   |       `-- Hyperbolic orbits [empty]
|   |           |-- Nonspinning [empty]
|   |           `-- Spinning [empty]
|   |               |-- Nonprecessing [empty]
|   |               |   |-- With tidal effects [empty]
|   |               |   `-- Without tidal effects [empty]
|   |               `-- Precessing [empty]
|   |                   |-- With tidal effects [empty]
|   |                   `-- Without tidal effects [empty]
|   |-- Fluxes
|   |   |-- Center-of-mass position flux [empty]
|   |   |   |-- Circular orbits [empty]
|   |   |   |   |-- Spinning [empty]
|   |   |   |   |   |-- Nonprecessing
|   |   |   |   |   `-- Precessing
|   |   |   |   `-- Nonspinning [empty]
|   |   |   |-- Elliptic orbits [empty]
|   |   |   |   |-- Spinning [empty]
|   |   |   |   |   |-- Nonprecessing
|   |   |   |   |   `-- Precessing
|   |   |   |   `-- Nonspinning [empty]
|   |   |   `-- Hyperbolic orbits [empty]
|   |   |       |-- Spinning [empty]
|   |   |       |   |-- Nonprecessing
|   |   |       |   `-- Precessing
|   |   |       `-- Nonspinning [empty]
|   |   |-- Linear momentum flux [empty]
|   |   |   |-- Circular orbits [empty]
|   |   |   |   |-- Spinning [empty]
|   |   |   |   |   |-- Nonprecessing
|   |   |   |   |   `-- Precessing
|   |   |   |   `-- Nonspinning [empty]
|   |   |   |-- Elliptic orbits [empty]
|   |   |   |   |-- Spinning [empty]
|   |   |   |   |   |-- Nonprecessing
|   |   |   |   |   `-- Precessing
|   |   |   |   `-- Nonspinning [empty]
|   |   |   `-- Hyperbolic orbits [empty]
|   |   |       |-- Spinning [empty]
|   |   |       |   |-- Nonprecessing
|   |   |       |   `-- Precessing
|   |   |       `-- Nonspinning [empty]
|   |   |-- Energy flux
|   |   |   |-- Circular orbits
|   |   |   |   |-- Nonspinning
|   |   |   |   `-- Spinning
|   |   |   |       |-- Nonprecessing
|   |   |   |       |   |-- With tidal effects
|   |   |   |       |   `-- Without tidal effects
|   |   |   |       `-- Precessing [empty]
|   |   |   |           |-- Without tidal effects [empty]
|   |   |   |           `-- With tidal effects [empty]
|   |   |   |-- Elliptic orbits
|   |   |   |   |-- Nonspinning
|   |   |   |   `-- Spinning [empty]
|   |   |   |       |-- Nonprecessing [empty]
|   |   |   |       |   |-- Without tidal effects [empty]
|   |   |   |       |   `-- With tidal effects [empty]
|   |   |   |       `-- Precessing [empty]
|   |   |   |           |-- Without tidal effects [empty]
|   |   |   |           `-- With tidal effects [empty]
|   |   |   `-- Hyperbolic orbits [empty]
|   |   |       |-- Nonspinning [empty]
|   |   |       `-- Spinning [empty]
|   |   |           |-- Nonprecessing [empty]
|   |   |           |   |-- Without tidal effects [empty]
|   |   |           |   `-- With tidal effects [empty]
|   |   |           `-- Precessing [empty]
|   |   |               |-- Without tidal effects [empty]
|   |   |               `-- With tidal effects [empty]
|   |   `-- Angular momentum flux
|   |       |-- Circular orbits
|   |       |   |-- Nonspinning
|   |       |   `-- Spinning
|   |       |       |-- Nonprecessing
|   |       |       |   |-- With tidal effects
|   |       |       |   `-- Without tidal effects
|   |       |       `-- Precessing [empty]
|   |       |           |-- With tidal effects [empty]
|   |       |           `-- Without tidal effects [empty]
|   |       |-- Hyperbolic orbits [empty]
|   |       |   |-- Spinning [empty]
|   |       |   |   |-- Nonprecessing
|   |       |   |   `-- Precessing
|   |       |   `-- Nonspinning [empty]
|   |       `-- Elliptic orbits
|   |           |-- Spinning [empty]
|   |           |   |-- Nonprecessing
|   |           |   `-- Precessing
|   |           `-- Nonspinning
|   `-- Waveform
|       |-- Circular orbits
|       |   |-- Nonspinning
|       |   `-- Spinning
|       |       |-- Nonprecessing
|       |       |   |-- With tidal effects
|       |       |   `-- Without tidal effects
|       |       `-- Precessing [empty]
|       |           |-- Without tidal effects [empty]
|       |           `-- With tidal effects [empty]
|       |-- Elliptic orbits [empty]
|       |   |-- Nonspinning [empty]
|       |   `-- Spinning [empty]
|       |       |-- Nonprecessing [empty]
|       |       |   |-- Without tidal effects [empty]
|       |       |   `-- With tidal effects [empty]
|       |       `-- Precessing [empty]
|       |           |-- Without tidal effects [empty]
|       |           `-- With tidal effects [empty]
|       `-- Hyperbolic orbits [empty]
|           |-- Nonspinning [empty]
|           `-- Spinning [empty]
|               |-- Nonprecessing [empty]
|               |   |-- Without tidal effects [empty]
|               |   `-- With tidal effects [empty]
|               `-- Precessing [empty]
|                   |-- Without tidal effects [empty]
|                   `-- With tidal effects [empty]
`-- Raw data
    |-- publications
    `-- arXiv
        |-- 0711.0302v2
        |-- 2304.11186v4
        |-- 2304.11185v4
        `-- 2410.16373v1
```
