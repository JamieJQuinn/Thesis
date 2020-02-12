# Anisotropic viscosity with astrophysical applications

## Running analysis

Ensure we have the SDF submodule and all its submodules:
`git submodule update --init --recursive`

Install virtual environment via pipenv:
`cd images`
`pipenv install`

Make SDF's C interface:
`cd SDF/C; make`

Make and link SDF's python interface inside pipenv environment:
`pipenv shell`
`cd SDF/utilities`
`./build -3 -s`

Run notebook:
`cd images`
`pipenv run jupyter notebook`

## TODO
- [ ] Get Athena style wave tests working in lare3d (continue from work done during NAM)
- [ ] Run parameter studies of static kink including anomalous resistivity
- [ ] Run parameter studies of static kink comparing shock with iso and aniso to provide some recommendation
- [ ] Run parameter studies of KHI null point using background & anomalous resist

## Folder Layout

This is the working repo for my thesis, structured in the following way:

- `build/` containing the final latex outputs
- `imagea/` containing the build scripts and outputs for all images
- `chapters/` containing the .tex files for each chapter
- `data-analysis/` containing any miscellaneous analysis not resulting in images
- `makefile` the global makefile for the entire thesis

## Data analysis (not resulting in images)

For analysis outputs like tables that do no result in images, the code will reside in the data-analysis folder within this repo.

## Images

The images will be structured in the following way:

Inside the `images` folder will be a list of projects, roughly corresponding to chapters:
`images/lare1d`, `images/slow_null_point`, `images/kink_instability`, etc

Within the projects involving Lare3d there will be the following:
- Slices produced via python inside one jupyter notebook per chapter-project
`images/project/slices`
- Energy graphs produced via python inside one jupyter notebook per chapter
`images/project/energy_plots`
- 3D field line plots produced via visit inside one folder per simulation
`images/project/3d_field_line_plots`
- Field line integrations produced via Mayavi + python
`images/field_line_integrations`
- Other plots based on the .sdf or en.dat files will be produced in a misc folder with corresponding source code
`images/project/misc`

`images/project/diagrams`

Other outputs for non-lare focussed chapters may have other outputs and folders will be assigned appropriately. Diagrams produced without some source code (e.g. from inkscape) will be held in a single folder per project.

### Guidance on committing

While editing the images, only the source code should be committed. Once the image is finalised in some respect (still to be finalised by supervisors), it should be committed. Since nearly every image should be reproducible, it's not a huge loss to totally remove an image from git's memory if it becomes too large. An exception should be made for those images that are not reproducible, i.e. those produced as diagrams in inkscape.

### On including images

The images may be included with more simplicity using the `graphicspath` command in latex:

```
\graphicspath{{images/}{images/project/}{etc.}}
```

## On intermediate files

The plotting field line integrations requires some intermediate files to be created. Firstly, the raw SDF files must be sliced to reduce file sizes and increase the resolution of the integrator. Secondly, the integrator itself produces a 2D dataset that can be analysed as secondary data. I believe these files should be retained in folders with a structure that matches that of the original data files and all these files stored adjacent to the raw files, probably in a structure that looks like the following:

data/project/type_of_simulation/lare3d_folder/secondary_analysis/ --> sliced_sdfs
                                                                  \-> integrated_quantities

The production of these files will be controlled by the field line integrations images makefile. This allows multiple projects/presentations to access these files, keeps the makefiles simple by using relative paths but keeps the large secondary data files out of the images folder.

## Template

Thesis based on the UofG Science and Engineering Graduate School's [template](https://www.gla.ac.uk/colleges/scienceengineering/graduateschool/postgraduateresearchstudy/submitthesis/).

# Plan

# What needs to be written?

## Introduction

This requires in depth writing after the rest is written. Potentially up to half can be extracted from previous progress reviews, reports and papers.

## Literature Review
    - Role of viscosity in astrophysical MHD
    - Anisotropic nature of coronal viscosity
    - Previous work in anistropic viscosity
    - Previous switching models

This requires heavy editing but I've kept good track of references so it should be a simple affair to gather those. Actually knitting them into a good motivation is another thing.

## The anisotropic tensor
    - Exploration of Braginskii derivation
    - Expansion to switching model
    - Potentially toy problem as pedagogical tool

This is a reimagining of work already published in Braginskii 1965 and MacTaggart et al 2017, so shouldn't be too difficult to rework.

## Experimental setup (2 chapters)
### Chp 1 - Implementation of 1D lare code to describe numerical method

    - Implementation & extrapolation to 3D

### Chp 2 - Implementation of Braginskii & switching models in 3D code

This is mostly done in previous reports. 

    - Testing
    - Comparison with Parrish's similar model in Athena
    - Potential comparison with Waikato's models

## Application to kink instability
    - Motivation
    - Model
    - Analysis is in depth, employing new tools, these will be described

    - TODO Anomolous resistivity

## Application to twisted null point
    - Motivation
    - Model
    - Heating results

The results are written up in MacTaggart, so I'll expand on these, potentially including the work on the Kelvin-Helmholtz in the null point. The KHI work is not complete at this point.

## Application to the KHI in a twisted null point

Work TODO

## Summary & future work suggestions

Not written until most of the rest is done.

## Appendix with tool development
    - KDE plotter
    - SDF splitter
    - field line integrator
    - python analysis scripts

# What needs to be done?

## What simulations & analyses need done?

- Null point Kelvin-Helmholtz runs
  - Parameter study changing visc and resist
  - Want to investigate effect of visc on
    - Reconnection rate
    - Onset time
    - Energy I/O rates

- Kink instability
  - Further work on static
    - Scaling the resultant Ohmic heating?
  - Anomalous resistivity (+ symmetry breaking?)
    - Are the results in any way similar to bg resistivity?
    - What can we say about reconn rate, heating, onset time and relaxation time?
  - Shock viscosity
    - What is the effect of shock vs isotropic?

- Testing of lare3d visc model
  - Similar run in Athena
  - Linear wave simulation

## Worst-case scenario

- KHI can be cut (1 month extra)
- linear waves & comparison to Athena can be cut (1 month extra) 

## Schedule

### Aug 2019 - Jan 2020

- run simulations for 
  - null point KHI
    - parameter study
    - high-res run of interesting case
  - dynamic kink
    - parameter study
    - high-res run of interesting case
  - add anomalous to dynamic and static kink
  - add shock visc to dynamic and static kink

### Feb 2020 - end of internship

- analysis of kink work (4w)

### Mar

- adapt kink paper (1w)
- write up kink work (3w)

### Apr

- analysis of KHI work (4w)

### May

- adapt null paper (1w)
- write up KHI work (3w)

### Jun

- Similar run in Athena (2w)
- linear waves (2w)

### Jul

- implementation in lare3d (2w)
- anisotropic tensor (2w)

### Aug

- lare1d (2w)
- Summary (2w)

### Sept

- tool appendix (4w)

### Oct

- introduction (1w)
- lit review (3w)
