# Thesis plan

Main bulk of writing to occur post-internship, starting Feb 2020. This leaves approx 9 months to complete earlier research on Kelvin-Helmholtz in null points, dynamically twisted coronal loops, or kink instability in a loop using anomolous resistivity.

---

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

### Chp 2 - Implementation of Braginskii & switching models in 3D code

This is mostly done in previous reports. 

    - Testing
    - Comparison with Parrish's similar model in Athena
    - Potential comparison with Waikato's models

## Application to kink instability
    - Motivation
    - Model
    - Analysis is in depth, employing new tools, these will be described

This is all written up for a journal article, so should require only some expansion in the right places.

    - Dynamic twisting of loop into unstable?
    - Anomolous resistivity
    - Shock viscosity

This work is semi-done but needs polishing and proper analysis.

## Application to twisted null point
    - Motivation
    - Model
    - Heating results
    - Kelvin-Helmholtz

The results are written up in MacTaggart, so I'll expand on these, potentially including the work on the Kelvin-Helmholtz in the null point. The KHI work is not complete at this point.

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
  - Dynamically twisted
    - Do the static results generalise?
    - Can we twist the field multiple times?
    - Is the driver in any way realistic?
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
