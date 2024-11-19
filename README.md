# GXBeamCS

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://byuflowlab.github.io/GXBeamCS.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://byuflowlab.github.io/GXBeamCS.jl/dev/)
[![Build Status](https://github.com/byuflowlab/GXBeamCS.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/byuflowlab/GXBeamCS.jl/actions/workflows/CI.yml?query=branch%3Amain)

*Composite beam cross sectional analysis in Julia*

Authors: Andrew Ning and Adam Cardoza

While this package can stand alone, it is designed as part of the [GXBeam](https://github.com/byuflowlab/GXBeam.jl) ecosystem.  This package computes compliance and inertial properties for composite beam cross sections, which are inputs to GXBeam.  After running GXBeam to compute forces and moments, GXBeamCS can map those loads back to strains, stresses, and Tsai-Wu failure criteria on the cross sections.

## Package Features

We provide alternative approaches to computing section properties depending on the desired level of fidelity.  All methods compute compliance matrices, mass matrices for pre-beam analysis, and compute strains, stresses, and Tsai-Wu failure criteria for post-beam analysis.  Each method also has some basic visualization capabilities.  Our methods have been written to be compatible with algorithmic differentiation.

### Finite Element Analysis

This approach has much greater accuracy than the classical laminate theory based methods, particularly for thicker structures and/or when greater anisotropy exists in the laminates.  Features of the methodology include:
- Quadrilateral finite elements (mesh generation discussed below) for better accuracy and cross coupling.
- General geometry definition with inhomogenous properties and anisotropic behavior (computes full 6x6 compliance matrix)
- Ply materials are general orthotropic (not plane orthotropic as in CLT)
- Strain recovery back to mesh

For mesh generation we provide several options.

- A built in mesh tool for airfoils.
- (not yet implemented) a wrapper to Airfoil2Becas (requires obtaining a separate license to Becas)
- (not yet implemented) a wrapper to PreVABS

The PreVABS mesher is the best of these.  It is freely available, but it requires compiling C++ code, and as such is not differentiable for optimization studies.  Airfoil2Becas is a simpler mesher, not quite as accurate, but is reasoanbly robust.  It is written in Python so will also not be differentiable for optimization studies.  The built-in mesher is written in pure Julia and is AD-compatible.  It however is not nearly as robust at this time.

In future work we intend to integrate a free form deformation approach so that the cross-sectional mesh can be changed in an optimization application without mesh regeneration.  This will enable any of the above approaches to be used while providing derivatives.

### Classical Laminate Theory (work in progress)

Classical laminate theory is a simpler alternative.  It assumes that the laminates are very thin compared to its other dimensions so that a plane stress assumption can be justified.  It is less accurate than the finite element based approach, but is faster and potentially more robust and so may be appropriate depending on the level of fidelity needed.  It is a reasonable approach when the laminates are thin and when the primary stiffness terms of interest are axial, bending, and torsion (it is typically much less accurate in shear and in cross-coupled terms).

We provide two options.  First is a Julia implementation of [PreComp](https://www.nrel.gov/wind/nwtc/precomp.html), which is a Fortran code developed by NREL.  PreComp only provides the pre-beam analysis but not post-beam (no mapping back to strains, stresses, and failure criteria).

The second is an in-house implementation of CLT.  An alternative CLT was pursued because of the above limitations (and limited documentation to easily extend the existing code).  The CLT theory for laminates is standard (theory derived in a supplemental document), and the CLT for beams is based on the method outlined by Koller and Springer [ref].

<!-- ## Installation

Enter the package manager by typing `]` and then run the following:

```julia
pkg> add GXBeamCS
```
-->

## Usage

See the [documentation](https://flow.byu.edu/GXBeamCS.jl/stable)

## Related Codes

The finite element cross sectional analysis uses the same underlying theory as in [BECAS](https://becas.dtu.dk), but was written independently to be fast and optimization-friendly.  [PreComp](https://www.nrel.gov/wind/nwtc/precomp.html) is an open source classical lamiante theory implementation and was rewritten in Julia and wrapped as an option in this package.  [VABS](https://analyswift.com/vabs-cross-sectional-analysis-tool-for-composite-beams/) is another popular tool, and is the best of these in terms of analysis accuracy, though is not freely available.

