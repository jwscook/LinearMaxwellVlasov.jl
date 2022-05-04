# LinearMaxwellVlasov.jl

![CI](https://github.com/jwscook/LinearMaxwellVlasov.jl/workflows/CI/badge.svg)
[![codecov.io](http://codecov.io/github/jwscook/LinearMaxwellVlasov.jl/coverage.svg?branch=main)](http://codecov.io/github/jwscook/LinearMaxwellVlasov.jl?branch=main)

Solutions to the homogeneous linear Vlasov-Maxwell equations.

This code solves for the dispersion relation of the linearised Maxwell-Vlasov equations for an infinite, spatially homogenous plasma. Multiple models of plasma species are available: 
 1. cold fluid
 1. warm fluid with optionally distinct parallel and perpendicular sound speeds
 1. kinetic (bi-)Maxwellian with optional parallel drift
 1. kinetic parallel Maxwellian with optional drift with "ring" perpendicular drift
 1. arbitrary decoupled parallel and perpendicular distribution functions. 
 1. relativistic species (not battle tested)

It is possible to solve for complex wavenumbers indicative of convective instabilities.

<img src="/misc/equations/LinearisedMaxwellValasov.png" height="50" />

which is far to big to display here. The relationships between the species contribution and the dielectric tensor and perturbed current are

<img src="/misc/equations/Relationship1.png" height="50" />
<img src="/misc/equations/Relationship2.png" height="50" />

References:

Books: Stix, Melrose, Brambilla

Particularly useful and succinct:
Chapter 15, "Electromagnetic Waves in Plasma" by Takayuki Umeda, Nagoya University Japan in book "Wave Propagation", edited by Andrey Petrin, published March 16th 2011 by IntechOpen
https://www.intechopen.com/books/wave-propagation
