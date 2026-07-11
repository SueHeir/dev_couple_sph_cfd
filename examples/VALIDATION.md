# Validation index

## Coupled SPH-CFD plume/surface minimum fluidization

[The plume_surface validation](plume_surface/README.md) runs the live SPH-CFD
coupling at minimum fluidization. Its committed
[result figure](plume_surface/plots/plume_surface_validation.svg) compares the
measured SPH `U_mf` with the independent Wen-Yu correlation and a resolved
DEM-CFD cross-method reference, with their tolerance bands shown. It also
shows two executable fault controls (omitted pressure gradient and an
incorrect voidage exponent) outside those bands, plus the dynamic
pressure/onset gate above and below `U_mf`.

The check is a finite, coarse regression case rather than experimental proof
of plume-surface interaction; its reference and model limitations are stated
with the case configuration and source output.
