# interp-marcs
Python-wrapped, modified version of the MARCS Fortran interpolator (by Thomas Masseron) designed for 4D interpolation/extrapolation including [alpha/Fe]

All modifications in Fortran code to accomodate interpolation/extrapolation of [alpha/Fe] are attributed to I. Escala (Carnegie).

MARCS models (including Fortran code documentation): https://marcs.astro.uu.se/

# Usage #

Prior to usage, I recommend reaading the manual for the Fortran MARCS interpolator, particularly Sections 4.4 and 4.5.

Example code to generate grid of MARCS model atmospheres contained in testmies.py. The Python function interp_marcs takes the parameters (Teff, Logg, Feh, Alphafe) of the desired interpolated/extrapolated model as input, and finds an appropriate selection of default MARCS models to utilize in generation of the user-specified interpolated/extrapolated model. These selected default models are then fed into the Python-wrapped Fortran interpolator.

Models are interpolated/extrapolated from the default grid with a minimum Fe/H spacing of 0.5 dex and [alpha/Fe] spacing of 0.4 dex. Note that the user must use their own discretion when evaluating the performance of the interpolator/extrapolator outside of the bounds tested by T. Masseron. Also note that this code has been extensively tested for interpolation/extrapolation simultaneously in Fe/H and Alpha/Fe, but keeping Teff and Logg fixed to default MARCS grid points. In principle simultaneous extrapolation in 3 or 4 parameters should work, but I have not explicitly tested this.

<!--NOTE!!! It is recommended to install my forked version of Turbospectrum2019 (https://github.com/iaescala/Turbospectrum2019) for use with the output ".mod" files with marcsfile=True. This modified version can read both default and interpolated MARCS models with this flag.-->

# Assumptions #

Linear interpolation only. Optimal interpolation (see T. Masseron's documentation) incomptaible with 4D interpolation that includes [alpha/Fe]. Spherical model atmospheres assumed for Logg <= 3.5 unless otherwise specified by user via the optional keyword argument "geometry". For spherical geometry, model mass of 1.0 Msun and microturbulent velocity of 2 km/s is adopted. For plane parallel geometry, micoturbulent velocity of 1 km/s is adopted. However, both model_mass and micro_turb_vel can be explicitly specified by the user. Additionally, by default the code checks whether a file exists in the default MARCS grid, or in the directory where the user is saving the newly generated models, before proceeding. The code extrapolates by default, but the user may restrict the program to interpolation only.

# System Requirements #

Depending on the user's operating system, it may be necessary to re-wrap the Fortran code. A bash script (pyrap_alpha.sh) can be run in the same directory as the Fortran files to accomplish this. A Python-wrapped version of the Fortran code comptabible with Linux is provided in this repository. (I have also Python-wrapped the Fortran code with success on macOS.)

# Citation #

For now, if you use this code in your research, please link to this repository, and state that it is I. Escala's modification of T. Masseron's code


