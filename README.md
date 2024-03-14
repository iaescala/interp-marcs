# interp-marcs
Python-wrapped, modified version of the MARCS Fortran interpolator (by Thomas Masseron) designed for 4D interpolation/extrapolation including [alpha/Fe]

All modifications in Fortran code to accomodate interpolation of [alpha/Fe] are attributed to I. Escala (Carnegie).

MARCS models (including Fortran code documentation): https://marcs.astro.uu.se/

Much more complete/regularly spaced APOGEE grid of MARCS models (Jonsson+20, by B. Edvardsson), under the "model atmospheres" section (recommended, set grid_type = 'apogee'): https://www.sdss4.org/dr17/irspec/apogee-libraries/

Mar 14 2024 NOTE!!! KNOWN ISSUE INTERPOLATING PLANE-PARALLEL MODELS.

# Usage #

Prior to usage, I recommend reading the manual for the Fortran MARCS interpolator, particularly Sections 4.4 and 4.5.

Example code to generate grid of MARCS model atmospheres contained in testmies.py. The Python function interp_marcs takes the parameters (Teff, Logg, Feh, Alphafe) of the desired interpolated model as input, and finds an appropriate selection of default MARCS models to utilize in generation of the user-specified interpolated model. These selected default models are then fed into the Python-wrapped Fortran interpolator.

Models are interpolated from the default grid with a maximum Fe/H spacing of 0.5 dex and [alpha/Fe] spacing of 0.4 dex (the user can specify these grid spacings). Note that the user must use their own discretion when evaluating the performance of the interpolator outside of the bounds tested by T. Masseron. Also note that this code has been extensively tested for interpolation in Fe/H and Alpha/Fe, but keeping Teff and Logg fixed to non-interpolated grid points. In principle simultaneous extrapolation in 3 or 4 parameters should work, but I have not explicitly tested this.

<!--NOTE!!! It is recommended to install my forked version of Turbospectrum2019 (https://github.com/iaescala/Turbospectrum2019) for use with the output ".mod" files with marcsfile=True. This modified version can read both default and interpolated MARCS models with this flag.-->

NOTE!! Although extrapolation is possible with this code, I do NOT recommend it. Extrapolate at your own risk!!

# Assumptions #

Linear interpolation only. Optimal interpolation (see T. Masseron's documentation) incomptaible with 4D interpolation that includes [alpha/Fe]. Spherical model atmospheres assumed for Logg <= 3.5 unless otherwise specified by user via the optional keyword argument "geometry". For spherical geometry, model mass of 1.0 Msun and microturbulent velocity of 2 km/s is adopted. For plane parallel geometry, micoturbulent velocity of 1 km/s is adopted. However, both model_mass and micro_turb_vel can be explicitly specified by the user. Additionally, by default the code checks whether a file exists in the default MARCS grid (from which the models are interpolated), or in the directory where the user is saving the newly generated models, before proceeding. The code does NOT extrapolate by default, but the user may enable extrpolation using the boolean keyword "extrapol".

# System Requirements #

Depending on the user's operating system, it may be necessary to re-wrap the Fortran code. A bash script (pyrap_alpha.sh) can be run in the same directory as the Fortran files to accomplish this. A Python-wrapped version of the Fortran code comptabible with Linux is provided in this repository. (I have also Python-wrapped the Fortran code with success on macOS.)

# Citation #

For now, if you use this code in your research, please link to this repository, and state that it is I. Escala's modification of T. Masseron's code (Escala et al. 2023, in preparation)


