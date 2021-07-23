rm *.pyf
python -m numpy.f2py interpolie_modeles_alpha.f interpolie_funcs_alpha.f -m interpoliepyc_alpha -h interpoliepyc_alpha.pyf
python -m numpy.f2py -c interpoliepyc_alpha.pyf interpolie_modeles_alpha.f interpolie_funcs_alpha.f
