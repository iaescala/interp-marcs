from interp_marcs_alpha_v6 import interp_marcs
import numpy as np

input_model_path='/home/iescala/data1/atmospheres'
output_model_path='GridLinSphPP4'

teff_arr = [3200,3300,3400,3500,3600,3700,3800,3900,4000,4250,4500,4750,
            5000]
logg_arr = np.arange(0., 5.5, 0.5)

feh_arr = np.arange(-4., 1.5, 0.5)
alphafe_arr = [-0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.,
                0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]

for teff in teff_arr:
	for logg in logg_arr:
		for feh in feh_arr:
			for alphafe in alphafe_arr:

				print(teff, logg, feh, alphafe)

				interp_marcs(teff, logg, feh, alphafe,
				output_model_path=output_model_path,
				input_model_path=input_model_path,
				check_file_exists=True, extrapol=False,
				geometry='sph')
