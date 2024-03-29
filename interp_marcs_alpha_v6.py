"""
@author: Ivanna Escala (Princeton/Carnegie) iescala@carnegiescience.edu

Function: Python wrapper for interpol_modeles.f Fortran routine
(Thomas Masseron) to interpolate Teff, log g, [Fe/H], and [alpha/Fe] from the
MARCS model atmosphere grid (https://marcs.astro.uu.se/).
"""

#from interpoliepyc import interpolie
from interpoliepyc_alpha import interpolie
import numpy as np
import copy
import os

def interp_marcs(teff, logg, feh, alphafe, geometry='sph', model_mass=1.0,
output_model_path='.', input_model_path='.', micro_turb_vel=2.,
check_file_exists=True, optimize=False, extrapol=False, max_step=[0.5, 0.4],
grid_type = 'default'):

  """
  Input parameters:
  ----------------
  teff: float, spectroscopic effective temperature (K)
  logg: float, log surface gravity (cgs)
  feh: float, metallicity (dex)
  alphafe: float, alpha-to-iron ratio (dex)
  geometry: string (optional), either spherical ('sph', valid for logg between
            -1. and 3.5, or giant stars) or plane parallel ('pp',
            valid for logg > 3.0 to 5. or 5.5 for cooler models).
            Default set to spherical for giant stars and plane parallel for
            dwarfs.
  model_mass: float (optional), stellar mass to assume if geometry = 'sph'.
              Default set to the standard model mass of 1.0 Msun.
  micro_turb_vel: float (optional), microturbulent velocity.
                  Default set to 2 km/s (for giant stars).
  output_model_path: string (optional), filepath to output directory for
                     interpolated model
  input_model_path: string (optional), filepath to directory containing
                    default grid of MARCS models for input
  check_file_exists: boolean (optional), check if interpolated model already
                     exists, and if so, do not re-generate
  optimize: boolean (optional), linear or optimal interpolation for the
            Fortran interpolator (see documentation by T. Masseron).
            Default set to linear (optimize = False). Note that optimal
            interpolation will *NOT* work for interpolating in alpha/Fe.
  extrapol: boolean (optional) If True, extrapolate. Default False.
  max_step: 2D array (optional), maximum step in [Fe/H] and [alpha/Fe] for
            separation between grid points for interpolation/extrapolation.
            Default [0.5, 0.4]. See the manual by T. Masseron. Note that
            interpolation error increases for wider step sizes.
  grid_type: string (optional), whether the program is interpolating/extrapolating
             off the default MARCS grid ('default') or the MARCS-APOGEE grid
             ('apogee'). Default is 'default'.

  Output:
  -------
  """

  #Define the MARCS model grid spacing in each dimension
  if grid_type == 'default':
      teff_arr = np.arange(2500., 4000., 100.).tolist() +\
                 np.arange(4000., 8000.+250., 250.).tolist() #THIS SHOULD BE 250
      teff_arr = np.array(teff_arr)
      logg_arr = np.arange(-0.5, 5.5+0.5, 0.5) #THIS SHOULD BE 0.5
      feh_arr = np.array([-5., -4., -3., -2.5, -2., -1.5, -1., -0.75, -0.5, -0.25,
                           0., 0.25, 0.5, 0.75, 1.])
      alphafe_arr = np.array([-0.4, 0., 0.4])

  #Define the MARCS-APOGEE grid spacing in each dimension (except C/M)
  if grid_type == 'apogee':
      teff_arr = np.arange(2800., 4000., 100.).tolist() +\
                 np.arange(4000, 8000.+250., 250.).tolist()
      teff_arr = np.array(teff_arr)
      logg_arr = np.arange(-0.5, 5.5+0.5, 0.5)
      feh_arr = np.round(np.arange(-2.5, 1.+0.25, 0.25), decimals=2)
      alphafe_arr = np.round(np.arange(-1., 1.+0.25, 0.25), decimals=2)

  #If user specified, check if the model exists. If so, skip the interpolation
  #Check if it is either a default model or already generated by user
  if check_file_exists:
      if does_file_exist(teff, logg, feh, alphafe, input_model_path=input_model_path,
      output_model_path=output_model_path, geometry=geometry, grid_type=grid_type,
      micro_turb_vel=micro_turb_vel):
        return

  ####################################################
  #### Select the sixteen nearest atmospheric models ###
  ####################################################

  #Interpolator documentation recommends max 500 K spacing for Teff interpolation,
  #0.5 dex for log g (default), and 0.5 dex for [Fe/H]
  #Use 0.4 dex steps for [alpha/Fe]
  step = np.array([500., 0.5, max_step[0], max_step[1]])
  params = np.array([teff, logg, feh, alphafe])
  params_grid = np.array([teff_arr, logg_arr, feh_arr, alphafe_arr], dtype='object')

  #Check if interpolation point is within bounds of grid
  in_grid = enforce_grid_check(params, params_grid)
  #if user does not desire extrapolation, then abort
  if not in_grid and not extrapol:
      print(f"Input values outside of grid boundaries")
      return

  #Determine the model grid points from which to interpolate
  iparams = select_grid_points(params, params_grid, step)

  ############################################################################
  #### Construct the filenames for the grid of models to be used in the
  #### interpolation
  ############################################################################

  Tefflo, Teffhi = teff_arr[iparams[0]]
  Logglo, Logghi = logg_arr[iparams[1]]
  Fehlo, Fehhi = feh_arr[iparams[2]]
  Alphalo, Alphahi = alphafe_arr[iparams[3]]

  kwargs = (geometry, model_mass, micro_turb_vel, input_model_path, grid_type)

  models = construct_all_model_filenames(Tefflo, Teffhi, Logglo, Logghi,
  Fehlo, Fehhi, Alphalo, Alphahi, kwargs)

  #####################################################################
  ##### Perform sanity checks on whether the grid points exist ########
  #####################################################################

  #Check if each grid point exists
  grid_point_exists = [does_gridpoint_exist(model) for model in models]
  interp_flag = False

  if any(elem is False for elem in grid_point_exists) and in_grid:

      #If only some of grid points don't exist, then
      #try to look for grid points that DO exist
      print("Iterating to find suitable grid point...")

      #If ONLY Fehlo is missing
      if all(elem is False for elem in grid_point_exists[0::2]) and\
         all(elem is True for elem in grid_point_exists[1::2]):

         iparams_new, interp_flag = increment_grid_point(2, iparams, params_grid,
         grid_point_exists, kwargs, kind='-', step=step)

      #If ONLY Fehhi is missing
      if all(elem is False for elem in grid_point_exists[1::2]) and\
         all(elem is True for elem in grid_point_exists[0::2]):

         iparams_new, interp_flag = increment_grid_point(2, iparams, params_grid,
         grid_point_exists, kwargs, kind='+', step=step)

      #If ONLY Alphalo is missing
      if all(elem is False for elem in grid_point_exists[:8]) and\
         all(elem is True for elem in grid_point_exists[8:]):

         iparams_new, interp_flag = increment_grid_point(3, iparams, params_grid,
         grid_point_exists, kwargs, kind='-', step=step)

      #If ONLY Alphahi is missing
      if all(elem is False for elem in grid_point_exists[8:]) and\
         all(elem is True for elem in grid_point_exists[:8]):

          iparams_new, interp_flag = increment_grid_point(3, iparams, params_grid,
          grid_point_exists, kwargs, kind='+', step=step)

      #If ALL grid points are missing
      if all(elem is False for elem in grid_point_exists):

          #Assume that both Fe points are present, but the alpha/Fe points
          #are missing
          if not interp_flag:
              iparams_new, interp_flag = increment_grid_point(3, iparams, params_grid,
              grid_point_exists, kwargs, kind='+-', step=step)

          #Assume that both Fe points are missing, but the alpha/Fe points
          #are present
          if not interp_flag:
              iparams_new, interp_flag = increment_grid_point(2, iparams, params_grid,
              grid_point_exists, kwargs, kind='+-', step=step)

  #Else If all the grid points already exist,
  if all(elem is True for elem in grid_point_exists) and in_grid:
      interp_flag = True
      iparams_new = np.copy(iparams)

  #Check whether the alpha grid points are the same
  #when they should not be, before proceeding with interpolation
  if interp_flag:
      if (Alphalo == Alphahi) and (Alphalo != alphafe):
          #If they are the same in this case, interpolation failed
          interp_flag = False

  ############################################
  ######### Extrapolation ####################
  ############################################

  #EXTRAPOLATE AT YOUR OWN RISK!!!!!

  #extrapolation not allowed, but cannot find necessary grid points
  #abort interpolation
  if (not extrapol) and (not interp_flag):
      print('Interpolation FAILED: Necessary grid points do not exist')
      return

  #extrapolation allowed, and interplation will fail if attempted
  if extrapol and (not interp_flag):

      print('Extrapolating...')

      ######## Try extrapolating in Fe first ############

      iparams_new, interp_flag = increment_grid_point(2, iparams, params_grid,
      grid_point_exists, kwargs, kind='+', extrapol=True, step=step)

      if not interp_flag:
         iparams_new, interp_flag = increment_grid_point(2, iparams, params_grid,
         grid_point_exists, kwargs, kind='-', extrapol=True, step=step)

      #If successful, then make sure that grid points in other dimension
      #are not identical
      if interp_flag:

          Alphalo, Alphahi = alphafe_arr[iparams_new[3]]
          if (Alphalo == Alphahi) and (alphafe != Alphalo):

              iparams_new, interp_flag = increment_grid_point(3, iparams, params_grid,
              grid_point_exists, kwargs, kind='+', extrapol=True, step=step)

              if not interp_flag:
                  iparams_new, interp_flag = increment_grid_point(3, iparams, params_grid,
                  grid_point_exists, kwargs, kind='-', extrapol=True, step=step)

      ########## If extrapolating in Fe doesn't initially work then try #######
      ########## extrapolating in alpha/Fe with Fe constant             #######

      else:

          iparams_new, interp_flag = increment_grid_point(3, iparams, params_grid,
          grid_point_exists, kwargs, kind='+', extrapol=True, step=step)

          if not interp_flag:
              iparams_new, interp_flag = increment_grid_point(3, iparams, params_grid,
              grid_point_exists, kwargs, kind='-', extrapol=True, step=step)

          #If successful, then make sure that the Fe points are not
          #identical when trying to extrapolate in this dimension
          #In practice I don't think this happens, but good to have a check
          if interp_flag:

              Fehlo, Fehhi = feh_arr[iparams_new[2]]
              if (Fehlo == Fehhi) and (Fehlo != feh):

                   iparams_new, interp_flag = increment_grid_point(2, iparams, params_grid,
                   grid_point_exists, kwargs, kind='+', extrapol=True, step=step)

                   if not interp_flag:
                      iparams_new, interp_flag = increment_grid_point(2, iparams, params_grid,
                      grid_point_exists, kwargs, kind='-', extrapol=True, step=step)

      #2D extrapolation time -- actually I think I don't need this??
      if not interp_flag:
         iparams_new, interp_flag = increment_grid_point_2D(iparams, params_grid,
            grid_point_exists, kwargs, step=step)

      #If still failed, then give up and abort
      if not interp_flag:
         print("Extrapolation FAILED: Necessary grid points do not exist")
         return

  ##########################################################
  ##### Final output names and interpolation ###############
  ##########################################################

  Tefflo, Teffhi = teff_arr[iparams_new[0]]
  Logglo, Logghi = logg_arr[iparams_new[1]]
  Fehlo, Fehhi = feh_arr[iparams_new[2]]
  Alphalo, Alphahi = alphafe_arr[iparams_new[3]]

  models = construct_all_model_filenames(Tefflo, Teffhi, Logglo, Logghi,
  Fehlo, Fehhi, Alphalo, Alphahi, kwargs)

  model1, model2, model3, model4, model5, model6, model7, model8, model9,\
  model10, model11, model12, model13, model14, model15, model16 = models

  #Output file names for Turbospec and ATLAS compatible models
  modelout1 = construct_model_filename(teff, logg, feh, alphafe,
              path=output_model_path, geometry=geometry, model_mass=model_mass,
              micro_turb_vel=micro_turb_vel, grid_type=grid_type,
              output=True)
  modelout2 = modelout1[:-4]+'.alt'

  interpolie(model1, model2, model3, model4, model5, model6, model7, model8,
  model9, model10, model11, model12, model13, model14, model15, model16,
  modelout1, modelout2, teff, logg, feh, alphafe, optimize)

  return


def increment_grid_point(index, iparams, params_grid, grid_point_exists,
kwargs, kind='-', extrapol=False, step=None):

    hit_grid_edge = False
    interp_flag = False

    iparams_0 = np.copy(iparams)
    iparams_f = np.copy(iparams)

    while not hit_grid_edge:

        if index == 2:
            low_condition = any(elem is False for elem in grid_point_exists[0::2])
            high_condition = any(elem is False for elem in grid_point_exists[1::2])
        if index == 3:
            low_condition = any(elem is False for elem in grid_point_exists[:8])
            high_condition = any(elem is False for elem in grid_point_exists[8:])

        #If just performing interpolation, then decrease the lower bound
        #and increase the upper bound when looking for existing grid points
        if not extrapol:

            if (kind == '-') or (kind == '+-') or (kind == '-+'):
                if low_condition:
                    iparams_f[index][0] -= 1
                    if iparams_f[index][0] < 0:
                        hit_grid_edge = True
                        continue
                    #if max step size exceed, abort while loop
                    if np.abs(np.diff(params_grid[index][iparams_f[index]])) > step[index]:
                        #iparams_f[index][0] -= 1
                        #if iparams_f[index][0] < 0:
                        #    hit_grid_edge = True
                        #    continue
                        continue

            if (kind == '+') or (kind == '+-') or (kind == '-+'):
                if high_condition:
                    iparams_f[index][1] += 1
                    if iparams_f[index][1] >= params_grid[index].size:
                        hit_grid_edge = True
                        continue
                    #if max step size exceed, abort while loop
                    if np.abs(np.diff(params_grid[index][iparams_f[index]])) > step[index]:
                        #iparams_f[index][1] += 1
                        #if iparams_f[index][1] >= params_grid[index].size:
                        #    hit_grid_edge = True
                        #    continue
                        continue

        #If extrapolation, then simulateneously increase or decrease the
        #upper and lower bounds to find an existing grid point
        if extrapol:

            if kind == '+' and (low_condition or high_condition):

                #Make sure the default selected grid points are not identical
                if np.diff(iparams_f[index]) != 0:
                    if low_condition: iparams_f[index][0] += 1
                    if high_condition: iparams_f[index][1] += 1
                else:
                    if low_condition: iparams_f[index][0] += 1
                    if high_condition: iparams_f[index][1] += 2

                #If the grid edge is hit, then abort while loop
                if iparams_f[index][1] >= params_grid[index].size:
                    hit_grid_edge = True
                    continue

                #if maximum step size is exceeded, abort while loop
                if np.abs(np.diff(params_grid[index][iparams_f[index]]))\
                   >= step[index]:
                    continue
                #else:
                #    iparams_f[index][1] += 1
                #    if iparams_f[index][1] >= params_grid[index].size:
                #        hit_grid_edge = True
                #        continue

            if kind == '-' and (low_condition or high_condition):

                #Make sure the default selected grid points are not identical
                if np.diff(iparams_f[index]) != 0:
                    if low_condition: iparams_f[index][0] -= 1
                    if high_condition: iparams_f[index][1] -= 1
                else:
                    if low_condition: iparams_f[index][0] -= 2
                    if high_condition: iparams_f[index][1] -= 1

                #If the grid edge is hit, then abort while loop
                if iparams_f[index][0] < 0:
                    hit_grid_edge = True
                    continue

                #if maximum step size is exceeded, abort while loop
                if np.abs(np.diff(params_grid[index][iparams_f[index]]))\
                   >= step[index]:
                    continue
                #else:
                #    iparams_f[index][0] -= 1
                #    if iparams_f[index][0] < 0:
                #        hit_grid_edge = True
                #        continue

        Tefflo, Teffhi = params_grid[0][iparams_f[0]]
        Logglo, Logghi = params_grid[1][iparams_f[1]]
        Fehlo, Fehhi = params_grid[2][iparams_f[2]]
        Alphalo, Alphahi = params_grid[3][iparams_f[3]]

        models = construct_all_model_filenames(Tefflo, Teffhi, Logglo, Logghi,
        Fehlo, Fehhi, Alphalo, Alphahi, kwargs)

        grid_point_exists = [does_gridpoint_exist(model) for model in models]
        if all(elem is True for elem in grid_point_exists):
            #Make sure that the grid points are distinct
            #for interpolation/extrapolation
            if np.diff(iparams_f[index]) != 0:
                interp_flag = True
            break

    if not interp_flag:
        iparams_f = iparams_0

    return iparams_f, interp_flag


def construct_all_model_filenames(Tefflo, Teffhi, Logglo, Logghi, Fehlo,
Fehhi, Alphalo, Alphahi, kwargs):

      model1 = construct_model_filename(Tefflo, Logglo, Fehlo, Alphalo, *kwargs)
      model2 = construct_model_filename(Tefflo, Logglo, Fehhi, Alphalo, *kwargs)
      model3 = construct_model_filename(Tefflo, Logghi, Fehlo, Alphalo, *kwargs)
      model4 = construct_model_filename(Tefflo, Logghi, Fehhi, Alphalo, *kwargs)
      model5 = construct_model_filename(Teffhi, Logglo, Fehlo, Alphalo, *kwargs)
      model6 = construct_model_filename(Teffhi, Logglo, Fehhi, Alphalo, *kwargs)
      model7 = construct_model_filename(Teffhi, Logghi, Fehlo, Alphalo, *kwargs)
      model8 = construct_model_filename(Teffhi, Logghi, Fehhi, Alphalo, *kwargs)

      model9 = construct_model_filename(Tefflo, Logglo, Fehlo, Alphahi, *kwargs)
      model10 = construct_model_filename(Tefflo, Logglo, Fehhi, Alphahi, *kwargs)
      model11 = construct_model_filename(Tefflo, Logghi, Fehlo, Alphahi, *kwargs)
      model12 = construct_model_filename(Tefflo, Logghi, Fehhi, Alphahi, *kwargs)
      model13 = construct_model_filename(Teffhi, Logglo, Fehlo, Alphahi, *kwargs)
      model14 = construct_model_filename(Teffhi, Logglo, Fehhi, Alphahi, *kwargs)
      model15 = construct_model_filename(Teffhi, Logghi, Fehlo, Alphahi, *kwargs)
      model16 = construct_model_filename(Teffhi, Logghi, Fehhi, Alphahi, *kwargs)

      models = [model1, model2, model3, model4, model5, model6, model7, model8,
                model9, model10, model11, model12, model13, model14, model15, model16]

      return models

def select_grid_points(params, params_grid, step):
    #Select the appropriate grid points to use for the interpolation
    #Essentially defining a 4D cube

    #Iterate over each parameter
    #define array to store model indices
    iparams = np.zeros((params_grid.size, 2)).astype('int')
    for i in range(params.size):

        #find the "bin" for the parameter
        w = np.digitize(params[i], params_grid[i])

        #If the input parameter value corresponds to a grid point
        #or is beyond the upper bound of the array (e.g. for a/Fe in the default
        #grid)
        if (params[i] in params_grid[i]) or (params[i] > params_grid[i].max()):
            iparam = [w-1, w-1]
        #If the input parameter is beyond the lower bound (e.g. for a/Fe in the
        #default grid)
        elif (params[i] < params_grid[i].min()) :
            iparam = [w, w]
        #If the input parameter is between grid points
        else:
            iparam = [w-1, w]

        iparams[i] = iparam

    return iparams

def construct_model_filename(teffi, loggi, fehi, alphai, geometry='sph',
model_mass=1.0, micro_turb_vel=2., path='.', grid_type='default',
output=False):

        #The user is allowed to define the geometry to some extent,
        #but there are limits as described above
        if loggi < 3.0: geometry = 'sph' #override user input if necessary
        if loggi > 3.5: geometry = 'pp' #only plane parallel models

        if geometry == 'sph':
            geostr = 's'
        if geometry == 'pp':
            geostr = 'p'
            model_mass = 0.0 #overwrite default parameters

        #if its the default MARCS grid, then check the zclass
        if grid_type == 'default':
            try:
                zclass = check_zclass(fehi, alphai)
            except:
                if alphai >= 0.4: zclass = 'ae'
                elif alphai <= -0.4: zclass = 'an'
                else: zclass = 'ap'

        #if using the MARCS-APOGEE grid, then just use 'x3'
        if grid_type == 'apogee':
            zclass = 'x3'

        if np.sign(loggi) == 1 or np.sign(loggi) == 0: gsign = '+'
        else: gsign = '-'

        if np.sign(fehi) == 1 or np.sign(fehi) == 0: fsign = '+'
        else: fsign = '-'

        if np.sign(alphai) == 1 or np.sign(alphai) == 0: asign = '+'
        else: asign = '-'

        feh_str = str(np.abs(fehi))
        if len(feh_str) != 4: feh_str += '0'

        alpha_str = str(np.abs(alphai))
        if len(alpha_str) != 4: alpha_str += '0'

        vt_str = f'0{int(micro_turb_vel)}'

        if grid_type == 'default':
            modelfn = (f"{path}/{geostr}{int(teffi)}_g{gsign}{np.abs(loggi)}"
                       f"_m{model_mass}_t{vt_str}_{zclass}_z{fsign}{feh_str}"
                       f"_a{asign}{alpha_str}_c+0.00_n+0.00_o{asign}{alpha_str}"
                       f"_r+0.00_s+0.00.mod")

        #Some of the models end in 'filled' for MARCS-APOGEE
        if grid_type == 'apogee':

            if not output: path = f"{path}/mod_z{fsign}{feh_str}"

            modelfn = (f"{path}/{geostr}{int(teffi)}_g{gsign}{np.abs(loggi)}"
                       f"_m{model_mass}_t{vt_str}_{zclass}_z{fsign}{feh_str}"
                       f"_a{asign}{alpha_str}_c+0.00_n+0.00_o{asign}{alpha_str}"
                       f"_r+0.00_s+0.00.mod")

            if (not output) and (not does_gridpoint_exist(modelfn)):
                modelfn = modelfn+".filled"

        return modelfn

def enforce_grid_check(params, params_grid):

    #Check that given stellar parameters are in range
    #Do not check for [alpha/Fe]

    #Define the rough bounds of the grid
    bounds = np.array([[params_grid[0].min(), params_grid[0].max()],
             [params_grid[1].min(), params_grid[1].max()],
             [params_grid[2].min(), params_grid[2].max()]])

    #First check that points are within these nominal bounds
    in_grid = True
    if (params[:-1] < bounds[:,0]).any() or (params[:-1] > bounds[:,1]).any():
        in_grid = False

    #Then do a separate check for varying bounds in the grid
    #missing = missing_model(params[0], params[1], params[2])
    #if missing: in_grid = False

    return in_grid

def check_zclass(fehi, alphai):
    #Determine the enrichment class of a set of Fe/H and alpha/Fe
    #If multiple, default to "standard" enrichment class
    if (fehi > -0.75) and (fehi <= 1.) and (alphai == 0.4):
        zclass = 'ae'
    if (fehi >= -2.) and (fehi <= 1.) and (alphai == -0.4):
        zclass = 'an'
    if (fehi >= -2.5) and (fehi <= -0.25) and (alphai == 0.):
        zclass = 'ap'
    if ((fehi >= 0.) and (alphai == 0.)) or \
       ((fehi == -0.25) and (alphai == 0.1)) or \
       ((fehi == -0.5) and (alphai == 0.2)) or \
       ((fehi == -0.75) and (alphai == 0.3)) or \
       ((fehi <= -1.) and (alphai == 0.4)):
        zclass = 'st'
    return zclass

def does_file_exist(teff, logg, feh, alphafe, input_model_path='.',
output_model_path='.', geometry='sph', grid_type='default',
micro_turb_vel=2.):

    file_exists = False

    #for now, check that model is within bounds of alpha/Fe b/c cannot
    #interpolate in this dimension
    model_in_inpath = construct_model_filename(teff, logg, feh, alphafe,
                path=input_model_path, geometry=geometry, micro_turb_vel=micro_turb_vel,
                grid_type=grid_type, output=False)
    model_in_outpath = construct_model_filename(teff, logg, feh, alphafe,
                path=output_model_path, geometry=geometry, micro_turb_vel=micro_turb_vel,
                grid_type=grid_type, output=True)

    if os.path.exists(model_in_inpath):
        print(f"Desired model part of default MARCS model atmosphere set: "
        f"skipping interpolation")
        file_exists = True

    if os.path.exists(model_in_outpath) or os.path.exists(model_in_outpath+'.gz'):
        print("Desired model already generated by user: skipping interpolation")
        file_exists = True

    return file_exists

def does_gridpoint_exist(modelfn):
    grid_point_exists = True
    if not os.path.exists(modelfn):
        grid_point_exists = False
    return grid_point_exists

def increment_grid_point_2D(iparams, params_grid, grid_point_exists_0,
 kwargs, step=None):

    #Increase Fe and alpha
    iparams_f, interp_flag = while_loop(iparams, params_grid, grid_point_exists_0,
    kwargs, kind_f='+', kind_a='+', step=step)

    #Incease Fe and decrease Alpha
    if not interp_flag:
        iparams_f, interp_flag = while_loop(iparams, params_grid, grid_point_exists_0,
         kwargs, kind_f='+', kind_a='-', step=step)

    #Decrease Fe and decrease alpha
    if not interp_flag:
        iparams_f, interp_flag = while_loop(iparams, params_grid, grid_point_exists_0,
        kwargs, kind_f='-', kind_a='-', step=step)

    #Next try decreasing metallicity and increasing [alpha/Fe]
    if not interp_flag:
        iparams_f, interp_flag = while_loop(iparams, params_grid, grid_point_exists_0,
        kwargs, kind_f='-', kind_a='+', step=step)

    return iparams_f, interp_flag

def while_loop(iparams, params_grid, grid_point_exists_0, kwargs,
               kind_f='+', kind_a='+', step=None):

    index_f = 2
    index_a = 3

    iparams_0 = np.copy(iparams)
    iparams_f = np.copy(iparams)

    Tefflo, Teffhi = params_grid[0][iparams_f[0]]
    Logglo, Logghi = params_grid[1][iparams_f[1]]

    is_same_a = np.diff(iparams_f[index_a])[0] == 0
    is_same_f = np.diff(iparams_f[index_f])[0] == 0

    lo_cond_f = any(elem is False for elem in grid_point_exists_0[0::2])|is_same_f
    hi_cond_f = any(elem is False for elem in grid_point_exists_0[1::2])|is_same_f
    lo_cond_a = any(elem is False for elem in grid_point_exists_0[:8])|is_same_a
    hi_cond_a = any(elem is False for elem in grid_point_exists_0[8:])|is_same_a

    interp_flag = False
    grid_edge_cond_f = False
    grid_edge_cond_a = False

    #Step through the Fe grid points, until the edge of the grid is hit
    while not grid_edge_cond_f:

        #Redefine the if_same conditions
        is_same_f = np.diff(iparams_f[index_f])[0] == 0

        #If the two grid points are the same, then first try changing
        #just one of them
        if is_same_f:
            if kind_f == '+':
                iparams_f[index_f][1] += 1
            if kind_f == '-':
                iparams_f[index_f][0] -= 1
        #If not, then change both of the grid points
        else:
            if lo_cond_f:
                    if kind_f == '+':
                        iparams_f[index_f][0] += 1
                    if kind_f == '-':
                        iparams_f[index_f][0] -= 1
            if hi_cond_f:
                if kind_f == '+':
                    iparams_f[index_f][1] += 1
                if kind_f == '-':
                    iparams_f[index_f][1] -= 1

        #Define the grid edge conditions for the extrapolation points
        if kind_f == '+':
            grid_edge_cond_f = iparams_f[index_f][1] >= params_grid[index_f].size
        if kind_f == '-':
            grid_edge_cond_f = iparams_f[index_f][0] < 0

        #If the edge of the grid is hit when searching for an extrapolation
        #point, then abort
        if grid_edge_cond_f:
            return iparams_0, interp_flag
        #If the alpha grid edge was hit in the last iteration of the alpha
        #while loop, but the Fe points exist, then select new Fe grid points
        #then re-evaluate the grid edge conditions
        else:
            if grid_edge_cond_a and (not lo_cond_f and not hi_cond_f):
                if kind_f == '+':
                    iparams_f[index_f] += 1
                    grid_edge_cond_f = iparams_f[index_f][1] >= params_grid[index_f].size
                if kind_f == '-':
                    iparams_f[index_f] -= 1
                    grid_edge_cond_f = iparams_f[index_f][0] < 0
                if grid_edge_cond_f:
                    return iparams_0, interp_flag

        #Ensure that the maximum grid step requirement is met when searching
        #for extrapolation points

        #while np.abs(np.diff(params_grid[index_f][iparams_f[index_f]])) < min_step[index_f]:
        #    if kind_f == '+':
        #        iparams_f[index_f][1] += 1
        #        grid_edge_cond_f = iparams_f[index_f][1] >= params_grid[index_f].size
        #    if kind_f == '-':
        #        iparams_f[index_f][0] -= 1
        #        grid_edge_cond_f = iparams_f[index_f][0] < 0
        #    if grid_edge_cond_f:
        #        return iparams_0, interp_flag

        if np.abs(np.diff(params_grid[index_f][iparams_f[index_f]])) > step[index_f]:
            return iparams_0, interp_flag

        #Determine the new extrapolation points
        Fehlo, Fehhi = params_grid[index_f][iparams_f[index_f]]

        grid_edge_cond_a = False
        #Step through alpha/Fe until the grid edge is hit (for fixed Fe points)
        while not grid_edge_cond_a:

            #Redefine the if_same conditions
            is_same_a = np.diff(iparams_f[index_a]) == 0

            #Do the same thing for the alpha points
            if is_same_a:
                if kind_a == '+':
                    iparams_f[index_a][1] += 1
                if kind_a == '-':
                    iparams_f[index_a][0] -= 1
            else:
                if lo_cond_a:
                    if kind_a == '+':
                        iparams_f[index_a][0] += 1
                    if kind_a == '-':
                        iparams_f[index_a][0] -= 1
                if hi_cond_a:
                    if kind_a == '+':
                        iparams_f[index_a][1] += 1
                    if kind_a == '-':
                        iparams_f[index_a][1] -= 1

            if kind_a == '+':
                grid_edge_cond_a = iparams_f[index_a][1] >= params_grid[index_a].size
            if kind_a == '-':
                grid_edge_cond_a = iparams_f[index_a][0] < 0

            #If the edge of the grid is hit, then reset the values of
            #the alpha grid points and exit the while loop
            if grid_edge_cond_a:
                iparams_f[index_a] = iparams_0[index_a]
                break

            #Do the same thing for alpha

            #while np.abs(np.diff(params_grid[index_a][iparams_f[index_a]])) < min_step[index_a]:
            #    if kind_a == '+':
            #        iparams_f[index_a][1] += 1
            #        grid_edge_cond_a = iparams_f[index_a][1] >= params_grid[index_a].size
            #    if kind_a == '-':
            #        iparams_f[index_a][0] -= 1
            #        grid_edge_cond_a = iparams_f[index_a][0] < 0
            #    if grid_edge_cond_a:
            #        iparams_f[index_a] = iparams_0[index_a]
            #        break

            if np.abs(np.diff(params_grid[index_a][iparams_f[index_a]])) > step[index_a]:
                break

            #Need an additional break statement, since the previous
            #one only exits the first while loop
            if grid_edge_cond_a:
                break

            #Determine the new alpha grid points
            Alphalo, Alphahi = params_grid[index_a][iparams_f[index_a]]

            #Construct the associated models
            models = construct_all_model_filenames(Tefflo, Teffhi, Logglo, Logghi,
            Fehlo, Fehhi, Alphalo, Alphahi, kwargs)

            #Evaluate whether these models exist, and if not, repeat loop
            grid_point_exists = [does_gridpoint_exist(model) for model in models]

            #Redefine the high and low conditions
            lo_cond_a = any(elem is False for elem in grid_point_exists[:8])
            hi_cond_a = any(elem is False for elem in grid_point_exists[8:])

            lo_cond_f = any(elem is False for elem in grid_point_exists[0::2])
            hi_cond_f = any(elem is False for elem in grid_point_exists[1::2])

            #If all the grid points exist, then proceed
            if all(elem is True for elem in grid_point_exists):
                interp_flag = True
                break

        if interp_flag:
            break

    return iparams_f, interp_flag
