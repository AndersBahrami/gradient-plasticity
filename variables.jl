module Variables

export __NUM_CELLS, __NUM_EQUATIONS, __half_NXC, __half_NYC, __DOMAIN, __MAX_WAVESPEED, __RHO, __YOUNG_MOD_E, __CFL, __DELTA_X, __DELTA_T, __NUM_TIMESTEPS, __modulus, eps_e11, eps_e12, eps_e22, eps_e33, eps11, eps12, eps22, eps33, u1, u2, p, p_til, w, p1, p2, __OUTPATH_TXT, __OUTPATH_VTK, __VTK_NAME, __OUTPATH_GIF, __GIF_NAME, __k, __OUTPATH_STRESS_STRAIN, __OUTPATH_PNG, __PATH, __OUTPATH, __OFFSET, __oneD, __INIT_CONDITIONS, __BOUND_CONDITIONS, __SOFTENING


### Changable Parameters
const __PATH = "output/final_test_sims/"
const __OUTPATH = __PATH*"sim8-3_1D_fast_gradient_softening_noBC_x250y100_t100mod5_ke-5/"

# Initial/Boundary Conditions
# 1: x impact
# 2: y impact
# 3: x pull test
# 4: y pull test
# 5: 
const __INIT_CONDITIONS = 1

# 1: pass through
# 2: reflective
# 3: 
const __BOUND_CONDITIONS = 1
const __SOFTENING = true
const __oneD = false

# Other Parameters
const __NUM_X_CELLS = 100 * 2^n
const __NUM_Y_CELLS = 100 * 2^n
const __X_DOMAIN = 1.0 * 2^n
const __Y_DOMAIN = 1.0 * 2^n
const __NUM_TIMESTEPS = 400 #floor(Int, __NUM_CELLS/5)
const __k = 1e-5
const n = 0
const __modulus = 5 #max(floor(Int, __NUM_TIMESTEPS/30),1) # plots the data every __modulus timesteps
const __OFFSET = 0 # stress strain
const __CFL = .45






# GIF
const __OUTPATH_GIF = __OUTPATH # "images/gifs/" # __SIM_NAME
const __GIF_NAME = "sim.gif"
# PNG
const __OUTPATH_PNG = __OUTPATH*"plot.png"# "images/still_images/plot.png"
# TXT
const __OUTPATH_TXT = __OUTPATH*"txt/" # eps11.txt etc
# VTK
const __OUTPATH_VTK = __OUTPATH*"vtk/" # 0000.vtk
const __VTK_NAME = "variables_"
const __OUTPATH_STRESS_STRAIN = __OUTPATH*"stress_strain.png" # "images/still_images/stress_strain_softening.png"

# Fixed Constants
const eps_e11 = 1
const eps_e12 = 2
const eps_e22 = 3
const eps_e33 = 4
const eps11 = 5
const eps12 = 6
const eps22 = 7
const eps33 = 8
const u1 = 9
const u2 = 10
const p = 11
const p_til = 12
const w = 13
const p1 = 14
const p2 = 15

const __NUM_EQUATIONS = 3 # unused in v4 and above. Maybe in v3. Needed to be 3 for v2/v1.
const __half_NXC = floor(Int, (__NUM_X_CELLS)/2)
const __half_NYC = floor(Int, (__NUM_Y_CELLS)/2)
const __RHO = 7800. # density of iron
const __YOUNG_MOD_E = 2.e11
const __MAX_WAVESPEED = (__YOUNG_MOD_E/__RHO)^0.5 *1.5 # square root of young_mod_e/rho # MAX_WAVESPEED = 6000.
const __DELTA_X = __X_DOMAIN/__NUM_X_CELLS
const __DELTA_Y = __Y_DOMAIN/__NUM_Y_CELLS
const __DELTA_T = __CFL * __DELTA_X / __MAX_WAVESPEED# min(__DELTA_X,__DELTA_Y)/(__MAX_WAVESPEED)



end