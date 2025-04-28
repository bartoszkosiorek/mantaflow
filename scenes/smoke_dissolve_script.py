
######################################################################
## LIBRARIES
######################################################################
from manta import *
from math import inf
import os.path, shutil, math, sys, gc, multiprocessing, platform, time

withMPBake = False # Bake files asynchronously
withMPSave = False # Save files asynchronously
isWindows = platform.system() != 'Darwin' and platform.system() != 'Linux'
# TODO(sebbas): Use this to simulate Windows multiprocessing (has default mode spawn)
#try:
#    multiprocessing.set_start_method('spawn')
#except:
#    pass

######################################################################
## VARIABLES
######################################################################

mantaMsg('Fluid variables')
dim_s7     = 3
res_s7     = 32
gravity_s7 = vec3(0.000000, 0.000000, -9.810000) # in SI unit (e.g. m/s^2)
gs_s7      = vec3(21, 21, 32)
maxVel_s7  = 0

domainClosed_s7     = False
boundConditions_s7  = 'xXyYz'
boundaryWidth_s7    = 1
deleteInObstacle_s7 = False

using_smoke_s7        = True
using_liquid_s7       = False
using_noise_s7        = True
using_adaptTime_s7    = True
using_obstacle_s7     = False
using_guiding_s7      = False
using_fractions_s7    = False
using_invel_s7        = False
using_outflow_s7      = False
using_sndparts_s7     = False
using_speedvectors_s7 = False
using_diffusion_s7    = False

# Fluid time params
timeScale_s7    = 1.000000
timeTotal_s7    = 0.000000
timePerFrame_s7 = 0.000000

# In Blender fluid.c: frame_length = DT_DEFAULT * (25.0 / fps) * time_scale
# with DT_DEFAULT = 0.1
frameLength_s7 = 0.104167
frameLengthUnscaled_s7 = frameLength_s7 / timeScale_s7
frameLengthRaw_s7 = 0.1 * 25 # dt = 0.1 at 25 fps

dt0_s7          = 0.104167
cflCond_s7      = 4.000000
timestepsMin_s7 = 1
timestepsMax_s7 = 4

# Start and stop for simulation
current_frame_s7 = 1
start_frame_s7   = 1
end_frame_s7     = 250

# Fluid diffusion / viscosity
domainSize_s7 = 6.000000 # longest domain side in meters
kinViscosity_s7 = 0.000001 / (domainSize_s7*domainSize_s7) # kinematic viscosity in m^2/s

# Factors to convert Blender units to Manta units
ratioMetersToRes_s7 = float(domainSize_s7) / float(res_s7) # [meters / cells]
mantaMsg('1 Mantaflow cell is ' + str(ratioMetersToRes_s7) + ' Blender length units long.')

ratioResToBLength_s7 = float(res_s7) / float(domainSize_s7) # [cells / blength] (blength: cm, m, or km, ... )
mantaMsg('1 Blender length unit is ' + str(ratioResToBLength_s7) + ' Mantaflow cells long.')

ratioBTimeToTimestep_s7 = float(1) / float(frameLengthRaw_s7) # the time within 1 blender time unit, see also fluid.c
mantaMsg('1 Blender time unit is ' + str(ratioBTimeToTimestep_s7) + ' Mantaflow time units long.')

ratioFrameToFramelength_s7 = float(1) / float(frameLengthUnscaled_s7 ) # the time within 1 frame
mantaMsg('frame / frameLength is ' + str(ratioFrameToFramelength_s7) + ' Mantaflow time units long.')

scaleAcceleration_s7 = ratioResToBLength_s7 * (ratioBTimeToTimestep_s7**2)# [meters/btime^2] to [cells/timestep^2] (btime: sec, min, or h, ...)
mantaMsg('scaleAcceleration is ' + str(scaleAcceleration_s7))

scaleSpeedFrames_s7 = ratioResToBLength_s7 * ratioFrameToFramelength_s7 # [blength/frame] to [cells/frameLength]
mantaMsg('scaleSpeed is ' + str(scaleSpeedFrames_s7))

gravity_s7 *= scaleAcceleration_s7 # scale from world acceleration to cell based acceleration

# OpenVDB options
vdbCompression_s7 = Compression_Blosc
# vdbPrecision_s7 = Precision_Half
vdbClip_s7 = 0.000001

# Cache file names
file_data_s7 = 'fluid_data'
file_noise_s7 = 'fluid_noise'
file_mesh_s7 = 'fluid_mesh'
file_meshvel_s7 = 'fluid_mesh'
file_particles_s7 = 'fluid_particles'
file_guiding_s7 = 'fluid_guiding'
mantaMsg('Smoke variables low')
preconditioner_s7    = PcMGStatic
using_colors_s7      = True
using_heat_s7        = True
using_fire_s7        = True
using_noise_s7       = True
vorticity_s7         = 0.000000
buoyancy_dens_s7     = float(1.000000) / float(6.000000)
buoyancy_heat_s7     = float(1.000000) / float(6.000000)
dissolveSpeed_s7     = 5
using_logdissolve_s7 = True
using_dissolve_s7    = False
flameVorticity_s7    = 0.500000
burningRate_s7       = 0.750000
flameSmoke_s7        = 1.000000
ignitionTemp_s7      = 1.500000
maxTemp_s7           = 3.000000
flameSmokeColor_s7   = vec3(0.700000,0.700000,0.700000)

mantaMsg('Fluid variables noise')
upres_sn7  = 2
gs_sn7     = vec3(upres_sn7*gs_s7.x, upres_sn7*gs_s7.y, upres_sn7*gs_s7.z)

mantaMsg('Smoke variables noise')
wltStrength_s7 = 1.000000
uvs_s7         = 2
uvs_offset_s7  = vec3(0, 0, 0)
octaves_s7     = int(math.log(upres_sn7) / math.log(2.0) + 0.5) if (upres_sn7 > 1) else 1

######################################################################
## SOLVERS
######################################################################

mantaMsg('Solver base')
s7 = Solver(name='solver_base7', gridSize=gs_s7, dim=dim_s7)

mantaMsg('Solver noise')
sn7 = Solver(name='solver_noise7', gridSize=gs_sn7)

######################################################################
## GRIDS
######################################################################

mantaMsg('Fluid alloc data')
flags_s7       = s7.create(FlagGrid, name='flags')
vel_s7         = s7.create(MACGrid, name='velocity', sparse=True)
velTmp_s7      = s7.create(MACGrid, name='velocity_previous', sparse=True)
x_vel_s7       = s7.create(RealGrid, name='x_vel')
y_vel_s7       = s7.create(RealGrid, name='y_vel')
z_vel_s7       = s7.create(RealGrid, name='z_vel')
pressure_s7    = s7.create(RealGrid, name='pressure')
phiObs_s7      = s7.create(LevelsetGrid, name='phi_obstacle')
phiSIn_s7      = s7.create(LevelsetGrid, name='phiSIn') # helper for static flow objects
phiIn_s7       = s7.create(LevelsetGrid, name='phi_inflow')
phiOut_s7      = s7.create(LevelsetGrid, name='phi_out')
forces_s7      = s7.create(Vec3Grid, name='forces')
x_force_s7     = s7.create(RealGrid, name='x_force')
y_force_s7     = s7.create(RealGrid, name='y_force')
z_force_s7     = s7.create(RealGrid, name='z_force')
obvel_s7       = None

# Set some initial values
phiObs_s7.setConst(9999)
phiSIn_s7.setConst(9999)
phiIn_s7.setConst(9999)
phiOut_s7.setConst(9999)

# Keep track of important objects in dict to load them later on
fluid_data_dict_final_s7  = { 'vel' : vel_s7 }
fluid_data_dict_resume_s7 = { 'phiObs' : phiObs_s7, 'phiIn' : phiIn_s7, 'phiOut' : phiOut_s7, 'flags' : flags_s7, 'velTmp' : velTmp_s7 }

mantaMsg('Smoke alloc')
shadow_s7     = s7.create(RealGrid, name='shadow', sparse=False)
emission_s7   = s7.create(RealGrid, name='emission', sparse=True)
emissionIn_s7 = s7.create(RealGrid, name='emissionIn')
density_s7    = s7.create(RealGrid, name='density', sparse=True)
densityIn_s7  = s7.create(RealGrid, name='density_inflow', sparse=True)
heat_s7       = None # allocated dynamically
heatIn_s7     = None
flame_s7      = None
fuel_s7       = None
react_s7      = None
fuelIn_s7     = None
reactIn_s7    = None
color_r_s7    = None
color_g_s7    = None
color_b_s7    = None
color_r_in_s7 = None
color_g_in_s7 = None
color_b_in_s7 = None

# Set some initial values
shadow_s7.setConst(-1)

# Keep track of important objects in dict to load them later on
smoke_data_dict_final_s7 = { 'density' : density_s7, 'shadow' : shadow_s7 }
smoke_data_dict_resume_s7 = { 'densityIn' : densityIn_s7, 'emission' : emission_s7 }

mantaMsg('Smoke alloc noise')
vel_sn7        = sn7.create(MACGrid, name='velocity_noise')
density_sn7    = sn7.create(RealGrid, name='density_noise', sparse=True)
phiIn_sn7      = sn7.create(LevelsetGrid, name='phiIn_noise')
phiOut_sn7     = sn7.create(LevelsetGrid, name='phiOut_noise')
phiObs_sn7     = sn7.create(LevelsetGrid, name='phiObs_noise')
flags_sn7      = sn7.create(FlagGrid, name='flags_noise')
tmpIn_sn7      = sn7.create(RealGrid, name='tmpIn_noise')
emissionIn_sn7 = sn7.create(RealGrid, name='emissionIn_noise')
energy_s7      = s7.create(RealGrid, name='energy')
tmpFlags_s7    = s7.create(FlagGrid, name='tmpFlags')
texture_u_s7   = s7.create(RealGrid, name='textureU')
texture_v_s7   = s7.create(RealGrid, name='textureV')
texture_w_s7   = s7.create(RealGrid, name='textureW')
texture_u2_s7  = s7.create(RealGrid, name='textureU2')
texture_v2_s7  = s7.create(RealGrid, name='textureV2')
texture_w2_s7  = s7.create(RealGrid, name='textureW2')
flame_sn7      = None
fuel_sn7       = None
react_sn7      = None
color_r_sn7    = None
color_g_sn7    = None
color_b_sn7    = None
wltnoise_sn7   = sn7.create(NoiseField, fixedSeed=265, loadFromFile=True)

mantaMsg('Initializing UV Grids')
uvGrid0_s7 = s7.create(VecGrid, name='uv_grid_0', sparse=False)
uvGrid1_s7 = s7.create(VecGrid, name='uv_grid_1', sparse=False)
resetUvGrid(target=uvGrid0_s7, offset=uvs_offset_s7)
resetUvGrid(target=uvGrid1_s7, offset=uvs_offset_s7)

# Sync UV and texture grids
copyVec3ToReal(source=uvGrid0_s7, targetX=texture_u_s7, targetY=texture_v_s7, targetZ=texture_w_s7)
copyVec3ToReal(source=uvGrid1_s7, targetX=texture_u2_s7, targetY=texture_v2_s7, targetZ=texture_w2_s7)

# Keep track of important objects in dict to load them later on
smoke_noise_dict_final_s7 = { 'density_noise' : density_sn7 }
smoke_noise_dict_resume_s7 = { 'uv0_noise' : uvGrid0_s7, 'uv1_noise' : uvGrid1_s7 }

# Sanity check, clear grids first
if 'color_r_sn7' in globals(): del color_r_sn7
if 'color_g_sn7' in globals(): del color_g_sn7
if 'color_b_sn7' in globals(): del color_b_sn7

mantaMsg('Allocating colors noise')
color_r_sn7 = sn7.create(RealGrid, name='color_r_noise', sparse=True)
color_g_sn7 = sn7.create(RealGrid, name='color_g_noise', sparse=True)
color_b_sn7 = sn7.create(RealGrid, name='color_b_noise', sparse=True)

# Add objects to dict to load them later on
if 'smoke_noise_dict_final_s7' in globals():
    smoke_noise_dict_final_s7.update(color_r_noise=color_r_sn7, color_g_noise=color_g_sn7, color_b_noise=color_b_sn7)

# Sanity check, clear grids first
if 'flame_sn7' in globals(): del flame_sn7
if 'fuel_sn7' in globals(): del fuel_sn7
if 'react_sn7' in globals(): del react_sn7

mantaMsg('Allocating fire noise')
flame_sn7 = sn7.create(RealGrid, name='flame_noise', sparse=True)
fuel_sn7  = sn7.create(RealGrid, name='fuel_noise', sparse=True)
react_sn7 = sn7.create(RealGrid, name='react_noise', sparse=True)

# Add objects to dict to load them later on
if 'smoke_noise_dict_final_s7' in globals():
    smoke_noise_dict_final_s7.update(flame_noise=flame_sn7)
if 'smoke_noise_dict_resume_s7' in globals():
    smoke_noise_dict_resume_s7.update(fuel_noise=fuel_sn7, react_noise=react_sn7)

# Sanity check, clear grids first
if 'heat_s7' in globals(): del heat_s7
if 'heatIn_s7' in globals(): del heatIn_s7

mantaMsg('Allocating heat')
heat_s7   = s7.create(RealGrid, name='temperature', sparse=True)
heatIn_s7 = s7.create(RealGrid, name='temperature_inflow', sparse=True)

# Add objects to dict to load them later on
if 'smoke_data_dict_final_s7' in globals():
    smoke_data_dict_final_s7.update(heat=heat_s7)
if 'smoke_data_dict_resume_s7' in globals():
    smoke_data_dict_resume_s7.update(heatIn=heatIn_s7)

# Sanity check, clear grids first
if 'color_r_s7' in globals(): del color_r_s7
if 'color_g_s7' in globals(): del color_g_s7
if 'color_b_s7' in globals(): del color_b_s7

mantaMsg('Allocating colors')
color_r_s7    = s7.create(RealGrid, name='color_r', sparse=True)
color_g_s7    = s7.create(RealGrid, name='color_g', sparse=True)
color_b_s7    = s7.create(RealGrid, name='color_b', sparse=True)
color_r_in_s7 = s7.create(RealGrid, name='color_r_inflow', sparse=True)
color_g_in_s7 = s7.create(RealGrid, name='color_g_inflow', sparse=True)
color_b_in_s7 = s7.create(RealGrid, name='color_b_inflow', sparse=True)

# Add objects to dict to load them later on
if 'smoke_data_dict_final_s7' in globals():
    smoke_data_dict_final_s7.update(color_r=color_r_s7, color_g=color_g_s7, color_b=color_b_s7)
if 'smoke_data_dict_resume_s7' in globals():
    smoke_data_dict_resume_s7.update(color_r_in=color_r_in_s7, color_g_in=color_g_in_s7, color_b_in=color_b_in_s7)

# Sanity check, clear grids first
if 'flame_s7' in globals(): del flame_s7
if 'fuel_s7' in globals(): del fuel_s7
if 'react_s7' in globals(): del react_s7
if 'fuelIn_s7' in globals(): del fuelIn_s7
if 'reactIn_s7' in globals(): del reactIn_s7

mantaMsg('Allocating fire')
flame_s7   = s7.create(RealGrid, name='flame', sparse=True)
fuel_s7    = s7.create(RealGrid, name='fuel', sparse=True)
react_s7   = s7.create(RealGrid, name='react', sparse=True)
fuelIn_s7  = s7.create(RealGrid, name='fuel_inflow', sparse=True)
reactIn_s7 = s7.create(RealGrid, name='react_inflow', sparse=True)

# Add objects to dict to load them later on
if 'smoke_data_dict_final_s7' in globals():
    smoke_data_dict_final_s7.update(flame=flame_s7)
if 'smoke_data_dict_resume_s7' in globals():
    smoke_data_dict_resume_s7.update(fuel=fuel_s7, react=react_s7, fuelIn=fuelIn_s7, reactIn=reactIn_s7)

# wavelet noise params
wltnoise_sn7.posScale = vec3(int(21), int(21), int(32)) * (1. / 2.000000)
wltnoise_sn7.timeAnim = 0.100000

######################################################################
## ADAPTIVE TIME
######################################################################

mantaMsg('Fluid adaptive time stepping')
s7.frameLength  = frameLength_s7
s7.timestepMin  = s7.frameLength / max(1, timestepsMax_s7)
s7.timestepMax  = s7.frameLength / max(1, timestepsMin_s7)
s7.cfl          = cflCond_s7
s7.timePerFrame = timePerFrame_s7
s7.timestep     = dt0_s7
s7.timeTotal    = timeTotal_s7
#mantaMsg('timestep: ' + str(s7.timestep) + ' // timPerFrame: ' + str(s7.timePerFrame) + ' // frameLength: ' + str(s7.frameLength) + ' // timeTotal: ' + str(s7.timeTotal) )

def fluid_adapt_time_step_7():
    mantaMsg('Fluid adapt time step')
    
    # time params are animatable
    s7.frameLength = frameLength_s7
    s7.cfl         = cflCond_s7
    s7.timestepMin  = s7.frameLength / max(1, timestepsMax_s7)
    s7.timestepMax  = s7.frameLength / max(1, timestepsMin_s7)
    
    # ensure that vel grid is full (remember: adaptive domain can reallocate solver)
    copyRealToVec3(sourceX=x_vel_s7, sourceY=y_vel_s7, sourceZ=z_vel_s7, target=vel_s7)
    maxVel_s7 = vel_s7.getMax() if vel_s7 else 0
    if using_adaptTime_s7:
        mantaMsg('Adapt timestep, maxvel: ' + str(maxVel_s7))
        s7.adaptTimestep(maxVel_s7)

######################################################################
## IMPORT
######################################################################

def fluid_file_import_s7(dict, path, framenr, file_format, file_name=None):
    mantaMsg('Fluid file import, frame: ' + str(framenr))
    try:
        framenr = fluid_cache_get_framenr_formatted_7(framenr)
        # New cache: Try to load the data from a single file
        loadCombined = 0
        if file_name is not None:
            file = os.path.join(path, file_name + '_' + framenr + file_format)
            if os.path.isfile(file):
                if file_format == '.vdb':
                    loadCombined = load(name=file, objects=list(dict.values()), worldSize=domainSize_s7)
                elif file_format == '.bobj.gz' or file_format == '.obj':
                    for name, object in dict.items():
                        if os.path.isfile(file):
                            loadCombined = object.load(file)
        
        # Old cache: Try to load the data from separate files, i.e. per object with the object based load() function
        if not loadCombined:
            for name, object in dict.items():
                file = os.path.join(path, name + '_' + framenr + file_format)
                if os.path.isfile(file):
                    loadCombined = object.load(file)
        
        if not loadCombined:
            mantaMsg('Could not load file ' + str(file))
    
    except Exception as e:
        mantaMsg('Exception in Python fluid file import: ' + str(e))
        pass # Just skip file load errors for now

def fluid_cache_get_framenr_formatted_7(framenr):
    return str(framenr).zfill(4) if framenr >= 0 else str(framenr).zfill(5)

def smoke_load_data_7(path, framenr, file_format, resumable):
    mantaMsg('Smoke load data')
    dict = { **fluid_data_dict_final_s7, **fluid_data_dict_resume_s7, **smoke_data_dict_final_s7, **smoke_data_dict_resume_s7 } if resumable else { **fluid_data_dict_final_s7, **smoke_data_dict_final_s7 }
    fluid_file_import_s7(dict=dict, path=path, framenr=framenr, file_format=file_format, file_name=file_data_s7)
    
    copyVec3ToReal(source=vel_s7, targetX=x_vel_s7, targetY=y_vel_s7, targetZ=z_vel_s7)

def smoke_load_noise_7(path, framenr, file_format, resumable):
    mantaMsg('Smoke load noise')
    dict = { **smoke_noise_dict_final_s7, **smoke_noise_dict_resume_s7 } if resumable else { **smoke_noise_dict_final_s7 } 
    fluid_file_import_s7(dict=dict, path=path, framenr=framenr, file_format=file_format, file_name=file_noise_s7)
    
    if resumable:
        # Fill up xyz texture grids, important when resuming a bake
        copyVec3ToReal(source=uvGrid0_s7, targetX=texture_u_s7, targetY=texture_v_s7, targetZ=texture_w_s7)
        copyVec3ToReal(source=uvGrid1_s7, targetX=texture_u2_s7, targetY=texture_v2_s7, targetZ=texture_w2_s7)

######################################################################
## PRE/POST STEPS
######################################################################

def fluid_pre_step_7():
    mantaMsg('Fluid pre step')
    
    phiObs_s7.setConst(9999)
    phiOut_s7.setConst(9999)
    
    # Main vel grid is copied in adapt time step function
    
    if using_obstacle_s7:
        # Average out velocities from multiple obstacle objects at one cell
        x_obvel_s7.safeDivide(numObs_s7)
        y_obvel_s7.safeDivide(numObs_s7)
        z_obvel_s7.safeDivide(numObs_s7)
        copyRealToVec3(sourceX=x_obvel_s7, sourceY=y_obvel_s7, sourceZ=z_obvel_s7, target=obvelC_s7)
    
    if using_invel_s7:
        copyRealToVec3(sourceX=x_invel_s7, sourceY=y_invel_s7, sourceZ=z_invel_s7, target=invelC_s7)
    
    if using_guiding_s7:
        weightGuide_s7.multConst(0)
        weightGuide_s7.addConst(alpha_sg7)
        interpolateMACGrid(source=guidevel_sg7, target=velT_s7)
        velT_s7.multConst(vec3(gamma_sg7))
    
    x_force_s7.multConst(scaleSpeedFrames_s7)
    y_force_s7.multConst(scaleSpeedFrames_s7)
    z_force_s7.multConst(scaleSpeedFrames_s7)
    copyRealToVec3(sourceX=x_force_s7, sourceY=y_force_s7, sourceZ=z_force_s7, target=forces_s7)
    
    # If obstacle has velocity, i.e. is a moving obstacle, switch to dynamic preconditioner
    if using_smoke_s7 and using_obstacle_s7 and obvelC_s7.getMax() > 0:
        mantaMsg('Using dynamic preconditioner')
        preconditioner_s7 = PcMGDynamic
    else:
        mantaMsg('Using static preconditioner')
        preconditioner_s7 = PcMGStatic

def fluid_post_step_7():
    mantaMsg('Fluid post step')
    
    # Copy vel grid to reals grids (which Blender internal will in turn use for vel access)
    copyVec3ToReal(source=vel_s7, targetX=x_vel_s7, targetY=y_vel_s7, targetZ=z_vel_s7)
    if using_guiding_s7:
        copyVec3ToReal(source=guidevel_sg7, targetX=x_guidevel_s7, targetY=y_guidevel_s7, targetZ=z_guidevel_s7)

######################################################################
## STEPS
######################################################################

def smoke_adaptive_step_7(framenr):
    mantaMsg('Manta step, frame ' + str(framenr))
    s7.frame = framenr
    
    fluid_pre_step_7()
    
    flags_s7.initDomain(boundaryWidth=0, phiWalls=phiObs_s7, outflow=boundConditions_s7)
    
    if using_obstacle_s7:
        mantaMsg('Extrapolating object velocity')
        # ensure velocities inside of obs object, slightly add obvels outside of obs object
        # extrapolate with phiObsIn before joining (static) phiObsSIn grid to prevent flows into static obs
        extrapolateVec3Simple(vel=obvelC_s7, phi=phiObsIn_s7, distance=6, inside=True)
        extrapolateVec3Simple(vel=obvelC_s7, phi=phiObsIn_s7, distance=3, inside=False)
        resampleVec3ToMac(source=obvelC_s7, target=obvel_s7)
        
        mantaMsg('Initializing obstacle levelset')
        phiObsIn_s7.join(phiObsSIn_s7) # Join static obstacle map
        phiObsIn_s7.floodFill(boundaryWidth=1)
        extrapolateLsSimple(phi=phiObsIn_s7, distance=6, inside=True)
        extrapolateLsSimple(phi=phiObsIn_s7, distance=3, inside=False)
        phiObs_s7.join(phiObsIn_s7)
        
        # Additional sanity check: fill holes in phiObs which can result after joining with phiObsIn
        phiObs_s7.floodFill(boundaryWidth=1)
        extrapolateLsSimple(phi=phiObs_s7, distance=6, inside=True)
        extrapolateLsSimple(phi=phiObs_s7, distance=3, inside=False)
    
    mantaMsg('Initializing fluid levelset')
    phiIn_s7.join(phiSIn_s7) # Join static flow map
    extrapolateLsSimple(phi=phiIn_s7, distance=6, inside=True)
    extrapolateLsSimple(phi=phiIn_s7, distance=3, inside=False)
    
    if using_outflow_s7:
        phiOutIn_s7.join(phiOutSIn_s7) # Join static outflow map
        phiOut_s7.join(phiOutIn_s7)
    
    setObstacleFlags(flags=flags_s7, phiObs=phiObs_s7, phiOut=phiOut_s7, phiIn=phiIn_s7, boundaryWidth=1)
    flags_s7.fillGrid()
    
    # reset emission accumulation at the beginning of an adaptive frame
    if not s7.timePerFrame:
        emission_s7.setConst(0.)
    # accumulate emission value per adaptive step for later use in noise computation
    emission_s7.join(emissionIn_s7)
    
    applyEmission(flags=flags_s7, target=density_s7, source=densityIn_s7, emissionTexture=emissionIn_s7, type=FlagInflow|FlagOutflow)
    if using_heat_s7:
        applyEmission(flags=flags_s7, target=heat_s7, source=heatIn_s7, emissionTexture=emissionIn_s7, type=FlagInflow|FlagOutflow)
    
    if using_colors_s7:
        applyEmission(flags=flags_s7, target=color_r_s7, source=color_r_in_s7, emissionTexture=emissionIn_s7, type=FlagInflow|FlagOutflow)
        applyEmission(flags=flags_s7, target=color_g_s7, source=color_g_in_s7, emissionTexture=emissionIn_s7, type=FlagInflow|FlagOutflow)
        applyEmission(flags=flags_s7, target=color_b_s7, source=color_b_in_s7, emissionTexture=emissionIn_s7, type=FlagInflow|FlagOutflow)
    
    if using_fire_s7:
        applyEmission(flags=flags_s7, target=fuel_s7, source=fuelIn_s7, emissionTexture=emissionIn_s7, type=FlagInflow|FlagOutflow)
        applyEmission(flags=flags_s7, target=react_s7, source=reactIn_s7, emissionTexture=emissionIn_s7, type=FlagInflow|FlagOutflow)
    
    mantaMsg('Smoke step / s7.frame: ' + str(s7.frame))
    if using_fire_s7:
        process_burn_7()
    smoke_step_7()
    if using_fire_s7:
        update_flame_7()
    
    s7.step()
    
    fluid_post_step_7()

def smoke_step_7():
    mantaMsg('Smoke step low')
    
    # save original state for later (used during noise creation)
    velTmp_s7.copyFrom(vel_s7)
    
    if using_dissolve_s7:
        mantaMsg('Dissolving smoke')
        dissolveSmoke(flags=flags_s7, density=density_s7, heat=heat_s7, red=color_r_s7, green=color_g_s7, blue=color_b_s7, speed=dissolveSpeed_s7, logFalloff=using_logdissolve_s7)
    
    mantaMsg('Advecting density')
    advectSemiLagrange(flags=flags_s7, vel=vel_s7, grid=density_s7, order=2)
    
    if using_heat_s7:
        mantaMsg('Advecting heat')
        advectSemiLagrange(flags=flags_s7, vel=vel_s7, grid=heat_s7, order=2)
    
    if using_fire_s7:
        mantaMsg('Advecting fire')
        advectSemiLagrange(flags=flags_s7, vel=vel_s7, grid=fuel_s7, order=2)
        advectSemiLagrange(flags=flags_s7, vel=vel_s7, grid=react_s7, order=2)
    
    if using_colors_s7:
        mantaMsg('Advecting colors')
        advectSemiLagrange(flags=flags_s7, vel=vel_s7, grid=color_r_s7, order=2)
        advectSemiLagrange(flags=flags_s7, vel=vel_s7, grid=color_g_s7, order=2)
        advectSemiLagrange(flags=flags_s7, vel=vel_s7, grid=color_b_s7, order=2)
    
    mantaMsg('Advecting velocity')
    advectSemiLagrange(flags=flags_s7, vel=vel_s7, grid=vel_s7, order=2)
    
    if not domainClosed_s7 or using_outflow_s7:
        resetOutflow(flags=flags_s7, real=density_s7)
    
    mantaMsg('Vorticity')
    if using_fire_s7:
        flame_s7.copyFrom(fuel_s7) # temporarily misuse flame grid as vorticity storage
        flame_s7.multConst(flameVorticity_s7)
    vorticityConfinement(vel=vel_s7, flags=flags_s7, strength=vorticity_s7, strengthCell=flame_s7 if using_fire_s7 else None)
    
    if using_heat_s7:
        mantaMsg('Adding heat buoyancy')
        addBuoyancy(flags=flags_s7, density=heat_s7, vel=vel_s7, gravity=gravity_s7, coefficient=buoyancy_heat_s7, scale=False)
    mantaMsg('Adding buoyancy')
    addBuoyancy(flags=flags_s7, density=density_s7, vel=vel_s7, gravity=gravity_s7, coefficient=buoyancy_dens_s7, scale=False)
    
    mantaMsg('Adding forces')
    addForceField(flags=flags_s7, vel=vel_s7, force=forces_s7)
    
    # Cells inside obstacle should not contain any density, fire, etc.
    if deleteInObstacle_s7:
        resetInObstacle(flags=flags_s7, density=density_s7, vel=vel_s7, heat=heat_s7, fuel=fuel_s7, flame=flame_s7, red=color_r_s7, green=color_g_s7, blue=color_b_s7)
    
    # add initial velocity
    if using_invel_s7:
        # Using cell centered invels, will be converted to MAC within the function
        setInitialVelocity(flags=flags_s7, vel=vel_s7, invel=invelC_s7)
    
    mantaMsg('Walls')
    setWallBcs(flags=flags_s7, vel=vel_s7, obvel=obvel_s7 if using_obstacle_s7 else None)
    
    preconditioner_s7 = PcMGDynamic if using_obstacle_s7 and obvel_s7.getMax() > 0 else PcMGStatic
    mantaMsg('Using preconditioner: ' + str(preconditioner_s7))
    if using_guiding_s7:
        mantaMsg('Guiding and pressure')
        PD_fluid_guiding(vel=vel_s7, velT=velT_s7, flags=flags_s7, weight=weightGuide_s7, blurRadius=beta_sg7, pressure=pressure_s7, tau=tau_sg7, sigma=sigma_sg7, theta=theta_sg7, preconditioner=preconditioner_s7, zeroPressureFixing=domainClosed_s7)
    else:
        mantaMsg('Pressure')
        solvePressure(flags=flags_s7, vel=vel_s7, pressure=pressure_s7, preconditioner=preconditioner_s7, zeroPressureFixing=domainClosed_s7) # closed domains require pressure fixing

def process_burn_7():
    mantaMsg('Process burn')
    processBurn(fuel=fuel_s7, density=density_s7, react=react_s7, red=color_r_s7, green=color_g_s7, blue=color_b_s7, heat=heat_s7, burningRate=burningRate_s7, flameSmoke=flameSmoke_s7, ignitionTemp=ignitionTemp_s7, maxTemp=maxTemp_s7, flameSmokeColor=flameSmokeColor_s7)

def update_flame_7():
    mantaMsg('Update flame')
    updateFlame(react=react_s7, flame=flame_s7)

def smoke_step_noise_7(framenr):
    mantaMsg('Manta step noise, frame ' + str(framenr))
    sn7.frame = framenr
    
    copyRealToVec3(sourceX=texture_u_s7, sourceY=texture_v_s7, sourceZ=texture_w_s7, target=uvGrid0_s7)
    copyRealToVec3(sourceX=texture_u2_s7, sourceY=texture_v2_s7, sourceZ=texture_w2_s7, target=uvGrid1_s7)
    
    flags_sn7.initDomain(boundaryWidth=0, phiWalls=phiObs_sn7, outflow=boundConditions_s7)
    
    mantaMsg('Interpolating grids')
    # Join big obstacle levelset after initDomain() call as it overwrites everything in phiObs
    if using_obstacle_s7:
        phiIn_sn7.copyFrom(phiObsIn_s7) if upres_sn7 <= 1 else interpolateGrid(target=phiIn_sn7, source=phiObsIn_s7) # mis-use phiIn_sn
        phiObs_sn7.join(phiIn_sn7)
    if using_outflow_s7:
        phiOut_sn7.copyFrom(phiOut_s7) if upres_sn7 <= 1 else interpolateGrid(target=phiOut_sn7, source=phiOut_s7)
    phiIn_sn7.copyFrom(phiIn_s7) if upres_sn7 <= 1 else interpolateGrid(target=phiIn_sn7, source=phiIn_s7)
    vel_sn7.copyFrom(velTmp_s7) if upres_sn7 <= 1 else interpolateMACGrid(target=vel_sn7, source=velTmp_s7)
    
    setObstacleFlags(flags=flags_sn7, phiObs=phiObs_sn7, phiOut=phiOut_sn7, phiIn=phiIn_sn7, boundaryWidth=1)
    flags_sn7.fillGrid()
    
    # Interpolate emission grids and apply them to big noise grids
    tmpIn_sn7.copyFrom(densityIn_s7) if upres_sn7 <= 1 else interpolateGrid(source=densityIn_s7, target=tmpIn_sn7)
    emissionIn_sn7.copyFrom(emission_s7) if upres_sn7 <= 1 else interpolateGrid(source=emission_s7, target=emissionIn_sn7)
    
    # Higher-res noise grid needs scaled emission values
    tmpIn_sn7.multConst(float(upres_sn7))
    applyEmission(flags=flags_sn7, target=density_sn7, source=tmpIn_sn7, emissionTexture=emissionIn_sn7, type=FlagInflow|FlagOutflow)
    
    if using_colors_s7:
        tmpIn_sn7.copyFrom(color_r_in_s7) if upres_sn7 <= 1 else interpolateGrid(source=color_r_in_s7, target=tmpIn_sn7)
        applyEmission(flags=flags_sn7, target=color_r_sn7, source=tmpIn_sn7, emissionTexture=emissionIn_sn7, type=FlagInflow|FlagOutflow)
        tmpIn_sn7.copyFrom(color_g_in_s7) if upres_sn7 <= 1 else interpolateGrid(source=color_g_in_s7, target=tmpIn_sn7)
        applyEmission(flags=flags_sn7, target=color_g_sn7, source=tmpIn_sn7, emissionTexture=emissionIn_sn7, type=FlagInflow|FlagOutflow)
        tmpIn_sn7.copyFrom(color_b_in_s7) if upres_sn7 <= 1 else interpolateGrid(source=color_b_in_s7, target=tmpIn_sn7)
        applyEmission(flags=flags_sn7, target=color_b_sn7, source=tmpIn_sn7, emissionTexture=emissionIn_sn7, type=FlagInflow|FlagOutflow)
    
    if using_fire_s7:
        tmpIn_sn7.copyFrom(fuelIn_s7) if upres_sn7 <= 1 else interpolateGrid(source=fuelIn_s7, target=tmpIn_sn7)
        applyEmission(flags=flags_sn7, target=fuel_sn7, source=tmpIn_sn7, emissionTexture=emissionIn_sn7, type=FlagInflow|FlagOutflow)
        tmpIn_sn7.copyFrom(reactIn_s7) if upres_sn7 <= 1 else interpolateGrid(source=reactIn_s7, target=tmpIn_sn7)
        applyEmission(flags=flags_sn7, target=react_sn7, source=tmpIn_sn7, emissionTexture=emissionIn_sn7, type=FlagInflow|FlagOutflow)
    
    mantaMsg('Noise step / sn7.frame: ' + str(sn7.frame))
    if using_fire_s7:
        process_burn_noise_7()
    step_noise_7()
    if using_fire_s7:
        update_flame_noise_7()
    
    sn7.step()
    
    copyVec3ToReal(source=uvGrid0_s7, targetX=texture_u_s7, targetY=texture_v_s7, targetZ=texture_w_s7)
    copyVec3ToReal(source=uvGrid1_s7, targetX=texture_u2_s7, targetY=texture_v2_s7, targetZ=texture_w2_s7)

def step_noise_7():
    mantaMsg('Smoke step noise')
    
    if using_dissolve_s7:
        mantaMsg('Dissolving noise')
        dissolveSmoke(flags=flags_sn7, density=density_sn7, heat=None, red=color_r_sn7, green=color_g_sn7, blue=color_b_sn7, speed=dissolveSpeed_s7, logFalloff=using_logdissolve_s7)
    
    mantaMsg('Advecting UVs and updating UV weight')
    advectSemiLagrange(flags=flags_s7, vel=vel_s7, grid=uvGrid0_s7, order=2)
    updateUvWeight(resetTime=sn7.timestep*10.0 , index=0, numUvs=uvs_s7, uv=uvGrid0_s7, offset=uvs_offset_s7)
    advectSemiLagrange(flags=flags_s7, vel=vel_s7, grid=uvGrid1_s7, order=2)
    updateUvWeight(resetTime=sn7.timestep*10.0 , index=1, numUvs=uvs_s7, uv=uvGrid1_s7, offset=uvs_offset_s7)
    
    if not domainClosed_s7 or using_outflow_s7:
        resetOutflow(flags=flags_sn7, real=density_sn7)
    
    mantaMsg('Energy')
    computeEnergy(flags=flags_s7, vel=vel_s7, energy=energy_s7)
    
    tmpFlags_s7.copyFrom(flags_s7)
    extrapolateSimpleFlags(flags=flags_s7, val=tmpFlags_s7, distance=2, flagFrom=FlagObstacle, flagTo=FlagFluid)
    extrapolateSimpleFlags(flags=tmpFlags_s7, val=energy_s7, distance=6, flagFrom=FlagFluid, flagTo=FlagObstacle)
    computeWaveletCoeffs(energy_s7)
    
    sStr_s7 = 1.0 * wltStrength_s7
    sPos_s7 = 2.0
    
    mantaMsg('Applying noise vec')
    for o in range(octaves_s7):
        uvWeight_s7 = getUvWeight(uvGrid0_s7)
        applyNoiseVec3(flags=flags_sn7, target=vel_sn7, noise=wltnoise_sn7, scale=sStr_s7 * uvWeight_s7, scaleSpatial=sPos_s7 , weight=energy_s7, uv=uvGrid0_s7)
        uvWeight_s7 = getUvWeight(uvGrid1_s7)
        applyNoiseVec3(flags=flags_sn7, target=vel_sn7, noise=wltnoise_sn7, scale=sStr_s7 * uvWeight_s7, scaleSpatial=sPos_s7 , weight=energy_s7, uv=uvGrid1_s7)
        
        sStr_s7 *= 0.06 # magic kolmogorov factor 
        sPos_s7 *= 2.0 
    
    for substep in range(int(upres_sn7)):
        if using_colors_s7: 
            mantaMsg('Advecting colors noise')
            advectSemiLagrange(flags=flags_sn7, vel=vel_sn7, grid=color_r_sn7, order=2)
            advectSemiLagrange(flags=flags_sn7, vel=vel_sn7, grid=color_g_sn7, order=2)
            advectSemiLagrange(flags=flags_sn7, vel=vel_sn7, grid=color_b_sn7, order=2)
        
        if using_fire_s7: 
            mantaMsg('Advecting fire noise')
            advectSemiLagrange(flags=flags_sn7, vel=vel_sn7, grid=fuel_sn7, order=2)
            advectSemiLagrange(flags=flags_sn7, vel=vel_sn7, grid=react_sn7, order=2)
        
        mantaMsg('Advecting density noise')
        advectSemiLagrange(flags=flags_sn7, vel=vel_sn7, grid=density_sn7, order=2)

def process_burn_noise_7():
    mantaMsg('Process burn noise')
    processBurn(fuel=fuel_sn7, density=density_sn7, react=react_sn7, red=color_r_sn7, green=color_g_sn7, blue=color_b_sn7, burningRate=burningRate_s7, flameSmoke=flameSmoke_s7, ignitionTemp=ignitionTemp_s7, maxTemp=maxTemp_s7, flameSmokeColor=flameSmokeColor_s7)

def update_flame_noise_7():
    mantaMsg('Update flame noise')
    updateFlame(react=react_sn7, flame=flame_sn7)

######################################################################
## MAIN
######################################################################

# Helper function to call cache load functions
def load_data(frame, cache_resumable):
    smoke_load_data_7(os.path.join(cache_dir, 'data'), frame, file_format_data, cache_resumable)
    if using_noise_s7:
        smoke_load_noise_7(os.path.join(cache_dir, 'noise'), frame, file_format_data, cache_resumable)
    if using_guiding_s7:
        fluid_load_guiding_7(os.path.join(cache_dir, 'guiding'), frame, file_format_data)

# Helper function to call step functions
def step(frame):
    smoke_adaptive_step_7(frame)
    if using_noise_s7:
        smoke_step_noise_7(frame)

gui = None
if (GUI):
    gui=Gui()
    gui.show()
    gui.pause()

cache_resumable       = True
cache_dir             = '/home/z/dev/Downloads/cache_fluid_blender_mantatest'
file_format_data      = '.vdb'
file_format_mesh      = '.bobj.gz'

# How many frame to load from cache
from_cache_count = 100

loop_count = 0
while current_frame_s7 <= end_frame_s7:
    
    # Load already simulated data from cache:
    if loop_count < from_cache_count:
        load_data(current_frame_s7, cache_resumable)
    
    # Otherwise simulate new data
    else:
        while(s7.frame <= current_frame_s7):
            if using_adaptTime_s7:
                fluid_adapt_time_step_7()
            step(current_frame_s7)
    
    current_frame_s7 += 1
    loop_count += 1
    
    if gui:
        gui.pause()
