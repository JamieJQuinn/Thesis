import sys
import math

import sdf
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from plotting import latexify, save_plot
from plotting_parameters import *
from parameters import *

INCLUDE_CBARS=True
INCLUDE_TITLE=True
INCLUDE_AXIS_LABELS=True

def load_sdf_file(sdfFilename):
    return sdf.read(sdfFilename)

def get_dx_dy_dz(sdfFile, cell_centre=False):
    extents = sdfFile.Grid_Grid.extents
    dims = sdfFile.Grid_Grid.dims

    if cell_centre:
        dims = [dim - 1 for dim in dims]

    return [(extents[i+3] - extents[i])/dims[i] for i in range(3)]

def get_variable(sdfFile, variable_name):
    return getattr(sdfFile, variable_name)

def get_variable_data(sdfFile, variable_name):
    if variable_name == "magnitude_current_density":
        return get_magnitude_current(sdfFile)
    elif variable_name == "vorticity_density":
        return get_magnitude_vorticity(sdfFile)
    elif variable_name == "kinetic_energy":
        return calc_kinetic_energy(sdfFile)
    elif variable_name == "kinetic_energy_z":
        return calc_kinetic_energy_z(sdfFile)
    elif variable_name == "abs_Velocity_Vz":
        return np.abs(get_variable_data(sdfFile, "Velocity_Vz"))
    elif variable_name == "alfven_velocity":
        return calc_alfven_velocity(sdfFile)
    elif variable_name == "sound_speed":
        return calc_sound_speed(sdfFile)
    elif variable_name == "pressure":
        return calc_pressure(sdfFile)
    elif variable_name == "magnetic_pressure":
        return calc_magnetic_pressure(sdfFile)
    elif variable_name == "parallel_electric_field":
        return calc_parallel_electric_field(sdfFile)
    else:
        return get_variable(sdfFile, variable_name).data

def calc_centred_velocity(var):
    x_av = var[1:,:,:] + var[:-1,:,:]
    y_av = x_av[:,1:,:] + x_av[:,:-1,:]
    z_av = y_av[:,:,1:] + y_av[:,:,:-1]
    return z_av / 8.0

def calc_centred_velocity_slice(var):
    x_av = var[1:,:] + var[:-1,:]
    y_av = x_av[:,1:] + x_av[:,:-1]
    return y_av / 4.0

def calc_pressure(sdfFile):
    rho = get_variable_data(sdfFile, "Fluid_Rho")
    energy = get_variable_data(sdfFile, "Fluid_Energy")
    pressure = rho*energy * (GAMMA - 1.0)
    return pressure

def calc_magnetic_pressure(sdfFile):
    bx = get_variable_data(sdfFile, "Magnetic_Field_bx_centred")
    by = get_variable_data(sdfFile, "Magnetic_Field_by_centred")
    bz = get_variable_data(sdfFile, "Magnetic_Field_bz_centred")

    return np.power(bx, 2) + np.power(by, 2) + np.power(bz, 2)

def calc_mag_velocity2(sdfFile):
    vx = get_variable_data(sdfFile, "Velocity_Vx")
    vy = get_variable_data(sdfFile, "Velocity_Vy")
    vz = get_variable_data(sdfFile, "Velocity_Vz")

    return np.power(vx, 2) + np.power(vy, 2) + np.power(vz, 2)

def calc_parallel_electric_field(sdfFile):
    bx = get_variable_data(sdfFile, "Magnetic_Field_bx_centred")
    by = get_variable_data(sdfFile, "Magnetic_Field_by_centred")
    bz = get_variable_data(sdfFile, "Magnetic_Field_bz_centred")

    dx, dy, dz = get_dx_dy_dz(sdfFile, cell_centre=True)

    mag_mag = np.sqrt(np.power(bx,2) + np.power(by,2) + np.power(bz,2))

    gradbx = np.gradient(bx, dx, dy, dz)
    gradby = np.gradient(by, dx, dy, dz)
    gradbz = np.gradient(bz, dx, dy, dz)

    return (bx*(gradbz[1] - gradby[2]) + by*(gradbz[0] - gradbx[2]) + bz*(gradby[0] - gradbx[1]))/mag_mag

def calc_sound_speed(sdfFile):
    # From Newton-Laplace formula
    # sound speed = sqrt(temperature)
    return np.sqrt(calc_temperature(sdfFile))

def calc_temperature(sdfFile):
    # Nondimensional temperature is just 1-gamma * internal energy
    return (GAMMA-1.0)*get_variable_data(sdfFile, "Fluid_Energy")

def calc_kinetic_energy(sdfFile):
    magV2 = calc_centred_velocity(
        calc_mag_velocity2(sdfFile)
    )

    rho = get_variable_data(sdfFile, "Fluid_Rho")

    return 0.5*rho*magV2

def calc_kinetic_energy_z(sdfFile):
    vz = calc_centred_velocity(
        get_variable_data(sdfFile, "Velocity_Vz")
    )

    rho = get_variable_data(sdfFile, "Fluid_Rho")

    return 0.5*rho*np.power(vz, 2)

def calc_mag(vector):
    return np.power(vector[0], 2) + np.power(vector[1], 2) + np.power(vector[2], 2)

def calc_alfven_velocity(sdfFile):
    B = calc_mag(get_magnetic_field(sdfFile))
    rho = get_variable_data(sdfFile, "Fluid_Rho")
    return B/np.sqrt(rho)

def get_magnetic_field(sdfFile):
    mag_field = np.array([
        get_variable_data(sdfFile, variable) for variable in
        ['Magnetic_Field_bx_centred', 'Magnetic_Field_by_centred', 'Magnetic_Field_bz_centred']
    ])

    return mag_field

def attach_colorbar(axis, im, side='right'):
    divider = make_axes_locatable(axis)
    cax = divider.append_axes(side, size="5%", pad=0.05)
    if side == 'right' or side =='left':
        orientation = 'vertical'
    else:
        orientation = 'horizontal'
    plt.colorbar(im, cax=cax, orientation=orientation)

def plot_slice(sdfFile, variable_name, dimension, slice_loc,
               cbar=INCLUDE_CBARS,
               include_title=INCLUDE_TITLE,
               include_axis_labels=INCLUDE_AXIS_LABELS,
               xlim=False, ylim=False, title=False, cmap="viridis",
              vmax = None, vmin = None):
    velocity = get_slice(sdfFile, variable_name, dimension, slice_loc)
    extents = get_slice_extents(sdfFile, dimension)

    latexify(columns=1)
    fig, axis = plt.subplots()

    im = axis.imshow(velocity.T,\
                vmax=vmax, vmin=vmin,\
                interpolation='bilinear',\
                cmap=plt.get_cmap(cmap),
               extent=extents, origin='lower')

    if xlim:
        axis.set_xlim(xlim)
    if ylim:
        axis.set_ylim(ylim)

    if include_title:
        if title:
            axis.title.set_text(" ".join(title.split("_")))
        else:
            axis.title.set_text(" ".join(variable_name.split("_")))

    if include_axis_labels:
        labels = get_axis_labels(dimension)
        axis.set_xlabel(labels[0])
        axis.set_ylabel(labels[1])

    if cbar:
        attach_colorbar(axis, im)


def get_axis_labels(dimension):
    if type(dimension) is str:
        dimension = get_dimension_index(dimension)
    labels = ["x", "y", "z"]
    labels.pop(dimension)
    return labels


def get_slice_extents(sdfFile, dimension):
    if type(dimension) is str:
        dimension = get_dimension_index(dimension)

    extents = list(sdfFile.Grid_Grid.extents)
    extents.pop(dimension+3)
    extents.pop(dimension)
    extents[1], extents[2] = extents[2], extents[1]
    return extents

def get_extents(sdfFile):
    return list(sdfFile.Grid_Grid.extents)

def get_slice(sdfFile, variable_name, dimension, slice_loc):
    if type(dimension) is str:
        dimension = get_dimension_index(dimension)

    index = length_to_index(sdfFile, slice_loc, dimension)
    data = get_variable_data(sdfFile, variable_name)

    return np.take(data, index, axis=dimension)


def slice_variable(sdf_var, x_min=0, x_max=-1, y_min=0, y_max=-1, z_min=0, z_max=-1):
    if(x_max < 0):
        x_max = sdf_var.dims[0]
    if(y_max < 0):
        y_max = sdf_var.dims[1]
    if(z_max < 0):
        z_max = sdf_var.dims[2]
    
    return sdf_var.data[x_min:x_max, y_min:y_max, z_min:z_max]

def curl(bx, by, bz, dx, dy, dz):
    gradbx = np.gradient(bx, dx, dy, dz)
    gradby = np.gradient(by, dx, dy, dz)
    gradbz = np.gradient(bz, dx, dy, dz)

    return np.array([(gradbz[1] - gradby[2]), -(gradbz[0] - gradbx[2]), (gradby[0] - gradbx[1])])

def mag_curl(bx, by, bz, dx, dy, dz):
    gradbx = np.gradient(bx, dx, dy, dz)
    gradby = np.gradient(by, dx, dy, dz)
    gradbz = np.gradient(bz, dx, dy, dz)

    return np.sqrt(np.power(gradbz[1] - gradby[2], 2) + np.power(gradbz[0] - gradbx[2], 2) + np.power(gradby[0] - gradbx[1], 2))

def get_magnitude_vorticity(sdfFile):
    vx = get_variable_data(sdfFile, "Velocity_Vx")
    vy = get_variable_data(sdfFile, "Velocity_Vy")
    vz = get_variable_data(sdfFile, "Velocity_Vz")

    dx, dy, dz = get_dx_dy_dz(sdfFile)

    return mag_curl(vx, vy, vz, dx, dy, dz)

def centre_magnetic_field(bx, by, bz):
    bx = 0.5*(bx[:-1,:,:] + bx[1:,:,:])
    by = 0.5*(by[:,:-1,:] + by[:,1:,:])
    bz = 0.5*(bz[:,:,:-1] + bz[:,:,1:])

    return bx, by, bz

def get_magnitude_current(sdfFile):
    bx = get_variable_data(sdfFile, "Magnetic_Field_Bx")
    by = get_variable_data(sdfFile, "Magnetic_Field_By")
    bz = get_variable_data(sdfFile, "Magnetic_Field_Bz")

    bx, by, bz = centre_magnetic_field(bx, by, bz)

    dx, dy, dz = get_dx_dy_dz(sdfFile, cell_centre=True)

    return mag_curl(bx, by, bz, dx, dy, dz)

def get_magnitude_current_at(sdfFile, zSliceIdx, xLimits = (0,-1), yLimits=(0,-1)):
    bx = slice_variable(\
        get_variable(sdfFile, "Magnetic_Field_bx_centred"),\
        x_min=xLimits[0], x_max = xLimits[1],\
        y_min=yLimits[0], y_max = yLimits[1],\
        z_min=zSliceIdx - 1, z_max = zSliceIdx + 2)
    by = slice_variable(\
        get_variable(sdfFile, "Magnetic_Field_by_centred"),\
        x_min=xLimits[0], x_max = xLimits[1],\
        y_min=yLimits[0], y_max = yLimits[1],\
        z_min=zSliceIdx - 1, z_max = zSliceIdx + 2)
    bz = slice_variable(\
        get_variable(sdfFile, "Magnetic_Field_bz_centred"),\
        x_min=xLimits[0], x_max = xLimits[1],\
        y_min=yLimits[0], y_max = yLimits[1],\
        z_min=zSliceIdx - 1, z_max = zSliceIdx + 2)
    
    dx, dy, dz = get_dx_dy_dz(sdfFile)

    gradbx = np.gradient(bx, dx, dy, dz)
    gradby = np.gradient(by, dx, dy, dz)
    gradbz = np.gradient(bz, dx, dy, dz)

    current_density = np.sqrt(np.power(gradbz[1] - gradby[2], 2) + np.power(gradbz[0] - gradbx[2], 2) + np.power(gradby[0] - gradbx[1], 2))
    current_density = current_density[:,::-1,int(current_density.shape[2]/2)]
    return current_density.transpose()

def get_temperature_at(sdfFile, zSliceIdx, xLimits = (0,-1), yLimits=(0,-1)):
    # Temperature in nondim units is just \gamma - 1 times the internal energy
    temp = (GAMMA-1) * slice_variable(\
        get_variable(sdfFile, "Fluid_Energy"),\
        x_min=xLimits[0], x_max = xLimits[1],\
        y_min=yLimits[0], y_max = yLimits[1],\
        z_min=zSliceIdx, z_max = zSliceIdx+1)
    
    temp = temp[:,:].squeeze(axis=2).transpose()
    
    return temp

def get_dimension_index(dimension):
    if dimension == "x":
        return 0
    elif dimension == "y":
        return 1
    elif dimension == "z":
        return 2

def length_to_index(sdfFile, x, dimension):
    if type(dimension) is str:
        dimension = get_dimension_index(dimension)
    
    extents = sdfFile.Grid_Grid.extents
    dims = sdfFile.Grid_Grid.dims
    
    x0 = extents[dimension]
    xN = extents[dimension + 3]
    N = dims[dimension]
    
    return int((x - x0)/(xN-x0) * (N-1))

def dimensionalise_temperature(tempIn):
    B0 = 5
    mf = 1.2
    mh_si = 1.672621777
    kb_si = 1.3806488
    mu0_si = 4.0 * np.pi
    RHO0 = 1.67
    T0 = (B0*B0)*mf*mh_si/(kb_si*mu0_si*RHO0) * 1e9

    return T0*tempIn
