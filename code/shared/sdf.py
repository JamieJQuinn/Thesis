import sys
import math

import sdf
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

sys.path.insert(0, '..')
from plotting_functions import latexify, save_plot
from plotting_parameters import *

INCLUDE_CBARS=True
INCLUDE_TITLE=True
INCLUDE_AXIS_LABELS=True

def get_var_slice(sdfFilename, variable):
    sdfFile = sdf.read(sdfFilename)

    sdf_var = getattr(sdfFile, variable)
    var = sdf_var.data
    dims = sdf_var.dims

    # If z-dimensions even, average between midpoints
    average = (dims[2] % 2 == 0)

    if average:
        return 0.5*(var[:,:,int(dims[2]/2)-1] + var[:,:,int(dims[2]/2)]).transpose()
    else:
        return var[:,:,int(dims[2]/2)].transpose()

def get_dx_dy_dz(sdfFile, cell_centre=False):
    extents = sdfFile.Grid_Grid.extents
    dims = sdfFile.Grid_Grid.dims

    if cell_centre:
        dims = [dim - 1 for dim in dims]

    return [(extents[i+3] - extents[i])/dims[i] for i in range(3)]

def get_variable(sdfFile, varname):
    return getattr(sdfFile, varname)

def get_magnetic_field(sdfFilename):
    sdfFile = sdf.read(sdfFilename)

    mag_field = np.array([
        getattr(sdfFile, variable).data for variable in
        ['Magnetic_Field_bx_centred', 'Magnetic_Field_by_centred', 'Magnetic_Field_bz_centred']
    ])

    return mag_field

def plot_slice(sdfFilename, variable_name, dimension, slice_loc,
               cbar=INCLUDE_CBARS,
               include_title=INCLUDE_TITLE,
               include_axis_labels=INCLUDE_AXIS_LABELS,
               xlim=False, ylim=False, title=False):
    velocity = get_slice(sdfFilename, variable_name, dimension, slice_loc)
    extents = get_slice_extents(sdfFilename, dimension)

    latexify(columns=1)
    fig, axis = plt.subplots()

    im = axis.imshow(velocity.T,\
    #             vmax=v_limits, vmin=-v_limits,\
                interpolation='bilinear',\
                cmap=plt.get_cmap("viridis"),
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
        divider = make_axes_locatable(axis)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(im, cax=cax)


def get_axis_labels(dimension):
    if type(dimension) is str:
        dimension = get_dimension_index(dimension)
    labels = ["x", "y", "z"]
    labels.pop(dimension)
    return labels


def get_slice_extents(sdfFilename, dimension):
    if type(dimension) is str:
        dimension = get_dimension_index(dimension)

    sdfFile = sdf.read(sdfFilename)
    extents = list(sdfFile.Grid_Grid.extents)
    extents.pop(dimension+3)
    extents.pop(dimension)
    extents[1], extents[2] = extents[2], extents[1]
    return extents


def get_slice(sdfFilename, variable_name, dimension, slice_loc):
    sdfFile = sdf.read(sdfFilename)
    
    if variable_name == "magnitude_current_density":
        data = get_magnitude_current(sdfFile)
    elif variable_name == "vorticity_density":
        data = get_magnitude_vorticity(sdfFile)
    else:
        var = get_variable(sdfFile, variable_name)
        data = var.data
    
    if type(dimension) is str:
        dimension = get_dimension_index(dimension)
    
    index = length_to_index(sdfFile, slice_loc, dimension)
    
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
    vx = get_variable(sdfFile, "Velocity_Vx").data
    vy = get_variable(sdfFile, "Velocity_Vy").data
    vz = get_variable(sdfFile, "Velocity_Vz").data

    dx, dy, dz = get_dx_dy_dz(sdfFile)

    return mag_curl(vx, vy, vz, dx, dy, dz)

def centre_magnetic_field(bx, by, bz):
    bx = 0.5*(bx[:-1,:,:] + bx[1:,:,:])
    by = 0.5*(by[:,:-1,:] + by[:,1:,:])
    bz = 0.5*(bz[:,:,:-1] + bz[:,:,1:])

    return bx, by, bz

def get_magnitude_current(sdfFile):
    bx = get_variable(sdfFile, "Magnetic_Field_Bx").data
    by = get_variable(sdfFile, "Magnetic_Field_By").data
    bz = get_variable(sdfFile, "Magnetic_Field_Bz").data

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
    gamma = 5.0/3
    temp = (gamma-1) * slice_variable(\
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
