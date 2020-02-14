import sys
import math

import sdf
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

sys.path.insert(0, '..')
from plotting_functions import latexify, save_plot
from plotting_parameters import *

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

def get_dx_dy_dz(sdfFile):
    extents = sdfFile.Grid_Grid.extents
    dims = sdfFile.Grid_Grid.dims
    
    dims = [dim + 1 for dim in dims]
    
    return [(extents[i+3] - extents[i])/dims[i] for i in range(3)]

def get_magnetic_field(sdfFilename):
    sdfFile = sdf.read(sdfFilename)
    
    mag_field = np.array([
        getattr(sdfFile, variable).data for variable in 
        ['Magnetic_Field_bx_centred', 'Magnetic_Field_by_centred', 'Magnetic_Field_bz_centred']
    ])

    return mag_field

def get_variable(sdfFile, varname):
    return getattr(sdfFile, varname)
    
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

def length_to_index(length):
    # This is hardcoded - could be made better by including limits from the actual sdf file
    return int((length + 2)*350/4)

def dimensionalise_temperature(tempIn):
    B0 = 5
    mf = 1.2
    mh_si = 1.672621777
    kb_si = 1.3806488
    mu0_si = 4.0 * np.pi
    RHO0 = 1.67
    T0 = (B0*B0)*mf*mh_si/(kb_si*mu0_si*RHO0) * 1e9

    return T0*tempIn

def plot_current_density_and_velocity(sdfFile, outname, xLimits=(-2, 2), yLimits=(-2, 2),\
                                      max_current=-1, max_temperature=-1, print_maxes=False,\
                                     remove_axes=False):
    latexify(columns=2, square=True)
    fig,axis = plt.subplots()

    # if not remove_axes:
        # # Make room for labels
        # fig.subplots_adjust(left=.12, bottom=0.05, right=.965, top=.99)
    # else:
        # fig.subplots_adjust(left=0, bottom=0, right=1, top=1)

    xIndexLimits = tuple((length_to_index(length) for length in xLimits))
    yIndexLimits = tuple((length_to_index(length) for length in yLimits))
    
    current_density = get_magnitude_current_at(sdfFile, 350, xLimits=xIndexLimits, yLimits=yIndexLimits)
    if max_current < 0:
        max_current = int(current_density.max())
    
    axis.imshow(current_density,\
                extent=[xLimits[0], xLimits[1], yLimits[0], yLimits[1]],\
                vmax=max_current, vmin=0,\
                interpolation='bilinear',\
                cmap=plt.get_cmap("Blues"))
    if print_maxes:
        axis.text(-1, 0.9, "Max Current: " + str(max_current), color='w')
    
    print("Max current:", max_current)
    
    temperature = get_temperature_at(sdfFile, 350, xLimits=xIndexLimits, yLimits=yIndexLimits)
    print(dimensionalise_temperature(temperature[0,0]))
    if max_temperature < 0:
        max_temperature = int(dimensionalise_temperature(temperature.max()))
    
    axis.imshow(dimensionalise_temperature(temperature),\
                extent=[xLimits[0], xLimits[1], -yLimits[1], yLimits[0]],\
                vmax=max_temperature, vmin=2e4,\
                interpolation='bilinear',\
                cmap=plt.get_cmap("Reds"))
    
    if print_maxes:
        axis.text(-1, -0.95, "Max temp: " + str(int(max_temperature/1e6)), color='w')
    
    print(r"Max temperature $\times 10^{-6}$:", max_temperature/1e6)
            
    vx = slice_variable(\
        get_variable(sdfFile, "Velocity_Vx"),\
        z_min=350, z_max = 350 + 1)
    vy = slice_variable(\
        get_variable(sdfFile, "Velocity_Vy"),\
        z_min=350, z_max = 350 + 1)
    
    vx = 0.25*(vx[1:,1:] + vx[:-1,1:] + vx[:-1,:-1] + vx[1:,:-1])
    vy = 0.25*(vy[1:,1:] + vy[:-1,1:] + vy[:-1,:-1] + vy[1:,:-1])
    
    vectorIdx = 6
    
    vx = vx.squeeze(axis=2)[\
                            length_to_index(xLimits[0])+int(vectorIdx/2):length_to_index(xLimits[1]):vectorIdx,\
                            length_to_index(yLimits[1])-int(vectorIdx/2):length_to_index(yLimits[0]):-vectorIdx\
                           ].transpose()
    vy = -vy.squeeze(axis=2)[\
                            length_to_index(xLimits[0])+int(vectorIdx/2):length_to_index(xLimits[1]):vectorIdx,\
                            length_to_index(yLimits[1])-int(vectorIdx/2):length_to_index(yLimits[0]):-vectorIdx\
                           ].transpose()
    
    X, Y = np.meshgrid(\
                      np.linspace(xLimits[0], xLimits[1], vx.shape[1]),\
                      np.linspace(-yLimits[1], yLimits[0], vx.shape[0])\
                      )
    
    if remove_axes:
        axis.tick_params(axis='both', which='both', bottom=False, top=False, labelbottom=False, right=False, left=False, labelleft=False)
    
#     axis.streamplot(X, Y, vx, vy)
    axis.quiver(X, Y, vx, vy, units='xy', pivot='tail', width=0.007, color='k')
    
    save_plot(outname)
    plt.show()
    plt.close(fig)
