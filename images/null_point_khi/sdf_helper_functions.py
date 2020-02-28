import numpy as np
import sdf

def get_magnetic_field(sdfFilename):
    sdfFile = sdf.read(sdfFilename)
    
    mag_field = np.array([
        getattr(sdfFile, variable).data for variable in 
        ['Magnetic_Field_bx_centred', 'Magnetic_Field_by_centred', 'Magnetic_Field_bz_centred']
    ])

    return mag_field
