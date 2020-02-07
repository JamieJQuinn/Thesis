# Anisotropic viscosity with astrophysical applications

This is the working repo for my thesis, structured in the following way:

- `build/` containing the final latex outputs
- `imagea/` containing the build scripts and outputs for all images
- `chapters/` containing the .tex files for each chapter
- `data-analysis/` containing any miscellaneous analysis not resulting in images
- `makefile` the global makefile for the entire thesis

## Data analysis (not resulting in images)

For analysis outputs like tables that do no result in images, the code will reside in the data-analysis folder within this repo.

## Images

The images will be structured in the following way:

Inside the `images` folder will be a list of projects, roughly corresponding to chapters:
`images/lare1d`, `images/slow_null_point`, `images/kink_instability`, etc

Within the projects involving Lare3d there will be the following:
- Slices produced via python inside one jupyter notebook per chapter-project
`images/project/slices`
- Energy graphs produced via python inside one jupyter notebook per chapter
`images/project/energy_plots`
- 3D field line plots produced via visit inside one folder per simulation
`images/project/3d_field_line_plots`
- Field line integrations produced via Mayavi + python
`images/field_line_integrations`
- Other plots based on the .sdf or en.dat files will be produced in a misc folder with corresponding source code
`images/project/misc`

`images/project/diagrams`

Other outputs for non-lare focussed chapters may have other outputs and folders will be assigned appropriately. Diagrams produced without some source code (e.g. from inkscape) will be held in a single folder per project.

### Guidance on committing

While editing the images, only the source code should be committed. Once the image is finalised in some respect (still to be finalised by supervisors), it should be committed. Since nearly every image should be reproducible, it's not a huge loss to totally remove an image from git's memory if it becomes too large. An exception should be made for those images that are not reproducible, i.e. those produced as diagrams in inkscape.

### On including images

The images may be included with more simplicity using the `graphicspath` command in latex:

```
\graphicspath{{images/}{images/project/}{etc.}}
```

## On intermediate files

The plotting field line integrations requires some intermediate files to be created. Firstly, the raw SDF files must be sliced to reduce file sizes and increase the resolution of the integrator. Secondly, the integrator itself produces a 2D dataset that can be analysed as secondary data. I believe these files should be retained in folders with a structure that matches that of the original data files and all these files stored adjacent to the raw files, probably in a structure that looks like the following:

data/project/type_of_simulation/lare3d_folder/secondary_analysis/ --> sliced_sdfs
                                                                  \-> integrated_quantities

The production of these files will be controlled by the field line integrations images makefile. This allows multiple projects/presentations to access these files, keeps the makefiles simple by using relative paths but keeps the large secondary data files out of the images folder.
