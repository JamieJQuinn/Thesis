def start_document(f):
    f.write("\\documentclass[12pt]{article}\n")
    f.write("\\usepackage{graphicx}\n")
    f.write("\\usepackage[a4paper, margin=0cm]{geometry}\n")
    f.write("\\usepackage{grffile}\n")
    f.write("\\begin{document}\n")

def end_document(f):
    f.write("\\end{document}\n")

def write_tex(texname, slice_options, time_range, visc_range, resist_range):
    print("Writing", texname)
    with open(texname, "w") as f:
        start_document(f)
        for time in time_range:
            timedump = '{0:04d}'.format(time)
            for resist in resist_range:
                for visc in visc_range:
                    # for model in ["isotropic", "switching"]:
                    for slice_option in slice_options:
                        for variable_name in slice_option["variables"]:
                            slice_dim = slice_option["slice_dim"]
                            slice_loc = slice_option["slice_loc"]
                            f.write("\\begin{figure}\n")
                            run_folder = "v-" + str(visc) + "r-" + str(resist) + "-isotropic"
                            outname = run_folder + "_" + variable_name + "_" + slice_dim + "_" + str(slice_loc) + "_" + timedump + ".pdf"
                            f.write("\\includegraphics[width=0.5\\textwidth]{" + outname + "}\n")
                            run_folder = "v-" + str(visc) + "r-" + str(resist) + "-switching"
                            outname = run_folder + "_" + variable_name + "_" + slice_dim + "_" + str(slice_loc) + "_" + timedump + ".pdf"
                            f.write("\\includegraphics[width=0.5\\textwidth]{" + outname + "}\n")
                            f.write("\\end{figure}\n")
                            # f.write("![v](" + outname + ")\n")
                    f.write("\\clearpage\n")
        end_document(f)

visc_range = [3,4,5]
resist_range = [3,4,5]
time_range = range(1,16)

z_slices = {"slice_dim":"z",
            "slice_loc":"0.0",
            "variables":["Velocity_Vz"]
            # "variables":["Velocity_Vz", "vorticity_density", "magnitude_current_density", "Fluid_Energy"]
           }
x_0_0_slices = {"slice_dim":"x",
            "slice_loc":"0.0",
            "variables":["Velocity_Vx"]
            # "variables":["Velocity_Vx", "vorticity_density", "magnitude_current_density", "Fluid_Energy"]
           }
x_0_85_slices = {"slice_dim":"x",
            "slice_loc":"0.85",
            "variables":["Velocity_Vz", "Velocity_Vy"]
            # "variables":["Velocity_Vz", "vorticity_density", "magnitude_current_density", "Velocity_Vy", "Fluid_Energy"]
           }

slice_options = [z_slices, x_0_0_slices, x_0_85_slices]

write_tex("velocity_overview.tex", slice_options, time_range, visc_range, resist_range)

for so in slice_options:
    so["variables"] = ["vorticity_density"]

write_tex("vorticity_overview.tex", slice_options, time_range, visc_range, resist_range)

for so in slice_options:
    so["variables"] = ["magnitude_current_density"]

write_tex("magnitude_current_density_overview.tex", slice_options, time_range, visc_range, resist_range)

for so in slice_options:
    so["variables"] = ["Fluid_Energy"]

write_tex("energy_overview.tex", slice_options, time_range, visc_range, resist_range)
