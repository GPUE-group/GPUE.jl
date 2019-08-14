using GPUE
using Plots
using HDF5

GPUE.HDF5.set_complex_field_names("re", "im")

f = loadFileData()
gstate = read(attrs(f.file)["gstate"])

group = gstate ? f.wfc_const : f.wfc_ev

anim = @animate for name in names(group)
  dset = group[name]
  data = read(dset)
  heatmap(abs.(data), clims=(0, 2e5))
end

gif(anim, "./output/wfc.gif", fps=5)

GPUE.closeFile(f)
