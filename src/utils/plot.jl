using GPUE
using Plots
using HDF5

GPUE.HDF5.set_complex_field_names("re", "im")

f = loadFileData()
gstate = read(attrs(f.file)["gstate"])

group = gstate ? f.wfc_const : f.wfc_ev

anim = @animate for name in names(group)[2:end]
  println("Writing frame: $name")
  dset = group[name]
  data = read(dset)
  heatmap(abs.(data))
  annotate!([(0, -20, "$name")])
end

gif(anim, "./output/wfc.gif", fps=5)

GPUE.closeFile(f)
