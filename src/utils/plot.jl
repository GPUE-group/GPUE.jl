using GPUE
using Plots
GPUE.HDF5.set_complex_field_names("re", "im")

f = loadFileData()

anim = @animate for name in names(f.wfc_const)
  dset = f.wfc_const[name]
  wfc = read(dset)
  heatmap(abs.(wfc))
end

gif(anim, "./wfc.gif", fps=5)

