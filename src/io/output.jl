struct FileData
  file::HDF5File

  wfc::HDF5Group
  wfc_const::HDF5Group
  wfc_ev::HDF5Group
  v::HDF5Group
  k::HDF5Group

  domain::HDF5Group
  
  datasets::Dict{String, HDF5Dataset}
end

## Creating new output

function writeAttributes(f::FileData, par::Params)
  for sym in propertynames(par)
    val = getproperty(par, sym)
    if typeof(val) <: HDF5.HDF5Scalar
      attrs(f.file)[string(sym)] = val
    end
  end
end

function writeWfc(f::FileData, par::Params, opr::Operators, aux::Aux)
  if (par.gstate)
    group = f.wfc_const
  else
    group = f.wfc_ev
  end
  
  group[string(aux.i), "chunk", par.chunks, "compress", par.compression] = Array(opr.wfc)
end

function writeV(f::FileData, par::Params, opr::Operators, aux::Aux)
  group = f.v
  
  group[string(aux.i), "chunk", par.chunks, "compress", par.compression] = Array(opr.V)
end

function writeK(f::FileData, par::Params, opr::Operators, aux::Aux)
  group = f.k

  group[string(aux.i), "chunk", par.chunks, "compress", par.compression] = Array(opr.K)
end

function writeAxes(f::FileData, par::Params)
  group = f.domain

  group["X"] = Array(par.x)
  group["Y"] = Array(par.y)
  group["Z"] = Array(par.z)
end

function initFileData()
  # WARNING: This has yet to be merged into HDF5.jl
  HDF5.set_complex_field_names("re", "im")

  file = h5open("./output.h5", "w")

  wfc = g_create(file, "WFC")
  wfc_const = g_create(wfc, "CONST")
  wfc_ev = g_create(wfc, "EV")

  v = g_create(file, "V")
  k = g_create(file, "K")

  domain = g_create(file, "DOMAIN")

  return FileData(file, wfc, wfc_const, wfc_ev, v, k, domain, Dict())
end


## Loading previous output

function loadGroup(f::FileData)

end

function loadAttribute(f::FileData)

end

function loadDataSet(f::FileData)

end

function loadFileData()

end


## Conclude operations

function closeFile(f::FileData)
  close(f.file)
end
