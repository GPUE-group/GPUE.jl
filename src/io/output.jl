struct FileData
  file::HDF5File

  wfc::HDF5Group
  wfc_const::HDF5Group
  wfc_ev::HDF5Group

  v::HDF5Group
  k::HDF5Group

  a::HDF5Group
  ax::HDF5Group
  ay::HDF5Group
  az::HDF5Group

  domain::HDF5Group
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

function getChunks(x)
  map(elem -> Int(floor(sqrt(elem))), size(x))
end

function writeWfc(f::FileData, par::Params, opr::Operators, aux::Aux)
  if (par.gstate)
    group = f.wfc_const
  else
    group = f.wfc_ev
  end
  
  group[string(aux.i), "chunk", getChunks(opr.wfc), "compress", par.compression] = Array(opr.wfc)
end

function writeV(f::FileData, par::Params, opr::Operators, aux::Aux)
  group = f.v
  
  group[string(aux.i), "chunk", getChunks(opr.V), "compress", par.compression] = Array(opr.V)
end

function writeK(f::FileData, par::Params, opr::Operators, aux::Aux)
  group = f.k

  group[string(aux.i), "chunk", getChunks(opr.K), "compress", par.compression] = Array(opr.K)
end

function writeGauge(f::FileData, par::Params, opr::Operators, aux::Aux)
  f.ax[string(aux.i), "chunk", getChunks(opr.Ax), "compress", par.compression] = Array(opr.Ax)
  f.ay[string(aux.i), "chunk", getChunks(opr.Ay), "compress", par.compression] = Array(opr.Ay)
  f.az[string(aux.i), "chunk", getChunks(opr.Az), "compress", par.compression] = Array(opr.Az)
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

  file = h5open("./output/output.h5", "w")

  wfc = g_create(file, "WFC")
  wfc_const = g_create(wfc, "CONST")
  wfc_ev = g_create(wfc, "EV")

  v = g_create(file, "V")
  k = g_create(file, "K")

  a = g_create(file, "A")
  ax = g_create(a, "AX")
  ay = g_create(a, "AY")
  az = g_create(a, "AZ")

  domain = g_create(file, "DOMAIN")

  return FileData(file, wfc, wfc_const, wfc_ev, v, k, a, ax, ay, az, domain)
end

## Loading previous output

function loadFileData()
  file = h5open("./output/output.h5", "r+")

  wfc = file["WFC"]
  wfc_const = wfc["CONST"]
  wfc_ev = wfc["EV"]

  v = file["V"]
  k = file["K"]

  a = file["A"]
  ax = a["AX"]
  ay = a["AY"]
  az = a["AZ"]

  domain = file["DOMAIN"]

  return FileData(file, wfc, wfc_const, wfc_ev, v, k, a, ax, ay, az, domain)
end


## Conclude operations

function closeFile(f::FileData)
  if isopen(f.file)
    close(f.file)
  end
end
