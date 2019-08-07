struct FileData
  file::HDF5File

  wfc::HDF5Group
  wfc_const::HDF5Group
  wfc_ev::HDF5Group
  v::HDF5Group
  k::HDF5Group

  domain::HDF5Group
  x::HDF5Group
  y::HDF5Group
  z::HDF5Group
  
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

function initFileData()
  file = h5open("./output.h5", "w")

  wfc = g_create(file, "WFC")
  wfc_const = g_create(wfc, "CONST")
  wfc_ev = g_create(wfc, "EV")

  v = g_create(file, "V")
  k = g_create(file, "K")

  domain = g_create(file, "DOMAIN")
  x = g_create(domain, "X")
  y = g_create(domain, "Y")
  z = g_create(domain, "Z")

  return FileData(file, wfc, wfc_const, wfc_ev, v, k, domain, x, y, z, Dict())
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
