"""
    FileData

GPUE.jl structure for file-related IO through HDF5.jl

Stores the file and the groups contained within it.
"""
struct FileData
  file::HDF5File

  wfc::HDF5Group
  wfc_const::HDF5Group
  wfc_ev::HDF5Group

  v::HDF5Group
  k::HDF5Group

  p::HDF5Group
  px::HDF5Group
  py::HDF5Group
  pz::HDF5Group

  a::HDF5Group
  ax::HDF5Group
  ay::HDF5Group
  az::HDF5Group

  domain::HDF5Group
end


## Creating new output

"""
    writeAttributes(f::FileData, par::Params)

Write all scalar simulation parameters to the file as
attributes on the file's root group.
"""
function writeAttributes(f::FileData, par::Params)
  for sym in propertynames(par)
    val = getproperty(par, sym)
    if typeof(val) <: HDF5.HDF5Scalar
      attrs(f.file)[string(sym)] = val
    end
  end
end

"""Helper to determine the size of the chunks given a nd-array."""
function getChunks(x)
  map(elem -> Int(floor(sqrt(elem))), size(x))
end

"""Helper to pad the output of the dataset index to a fixed number of characters."""
function zpad(i)
  lpad(i, 6, "0")
end

"""
    writeWfc(f::FileData, par::Params, opr::Operators, aux::Aux)

Write the wave function to the output file, either in group 
wfc_const (`/WFC/CONST/`)
or wfc_ev (`/WFC/EV/`)
for imaginary and real time respectively.

Each call writes the wave function as a new dataset with name corresponding
to the zero-padded iteration index (from 0).
"""
function writeWfc(f::FileData, par::Params, opr::Operators, aux::Aux)
  if (par.gstate)
    group = f.wfc_const
  else
    group = f.wfc_ev
  end
  println("Writing wfc $(aux.i) to file")
  group[zpad(aux.i), "chunk", getChunks(opr.wfc), "compress", par.compression] = Array(opr.wfc)
end

"""
    writeV(f::FileData, par::Params, opr::Operators, aux::Aux)

Write the real-space trapping potential values to the output file, in group
v (`/V/`).

Each call writes the values as a new dataset with name corresponding to the
zero-padded iteration index (from 0).
"""
function writeV(f::FileData, par::Params, opr::Operators, aux::Aux)
  group = f.v
  
  group[zpad(aux.i), "chunk", getChunks(opr.v), "compress", par.compression] = opr.v
end

"""
    writeK(f::FileData, par::Params, opr::Operators, aux::Aux)

Write the momentum-space values to the output file, in group k (`/K/`).
Write each of the momentum domains (px, py, pz) to the groups:
px (`/P/PX`), py (`/P/PY`), and px (`/P/PZ`).

Each is written as a new dataset with name corresponding to the zero-padded
iteration index (from 0).
"""
function writeK(f::FileData, par::Params, opr::Operators, aux::Aux)
  f.k[zpad(aux.i), "chunk", getChunks(opr.k), "compress", par.compression] = opr.k

  f.px[zpad(aux.i), "chunk", getChunks(par.px), "compress", par.compression] = Array(par.px)
  f.py[zpad(aux.i), "chunk", getChunks(par.py), "compress", par.compression] = Array(par.py)
  f.pz[zpad(aux.i), "chunk", getChunks(par.pz), "compress", par.compression] = Array(par.pz)
end

"""
    writeGauge(f::FileData, par::Params, opr::Operators, aux::Aux)

Write the gauge field values to the output file, in the groups:
ax (`/A/AX/`), ay (`/A/AY/`), az (`/A/AZ/`).

Each is written as a new dataset with the name corresponding to the
zero-padded iteration index (from 0).
"""
function writeGauge(f::FileData, par::Params, opr::Operators, aux::Aux)
  f.ax[zpad(aux.i), "chunk", getChunks(opr.ax), "compress", par.compression] = opr.ax
  f.ay[zpad(aux.i), "chunk", getChunks(opr.ay), "compress", par.compression] = opr.ay
  f.az[zpad(aux.i), "chunk", getChunks(opr.az), "compress", par.compression] = opr.az
end

"""
    writeAxes(f::FileData, par::Params)

Write the coordinate grid axes to the output file, in the group:
domain (`/DOMAIN/`)

Each is written as a new dataset with the name `X`, `Y`, or `Z`.
"""
function writeAxes(f::FileData, par::Params)
  group = f.domain

  group["X"] = Array(par.x)
  group["Y"] = Array(par.y)
  group["Z"] = Array(par.z)
end

"""
    initFileData(file_name="data.h5")

Initialize a new FileData structure with the appropriate groups.
The file is created with the given name in the output/ folder, and will
overwrite any existing file (to load a file, see loadFileData).
"""
function initFileData(file_name="data.h5")
  HDF5.set_complex_field_names("re", "im")

  file = h5open("./output/$file_name", "w")

  wfc = g_create(file, "WFC")
  wfc_const = g_create(wfc, "CONST")
  wfc_ev = g_create(wfc, "EV")

  v = g_create(file, "V")
  k = g_create(file, "K")

  p = g_create(file, "P")
  px = g_create(p, "PX")
  py = g_create(p, "PY")
  pz = g_create(p, "PZ")
  
  a = g_create(file, "A")
  ax = g_create(a, "AX")
  ay = g_create(a, "AY")
  az = g_create(a, "AZ")

  domain = g_create(file, "DOMAIN")

  return FileData(file, wfc, wfc_const, wfc_ev, v, k, p, px, py, pz, a, ax, ay, az, domain)
end


## Loading previous output

"""Helper to open a group if it exists"""
function openGroup(parent, name)
  if exists(parent, name)
    return parent[name]
  else
    return g_create(parent, name)
  end
end

"""
    loadFileData(file_name="data.h5")

Initialize a new FileData structure with the appropriate groups.
The file is loaded from that with the given name in the output folder,
and groups are opened, or created if they do not exist.
"""
function loadFileData(file_name="data.h5")
  HDF5.set_complex_field_names("re", "im")
  file = h5open("./output/$file_name", "r+")

  wfc = openGroup(file, "WFC")
  wfc_const = openGroup(wfc, "CONST")
  wfc_ev = openGroup(wfc, "EV")

  v = openGroup(file, "V")
  k = openGroup(file, "K")

  p = openGroup(file, "P")
  px = openGroup(p, "PX")
  py = openGroup(p, "PY")
  pz = openGroup(p, "PZ")

  a = openGroup(file, "A")
  ax = openGroup(a, "AX")
  ay = openGroup(a, "AY")
  az = openGroup(a, "AZ")

  domain = openGroup(file, "DOMAIN")

  return FileData(file, wfc, wfc_const, wfc_ev, v, k, p, px, py, pz, a, ax, ay, az, domain)
end

"""
    loadParams(f::FileData, par::Params)

Load the params from the data file
"""
function loadParams(f::FileData, par::Params)
  for name in names(attrs(f.file))
    val = f.file[name].value
    setproperty!(par, Symbol(name), val)
  end
end

"""Helper function to load the data from the latest dataset from a group"""
function loadDset(group)
  if length(names(group)) < 1
    return missing
  end
  iter_str = last(sort(names(group)))
  dset = group[iter_str]

  return read(dset)
end

"""
    loadWfc(f::FileData, par::Params)

Load the wave function from the data file
"""
function loadWfc(f::FileData, par::Params)
  group = par.gstate ? f.wfc_const : f.wfc_ev
  return loadDset(group)
end

"""
    loadV(f::FileData, par::Params)

Load the real-space trapping potential values from the data file
"""
function loadV(f::FileData, par::Params)
  return loadDset(f.v)
end

"""
    loadK(f::FileData, par::Params)

Load the momentum-space values from the data file
"""
function loadK(f::FileData, par::Params)
  return map(loadDset, [f.k, f.px, f.py, f.pz])
end

"""
    loadA(f::FileData, par::Params)

Load the gauge field values from the data file
"""
function loadA(f::FileData, par::Params)
  return map(loadDset, [f.ax, f.ay, f.az])
end

## Conclude operations

"""
    closeFile(f::FileData)

Close the file and terminate file operations.
"""
function closeFile(f::FileData)
  if isopen(f.file)
    close(f.file)
  end
end

