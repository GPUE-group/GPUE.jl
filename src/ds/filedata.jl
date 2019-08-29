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
  
  aux::HDF5Group
end

# Constructors

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

  aux = g_create(file, "AUX")

  return FileData(file, wfc, wfc_const, wfc_ev, v, k, p, px, py, pz, a, ax, ay, az, domain, aux)
end

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

  aux = openGroup(file, "AUX")

  return FileData(file, wfc, wfc_const, wfc_ev, v, k, p, px, py, pz, a, ax, ay, az, domain, aux)
end

# Destructor

"""
    closeFile(f::FileData)

Close the file and terminate file operations.
"""
function closeFile(f::FileData)
  if isopen(f.file)
    close(f.file)
  end
end

