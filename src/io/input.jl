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

