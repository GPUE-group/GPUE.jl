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

"""
    writeAux(f::FileData, aux::Aux)

Write all scalar auxiliary variables (namely the iteration counter `i`)
to the file as attributes on the file's `AUX` group.
"""
function writeAux(f::FileData, aux::Aux)
  for sym in propertynames(aux)
    val = getproperty(aux, sym)
    if typeof(val) <: HDF5.HDF5Scalar
      if exists(attrs(f.aux), string(sym))
        a_delete(f.aux, string(sym))
      end
      if string(sym) == "i"
        val += 1
      end
      attrs(f.aux)[string(sym)] = val
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

