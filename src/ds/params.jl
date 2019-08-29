"""
    Params

GPUE.jl structure for all readonly values

Entries include simulation constants, settings, and coordinate grids
"""
mutable struct Params
  dimnum::Integer

  xDim::Integer
  yDim::Integer
  zDim::Integer

  dx::Float64
  dy::Float64
  dz::Float64

  xMax::Float64
  yMax::Float64
  zMax::Float64

  omega::Float64
  omegaX::Float64
  omegaY::Float64
  omegaZ::Float64

  x::CuArray{Float64}
  y::CuArray{Float64}
  z::CuArray{Float64}

  nAtoms::Integer
  mass::Float64
  scatterLen::Float64
  g::Float64

  a0x::Float64
  a0y::Float64
  a0z::Float64

  Rxy::Float64

  winding::Float64

  px::CuArray{Complex{Float64}}
  py::CuArray{Complex{Float64}}
  pz::CuArray{Complex{Float64}}

  compression::Integer

  dt::Complex{Float64}

  gstate::Bool

  iterations::Integer
  printSteps::Integer
  writeOut::Bool
end

"""
    Params(; kwargs...)

Constructs the `Params` structure.
Takes only keyword arguments and returns a new `Params` object with the given and derived data.

# Arguments
- `xDim::Integer=256`: The width of the simulation.
- `yDim::Integer=256`: The height of the simulation. For 1D, set equal to 1.
- `zDim::Integer=1`: The depth of the simulation. For 2D or below, set equal to 1.
- `boxSize::Float64`: The size of the simulation boundary length. Will default to `6*Rxy*a0{x,y,z}`.
- `omega::Float64=0.0`: The rotation coefficient for gauge field simulations.
- `omegaX::Float64=2π`: The rotation coefficient for gauge field simulation in the X dimension.
- `omegaY::Float64=2π`: The rotation coefficient for gauge field simulation in the Y dimension.
- `omegaZ::Float64=2π`: The rotation coefficient for gauge field simulation in the Z dimension
- `winding::Float64=0.0`: The scaling factor for induced phase winding in the initial wave function.
- `compression::Integer=6` The compression level for HDF5 output. Must be in the range [0, 9].
- `dt::Float64=1e-4`: The timestep size. For groundstate simulation, divide `dt` by the imaginary unit.
- `nAtoms::Integer=1`: The number of atoms to simulate.
- `mass::Float64=1.4431607e-25`: The mass of the particle (in kg), defaults to that of Rubidium-87.
- `scatterLen::Float64=4.76e-9`: The scattering length of the particle.
- `iterations::Integer=1`: The number of simulation steps to run.
- `printSteps::Integer=100`: The frequency of writing the simulation data to file.
- `writeOut::Bool=true`: The condition for writing simulation data to file.

"""
function Params(; xDim=256, yDim=256, zDim=1, boxSize=0.0, omega=0.0, omegaX=2*pi, omegaY=2*pi, omegaZ=2*pi, winding=0.0, compression=6, dt=1e-4, nAtoms=1, mass=1.4431607e-25, scatterLen=4.76e-9, iterations=1, printSteps=100, writeOut=true)
  if yDim == zDim == 1
    dimnum = 1
  elseif zDim == 1
    dimnum = 2
  else
    dimnum = 3
  end

  a0x = sqrt(ħ / (2 * mass * omegaX))
  a0y = sqrt(ħ / (2 * mass * omegaY))
  a0z = sqrt(ħ / (2 * mass * omegaZ))

  Rxy = (15.0 * nAtoms * scatterLen * sqrt(mass * omegaZ / ħ)) ^ 0.2

  g = 4.0 * nAtoms * ħ * ħ * pi * scatterLen / mass
  if dimnum == 2
    g *= sqrt(mass * omegaZ / (2 * pi * ħ))
  end

  if boxSize > 0
    xMax = yMax = zMax = boxSize
  else
    xMax = 6 * Rxy * a0x
    yMax = 6 * Rxy * a0y
    zMax = 6 * Rxy * a0z
  end

  dx = dimnum < 1 ? 1 : 2 * xMax / xDim
  dy = dimnum < 2 ? 1 : 2 * yMax / yDim
  dz = dimnum < 3 ? 1 : 2 * zMax / zDim

  dpx = pi / xMax
  dpy = pi / yMax
  dpz = pi / zMax

  pxMax = dpx * (xDim / 2)
  pyMax = dpy * (yDim / 2)
  pzMax = dpz * (zDim / 2)

  gstate = real(dt) != dt

  x = dimnum < 1 ? [0] : reshape(collect(-xMax + dx/2 : dx : xMax), (xDim,))
  y = dimnum < 2 ? [0] : reshape(collect(-yMax + dy/2 : dy : yMax), (1, yDim))
  z = dimnum < 3 ? [0] : reshape(collect(-zMax + dz/2 : dz : zMax), (1, 1, zDim))

  px = dimnum < 1 ? [0] : reshape(vcat(0 : dpx : pxMax - dpx, -pxMax : dpx : -dpx / 2), (xDim,))
  py = dimnum < 2 ? [0] : reshape(vcat(0 : dpy : pyMax - dpy, -pyMax : dpy : -dpy / 2), (1, yDim))
  pz = dimnum < 3 ? [0] : reshape(vcat(0 : dpz : pzMax - dpz, -pzMax : dpz : -dpz / 2), (1, 1, zDim))

  return Params(dimnum,
            xDim, yDim, zDim,
            dx, dy, dz,
            xMax, yMax, zMax,
            omega, omegaX, omegaY, omegaZ,
            CuArray(x), CuArray(y), CuArray(z),
            nAtoms, mass, scatterLen, g,
            a0x, a0y, a0z, Rxy, winding,
            CuArray(complex(px)), CuArray(complex(py)), CuArray(complex(pz)),
            compression,
            dt,
            gstate,
            iterations, printSteps, writeOut)
end

"""Helper function to add an item to dictionary if it is not missing"""
function addIf!(dict, name, value, fallbackDict)
  if !ismissing(value)
    dict[name] = value
  else
    try
      dict[name] = fallbackDict[name]
    catch
    end
  end
end

"""
    Params(f::FileData; kwargs...)

Wrapper from the Params constructor, loading from file when possible
"""
function Params(f::FileData; xDim=missing, yDim=missing, zDim=missing, boxSize=missing, omega=missing, omegaX=missing, omegaY=missing, omegaZ=missing, winding=missing, compression=missing, dt=missing, nAtoms=missing, mass=missing, scatterLen=missing, iterations=missing, printSteps=missing, writeOut=missing)
  # Get all stored params from the file
  paramDict = loadParams(f)

  out = Dict()

  # Add every explicitly given kwarg into the kwargs dict, overwriting the loaded data from file
  for (expr, value) in [(:xDim, xDim), (:yDim, yDim), (:zDim, zDim), (:boxSize, boxSize), (:omega, omega), (:omegaX, omegaX), (:omegaY, omegaY), (:omegaZ, omegaZ), (:winding, winding), (:compression, compression), (:dt, dt), (:nAtoms, nAtoms), (:mass, mass), (:scatterLen, scatterLen), (:iterations, iterations), (:printSteps, printSteps), (:writeOut, writeOut)]
    addIf!(out, expr, value, paramDict)
  end

  return Params(; out...)
end

