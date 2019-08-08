module GPUE
  using CuArrays, CUDAnative, FFTW
  using HDF5

  const ħ = 1.0545718176461565e-34

  
  # Temporary hack until CUDAnative issue #430 gets fixed
  @inline CUDAnative.exp(x::Complex{Float64}) = CUDAnative.exp(x.re) * (CUDAnative.cos(x.im) + 1.0im * CUDAnative.sin(x.im))

  # Temporary hack until CUDAnative PR #443 gets merged
  @inline CUDAnative.atan(y, x) = CUDAnative.atan2(y, x)

  # Load data structures

  include("ds/params.jl")
  export Params

  include("ds/operators.jl")
  export Operators

  include("ds/aux.jl")
  export Aux

  # Load I/O functions

  include("io/output.jl")
  export FileData
  export initFileData, loadFileData

  # Load simulation functions

  include("sim/split_op.jl")

  include("sim/gauge.jl")

  include("sim/evolution.jl")
  export evolve

  # Load management functions

  include("sim/terminate.jl")
  export terminate
end

