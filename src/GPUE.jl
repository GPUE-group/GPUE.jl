"""
    GPUE.jl

Gross—Pitaevskii equation solver that is accelerated on GPU hardware with CUDA.

GPUE.jl has been developed as a maintainable alternative to [GPUE](https://github.com/GPUE-group/GPUE),
which functions similarly and is written in CUDA/C++.
For more on the original GPUE, see the [documentation](https://github.com/GPUE-group/GPUE)
"""
module GPUE
  # Required GPU libraries
  using CuArrays, CUDAnative

  # Import FFT definitions from CuArrays through AbstractFFTs
  using AbstractFFTs

  # Import I/O functionality through HDF5.jl
  using HDF5

  # Mathematical and Physical Constants
  const ħ = 1.0545718176461565e-34

  
  # Implement GPU wrappers for operations on complex numbers
  # Can be removed when CUDAnative issue #430 gets fixed
  @inline CUDAnative.exp(x::Complex{Float64}) = CUDAnative.exp(x.re) * (CUDAnative.cos(x.im) + 1.0im * CUDAnative.sin(x.im))
  @inline CUDAnative.abs(x::Complex{Float64}) = CUDAnative.hypot(x.re, x.im)
  @inline CUDAnative.abs2(x::Complex{Float64}) = x.re * x.re + x.im * x.im

  # Rename atan2 to atan as per julia standard,
  # Can be removed when CUDAnative PR #443 gets merged
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

