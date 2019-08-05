module GPUE
  using CuArrays, FFTW

  const ħ = 1.0545718176461565e-34


  # Load data structures

  include("ds/params.jl")
  export Params

  include("ds/operators.jl")
  export Operators

  include("ds/aux.jl")
  export Aux


  # Load simulation functions

  include("sim/split_op.jl")
  export split_op!

end

