mutable struct Aux

  forward_plan::CuArrays.CUFFT.cCuFFTPlan{Complex{Float64}}
  inverse_plan::CuArrays.CUFFT.cCuFFTPlan{Complex{Float64}}

  function Aux(forward_plan, inverse_plan)
    new(forward_plan, inverse_plan)
  end

end
