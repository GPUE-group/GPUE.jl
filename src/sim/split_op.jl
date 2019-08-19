function normalize!(par::Params, opr::Operators, aux::Aux)
  opr.wfc ./= sqrt(sum(mapreduce(abs2, +, opr.wfc; init=1.0im)) * par.dx * par.dy * par.dz)
end

function split_op!(par::Params, opr::Operators, aux::Aux)

  # Half-step in real space
  @. opr.wfc *= opr.V * exp(1.0im * par.g * abs2(opr.wfc) * par.dt / ħ)

  # FFT to momentum space
  aux.forward_plan_all * opr.wfc

  # Full step in momentum space
  @. opr.wfc *= opr.K

  # iFFT back
  aux.inverse_plan_all * opr.wfc

  # Final half-step in real space
  @. opr.wfc *= opr.V * exp(1.0im * par.g * abs2(opr.wfc) * par.dt / ħ)

  # Renormalize if we are in imaginary time
  if par.gstate
    normalize!(par, opr, aux)
  end
end
