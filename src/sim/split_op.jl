function normalize!(par::Params, opr::Operators, aux::Aux)
  opr.wfc ./= sqrt(sum(abs2.(opr.wfc)) * par.dx * par.dy * par.dz)
end

function split_op!(par::Params, opr::Operators, aux::Aux)
  density = abs2.(opr.wfc)
  density_opr = @. exp(im * par.g * density * par.dt / ħ)

  # Half-step in real space
  @. opr.wfc *= opr.V * density_opr

  # FFT to momentum space
  aux.forward_plan_all * opr.wfc

  # Full step in momentum space
  @. opr.wfc *= opr.K

  # iFFT back
  aux.inverse_plan_all * opr.wfc

  # Final half-step in real space
  @. opr.wfc *= opr.V * density_opr

  # Renormalize if we are in imaginary time
  if par.gstate
    normalize!(par, opr, aux)
  end
end
