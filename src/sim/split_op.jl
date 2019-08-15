function update_density!(par::Params, opr::Operators, aux::Aux)
  aux.density = sum(abs2.(opr.wfc))
end

function normalize!(par::Params, opr::Operators, aux::Aux)
  opr.wfc ./= sqrt(aux.density * par.dx * par.dy * par.dz)
end

function split_op!(par::Params, opr::Operators, aux::Aux)
  update_density!(par, opr, aux)
  density_opr = exp(im * par.g * aux.density * par.dt / ħ)

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
    update_density!(par, opr, aux)
    normalize!(par, opr, aux)
  end
end
