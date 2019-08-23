"""Normalize the wave function to maintain density stability"""
function normalize!(par::Params, opr::Operators, aux::Aux)
  opr.wfc ./= sqrt(sum(mapreduce(abs2, +, opr.wfc; init=1.0im)) * par.dx * par.dy * par.dz)
end

"""
    split_op!(par::Params, opr::Operators, aux::Aux)

Perform a single Split-Step Fourier Method (SSFM) iteration.
Procedure is as follows:
- Multiply by the half the position space operators (V and g|ψ|^2)
- Apply a Fourier transform to move to momentum space
- Multiply by the momentum space operator (K)
- Apply an inverse Fourier transform to move back into position space
- Multiply by the rest of the position space operators
- If running in imaginary time, normalize to keep density fixed

For more information, see the
[official gpue documentation](https://gpue-group.github.io/intro/),
or the [Algorithm Archive page on the Split-Operator Method]
(https://www.algorithm-archive.org/contents/split-operator_method/split-operator_method.html)
"""
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
