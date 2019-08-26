using GPUE

function main()
  f = initFileData()

  par = Params(xDim=512, yDim=512, printSteps=2500, iterations=50000, dt=1e-4im, omega=0.6, nAtoms=1e5, omegaZ=100)
  opr = Operators(par)
  aux = Aux(par, opr)

  try
    evolve!(f, par, opr, aux)
  finally
    terminate!(f, par, opr, aux)
  end
end

main()
