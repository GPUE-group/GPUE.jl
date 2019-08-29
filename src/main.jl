using GPUE

function main()
  f = loadFileData()

  par = Params(f; xDim=512, yDim=512, printSteps=2500, iterations=100000, dt=1e-4/im, omega=0.6, nAtoms=1e5, omegaZ=100)
  opr = Operators(f, par)
  aux = Aux(f, par, opr)

  try
    evolve!(f, par, opr, aux)
  finally
    terminate!(f, par, opr, aux)
  end
end

main()

