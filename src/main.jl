using GPUE


try
  global f = initFileData()
  global par = Params(iterations=10000)
  global opr = Operators(par)
  global aux = Aux(par, opr)

  evolve(f, par, opr, aux)
finally
  terminate(f, par, opr, aux)
end
