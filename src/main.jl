using GPUE

f = initFileData()
par = Params()
opr = Operators(par)
aux = Aux(opr)

try
  evolve(f, par, opr, aux)
finally
  terminate(f, par, opr, aux)
end
