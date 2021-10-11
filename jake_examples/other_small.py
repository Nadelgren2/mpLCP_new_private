from cypari2 import Pari
from joblib import Parallel, delayed

def AddOne(v):
    return v + pari.one()

pari = Pari()
vec = [pari('x_1'), pari('x_2')]
print(vec)

#works
newVec = Parallel(n_jobs=1)(delayed(AddOne)(i) for i in vec)
print(newVec)

print(pari.parapply(lambda i: i + pari.one(), vec))
