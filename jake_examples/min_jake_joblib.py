from cypari2 import Pari
from joblib import Parallel, delayed
import time

def AddOne(v):
    return v + pari.one()

pari = Pari()
vec = []
for i in range(20000):
    vec.append(pari('x_' + str(i+1)))
#print(vec)

#works
t = time.time()
newVec = Parallel(n_jobs=1)(delayed(AddOne)(i) for i in vec)
print(time.time() - t)
#print(newVec)

#doesn't work
t = time.time()
newVec2 = Parallel(n_jobs=-1, backend="threading")(delayed(AddOne)(i) for i in vec)
print(time.time() - t)
#print(newVec2)
