import sympy
from joblib import Parallel, delayed
import time

def AddOne(v):
    return v + 1

vec = []
for i in range(2000):
    vec.append(sympy.symbols('x_'+str(i)))
print(vec)

#works
t = time.time()
newVec = Parallel(n_jobs=1)(delayed(AddOne)(i) for i in vec)
print(time.time() - t)
print(newVec)

t = time.time()
newVec = Parallel(n_jobs=-1, backend="threading")(delayed(AddOne)(i) for i in vec)
print(time.time() - t)
print(newVec)
