import sympy
from joblib import Parallel, delayed
import time

def Div(num,den):
    return num/den

vec = []
for i in range(20000):
    vec.append(sympy.symbols('x_'+str(i)))
#print(vec)

#works
t = time.time()
newVec = Parallel(n_jobs=1)(delayed(Div)(i, vec[0]) for i in vec)
print(time.time() - t)
#print(newVec)

t = time.time()
newVec = Parallel(n_jobs=-1, backend="threading")(delayed(Div)(i, vec[0]) for i in vec)
print(time.time() - t)
print(newVec[0:10])

temp = vec[0]
t = time.time()
for i in range(len(vec)):
    vec[i] = vec[i]/temp
print(time.time() - t)
print(vec[0:10])
