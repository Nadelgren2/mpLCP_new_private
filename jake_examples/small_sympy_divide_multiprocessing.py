import sympy
import multiprocessing
from joblib import Parallel, delayed
import time

def Div(v):
    num, den = v
    return num/den

vec = []
for i in range(20000):
    vec.append(sympy.symbols('x_'+str(i)))
#print(vec)

t = time.time()
if __name__ == '__main__':
    pool = multiprocessing.Pool(processes = None) ## doesn't matter how many I use
    newVec = pool.map(Div, ((i, vec[0]) for i in vec))
    #print(newVec)
print(time.time() - t)

temp = vec[0]
t = time.time()
for i in range(len(vec)):
    vec[i] = vec[i]/temp
print(time.time() - t)
print(vec[0:10])
