from cypari2 import Pari
from joblib import Parallel, delayed
import multiprocessing
import time

def AddOne(v):
    return v + pari.one()

pari = Pari()
vec = []
for i in range(50000):
    vec.append(pari('x_' + str(i+1)))
#print(vec)

vec2 = []
t = time.time()
for v in vec:
    vec2.append(AddOne(v))
print("Serial For Loop")
print(time.time() - t)
#print(vec2)

#works
t = time.time()
newVec = Parallel(n_jobs=1)(delayed(AddOne)(i) for i in vec)
print("Joblib time -- 1 thread")
print(time.time() - t)
#print(newVec)

#doesn't work
t = time.time()
newVec2 = Parallel(n_jobs=-1, backend="threading")(delayed(AddOne)(i) for i in vec)
print("Joblib time -- all threads")
print(time.time() - t)
#print(newVec2)

#doesn't work
if __name__ == '__main__':
    t = time.time()
    pool = multiprocessing.Pool(processes = None) ## doesn't matter how many I use
    newVec = pool.map(AddOne, (i for i in vec))
    print("Multiprocessing time -- all threads")
    print(time.time() - t)
#    print(newVec)

if __name__ == '__main__':
    t = time.time()
    pool = multiprocessing.Pool(processes = None) ## doesn't matter how many I use
    r = pool.map_async(AddOne, (i for i in vec))
    print("Multiprocessing time (async) -- all threads")
    newVec = r.get()
    print(time.time() - t)
#    print(r.get())
