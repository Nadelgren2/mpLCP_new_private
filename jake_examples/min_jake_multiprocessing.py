from cypari2 import Pari
import multiprocessing
import time

def AddOne(v):
    return v + pari.one()

pari = Pari()
vec = []
for i in range(20000):
    vec.append(pari('x_' + str(i+1)))
#print(vec)

#doesn't work
if __name__ == '__main__':
    t = time.time()
    pool = multiprocessing.Pool(processes = None) ## doesn't matter how many I use
    newVec = pool.map(AddOne, (i for i in vec))
    print(time.time() - t)
#    print(newVec)
