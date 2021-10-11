import sympy
import multiprocessing
import time

def Div(num,den):
    return num/den

vec = []
for i in range(2):
    vec.append(sympy.symbols('x_'+str(i)))
#print(vec)

t = time.time()
if __name__ == '__main__':
    pool = multiprocessing.Pool(processes = None) ## doesn't matter how many I use
    newVec = pool.map(Div, ((i, vec[0]) for i in vec))
    print(newVec)
print(time.time() - t)
