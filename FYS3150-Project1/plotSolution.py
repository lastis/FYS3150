from scitools.std import *
from numpy import *
from StringIO import StringIO

sol1 = loadtxt('N10.txt')
sol2 = loadtxt('N100.txt')
sol3 = loadtxt('N1000.txt')

X = linspace(0, 1, 1000)
U = 1-(1-exp(-10))*X-exp(-10*X)

plot(sol1[:,0], sol1[:,1])
legend('numerical N=10')
hold('on')
plot(sol2[:,0], sol2[:,1])
legend('numerical N=100')
hold('on')
plot(sol3[:,0], sol3[:,1])
legend('numerical N=1000')
hold('on')
plot(X, U)
legend('exact')

hardcopy('graph.png')
