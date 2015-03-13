from scitools.std import *
from StringIO import StringIO

graph = loadtxt('../cpp/test.txt')

x = graph[:,0]
u1 = graph[:,1]
u2 = graph[:,2]
plot(x,u1,'-r')
hold('on')
plot(x,u2,'-b')
