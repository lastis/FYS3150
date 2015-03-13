from scitools.std import *
from StringIO import StringIO

graph = loadtxt('../cpp/test.txt')

plot(graph[:,0], graph[:,1])

hardcopy('test.png')
