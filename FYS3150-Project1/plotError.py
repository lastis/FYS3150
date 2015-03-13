from scitools.std import *
from numpy import *
from StringIO import StringIO

graph1 = loadtxt('errorN10.txt')
graph2 = loadtxt('errorN100.txt')
graph3 = loadtxt('errorN1000.txt')
graph4 = loadtxt('errorN10000.txt')
graph5 = loadtxt('errorN100000.txt')

plot(graph1[:,0], graph1[:,1])
legend('Relative error N=10')
hold('on')
plot(graph2[:,0], graph2[:,1])
legend('Relative error N=100')
hold('on')
plot(graph3[:,0], graph3[:,1])
legend('Relative error N=1000')
hold('on')
plot(graph4[:,0], graph4[:,1])
legend('Relative error N=10000')
hold('on')
plot(graph5[:,0], graph5[:,1])
legend('Relative error N=100000')

hardcopy('errorgraph.png')
