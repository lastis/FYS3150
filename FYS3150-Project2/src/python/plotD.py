from scitools.std import *
from StringIO import StringIO

graph1 = loadtxt('../cpp/WaveFunc1.txt')
graph2 = loadtxt('../cpp/WaveFunc2.txt')
graph3 = loadtxt('../cpp/WaveFunc3.txt')
graph4 = loadtxt('../cpp/WaveFunc4.txt')
for i in range(0,length(graph1[0,:])):
	graph1[:,1] = graph1[:,1]**2
	graph2[:,1] = graph2[:,1]**2
	graph3[:,1] = graph3[:,1]**2
	graph4[:,1] = graph4[:,1]**2
plot(graph1[:,0], graph1[:,1])
legend('Omega = 0.01')
hold('on')
plot(graph2[:,0], graph2[:,1])
legend('Omega = 0.50')
hold('on')
plot(graph3[:,0], graph3[:,1])
legend('Omega = 1.00')
hold('on')
plot(graph4[:,0], graph4[:,1])
legend('Omega = 5.00')


hardcopy('WaveFunc.png')
