from scitools.std import *
from StringIO import StringIO

graph = loadtxt('../cpp/SolarSystem.txt')

t  = graph[:,0]
x1 = graph[:,1]
y1 = graph[:,2]
x2 = graph[:,3]
y2 = graph[:,4]
x3 = graph[:,5]
y3 = graph[:,6]
x4 = graph[:,7]
y4 = graph[:,8]
x5 = graph[:,9]
y5 = graph[:,10]
x6 = graph[:,11]
y6 = graph[:,12]
x7 = graph[:,13]
y7 = graph[:,14]
x8 = graph[:,15]
y8 = graph[:,16]
x9 = graph[:,17]
y9 = graph[:,18]
x10 = graph[:,19]
y10 = graph[:,20]
x11 = graph[:,21]
y11 = graph[:,22]

plot(x1,y1,legend='Sun')
hold('on')
"""
plot(x2,y2,legend='mercury')
hold('on')
plot(x3,y3,legend='venus')
hold('on')
"""
plot(x4,y4,legend='Earth')
hold('on')
"""
plot(x5,y5,legend='moon')
hold('on')
plot(x6,y6,legend='mars')
hold('on')
"""
plot(x7,y7,legend='jupiter')
hold('on')
"""
plot(x8,y8,legend='saturn')
hold('on')
plot(x9,y9,legend='uranus')
hold('on')
plot(x10,y10,legend='neptune')
hold('on')
plot(x11,y11,legend='pluto')
"""

axis('equal')
#axis([-50,50,-50,50])
xlabel('Distance (AU)')
ylabel('Distance (AU)')

hardcopy('SolarSystem.png')
