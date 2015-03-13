from scitools.std import *
from StringIO import StringIO

T0_01 = loadtxt('../cpp/2DT0_01ProbWalk.dat')
T0_1 	= loadtxt('../cpp/2DT0_1ProbWalk.dat')
T0_5 	= loadtxt('../cpp/2DT0_5ProbWalk.dat')

x1,probx1,proby1 = T0_01.T
x2,probx2,proby2 = T0_1.T
x3,probx3,proby3 = T0_5.T

probx1 = probx1/max(probx1)
probx2 = probx2/max(probx2)
probx3 = probx3/max(probx3)
proby1 = proby1/max(proby1)
proby2 = proby2/max(proby2)
proby3 = proby3/max(proby3)

plot(x1,probx1)
hold('on')
plot(x1,proby1)
title('Diffusion Walker 2D T = 0.01')
legend('x axis','y axis')
xlabel('Position')
ylabel('Solution')
hardcopy('Walker2DT0_01.png')
hold('off')

plot(x2,probx2)
hold('on')
plot(x2,proby2)
title('Diffusion Walker 2D T = 0.1')
legend('x axis','y axis')
xlabel('Position')
ylabel('Solution')
hardcopy('Walker2DT0_1.png')
hold('off')

plot(x3,probx3)
hold('on')
plot(x3,proby3)
title('Diffusion Walker 2D T = 0.5')
legend('x axis','y axis')
xlabel('Position')
ylabel('Solution')
hardcopy('Walker2DT0_5.png')
hold('off')
