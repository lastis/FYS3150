from scitools.std import *
from StringIO import StringIO

T0_01 = loadtxt('../cpp/1DT0_01ProbWalk.dat')
T0_1 	= loadtxt('../cpp/1DT0_1ProbWalk.dat')
T0_5 	= loadtxt('../cpp/1DT0_5ProbWalk.dat')
T0_01G = loadtxt('../cpp/1DT0_01GaussProbWalk.dat')
T0_1G 	= loadtxt('../cpp/1DT0_1GaussProbWalk.dat')
T0_5G 	= loadtxt('../cpp/1DT0_5GaussProbWalk.dat')

x1,probx1,proby1 = T0_01.T
x2,probx2,proby2 = T0_1.T
x3,probx3,proby3 = T0_5.T
x1G,probx1G,proby1G = T0_01G.T
x2G,probx2G,proby2G = T0_1G.T
x3G,probx3G,proby3G = T0_5G.T

probx1 = probx1/max(probx1)
probx2 = probx2/max(probx2)
probx3 = probx3/max(probx3)
probx1G = probx1G/max(probx1G)
probx2G = probx2G/max(probx2G)
probx3G = probx3G/max(probx3G)

def solution(x,t,N):
	pi = 3.14
	num = 0
	sum = 0
	for n in xrange(1,N):
		num = 1.0/n*sin(n*pi*x)*exp(-(n**2)*(pi**2)*t)
		num = 2*num/pi
		sum = sum + num
	return 1-x-sum


plot(x1,probx1)
hold('on')
plot(x1,solution(x1,0.01,100))
title('Diffusion Walker 1D T = 0.01')
legend('Numerical Solution','Analytical Solution')
xlabel('Position')
ylabel('Solution')
hardcopy('Walker1DT0_01.png')
hold('off')

plot(x2,probx2)
hold('on')
plot(x2,solution(x2,0.1,100))
title('Diffusion Walker 1D T = 0.1')
legend('Numerical Solution','Analytical Solution')
xlabel('Position')
ylabel('Solution')
hardcopy('Walker1DT0_1.png')
hold('off')

plot(x3,probx3)
hold('on')
plot(x3,solution(x3,0.5,100))
title('Diffusion Walker 1D T = 0.5')
legend('Numerical Solution','Analytical Solution')
xlabel('Position')
ylabel('Solution')
hardcopy('Walker1DT0_5.png')
hold('off')

# Gauss plots
plot(x1G,probx1G)
hold('on')
plot(x1G,solution(x1G,0.01,100))
title('Diffusion Walker Gaussian Step 1D T = 0.01')
legend('Numerical Solution','Analytical Solution')
xlabel('Position')
ylabel('Solution')
hardcopy('Walker1DGaussT0_01.png')
hold('off')

plot(x2G,probx2G)
hold('on')
plot(x2G,solution(x2G,0.1,100))
title('Diffusion Walker Gaussian Step 1D T = 0.1')
legend('Numerical Solution','Analytical Solution')
xlabel('Position')
ylabel('Solution')
hardcopy('Walker1DGaussT0_1.png')
hold('off')

plot(x3G,probx3G)
hold('on')
plot(x3G,solution(x3G,0.5,100))
title('Diffusion Walker Gaussian Step 1D T = 0.5')
legend('Numerical Solution','Analytical Solution')
xlabel('Position')
ylabel('Solution')
hardcopy('Walker1DGaussT0_5.png')
hold('off')

# Error plots


diff  = (probx1  - solution(x1,0.01,100))
diffG = (probx1G - solution(x1,0.01,100))

plot(x1,diff)
hold('on')
plot(x1,diffG)
title('Diffusion Error Walker 1D T = 0.01')
legend('Constant step','Gaussian step')
xlabel('Position')
ylabel('Solution')
hardcopy('ErrorWalker1DT0_01.png')
hold('off')

diff  = (probx2  - solution(x2,0.1,100))
diffG = (probx2G - solution(x2,0.1,100))

plot(x2,diff)
hold('on')
plot(x2,diffG)
title('Diffusion Error Walker 1D T = 0.1')
legend('Constant step','Gaussian step')
xlabel('Position')
ylabel('Solution')
hardcopy('ErrorWalker1DT0_1.png')
hold('off')

diff  = (probx3  - solution(x3,0.5,100))
diffG = (probx3G - solution(x3,0.5,100))

plot(x3,diff)
hold('on')
plot(x3,diffG)
title('Diffusion Error Walker 1D T = 0.5')
legend('Constant step','Gaussian step')
xlabel('Position')
ylabel('Solution')
hardcopy('ErrorWalker1DT0_5.png')
hold('off')
