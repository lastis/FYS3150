from scitools.std import *
from StringIO import StringIO

ImpT0_01 = loadtxt("../cpp/2DT0_01Imp.dat")
ImpT0_1 = loadtxt("../cpp/2DT0_1Imp.dat")
ImpT0_5 = loadtxt("../cpp/2DT0_5Imp.dat")


x1 = linspace(0,1,len(ImpT0_01))
x2 = linspace(0,1,len(ImpT0_1))
x3 = linspace(0,1,len(ImpT0_5))

def solution(x,y,t):
	return (1-y)*exp(x+t)

# Where T = 0.01
xMat, yMat = ndgrid(x1,x1)
solMat = solution(xMat,yMat,0.01)

pcolor(x1,x1,solMat,shading='flat')
colorbar()
axis('equal')
title('Diffusion Analytical 2D T = 0.01')
xlabel('x position')
ylabel('y position')
hardcopy('Analytical2DT0_01.png')
hold('off')

pcolor(x1,x1,ImpT0_01,shading='flat')
colorbar()
axis('equal')
title('Diffusion Implicit Scheme 2D T = 0.01')
xlabel('x position')
ylabel('y position')
hardcopy('Implicit2DT0_01.png')
hold('off')

diff = abs(ImpT0_01 - solMat)
pcolor(x1,x1,diff,shading='flat')
colorbar()
axis('equal')
title('Diffusion Error Implicit Scheme 2D T = 0.01')
xlabel('x position')
ylabel('y position')
hardcopy('ErrorImplicit2DT0_01.png')
hold('off')

# Where T = 0.1
xMat, yMat = ndgrid(x2,x2)
solMat = solution(xMat,yMat,0.1)

pcolor(x2,x2,solMat,shading='flat')
colorbar()
axis('equal')
title('Diffusion Analytical 2D T = 0.1')
xlabel('x position')
ylabel('y position')
hardcopy('Analytical2DT0_1.png')
hold('off')

pcolor(x2,x2,ImpT0_1,shading='flat')
colorbar()
axis('equal')
title('Diffusion Implicit Scheme 2D T = 0.1')
xlabel('x position')
ylabel('y position')
hardcopy('Implicit2DT0_1.png')
hold('off')

diff = abs(ImpT0_1 - solMat)
pcolor(x2,x2,diff,shading='flat')
colorbar()
axis('equal')
title('Diffusion Error Implicit Scheme 2D T = 0.1')
xlabel('x position')
ylabel('y position')
hardcopy('ErrorImplicit2DT0_1.png')
hold('off')

# Where T = 0.5
xMat, yMat = ndgrid(x3,x3)
solMat = solution(xMat,yMat,0.5)

pcolor(x3,x3,solMat,shading='flat')
colorbar()
axis('equal')
title('Diffusion Analytical 2D T = 0.5')
xlabel('x position')
ylabel('y position')
hardcopy('Analytical2DT0_5.png')
hold('off')

pcolor(x3,x3,ImpT0_5,shading='flat')
colorbar()
axis('equal')
title('Diffusion Implicit Scheme 2D T = 0.5')
xlabel('x position')
ylabel('y position')
hardcopy('Implicit2DT0_5.png')
hold('off')

diff = abs(ImpT0_5 - solMat)
pcolor(x3,x3,diff,shading='flat')
colorbar()
axis('equal')
title('Diffusion Error Implicit Scheme 2D T = 0.5')
xlabel('x position')
ylabel('y position')
hardcopy('ErrorImplicit2DT0_5.png')
hold('off')
