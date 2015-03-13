from scitools.std import *
from StringIO import StringIO

ExpT0_01 = loadtxt("../cpp/2DT0_01Exp.dat")
ExpT0_1 = loadtxt("../cpp/2DT0_1Exp.dat")
ExpT0_5 = loadtxt("../cpp/2DT0_5Exp.dat")


x1 = linspace(0,1,len(ExpT0_01))
x2 = linspace(0,1,len(ExpT0_1))
x3 = linspace(0,1,len(ExpT0_5))

def solution(x,y,t):
	return (1-y)*exp(x+t)

# Difference from T 0.01 to T 0.05
xMat, yMat = ndgrid(x1,x1)
solMat1 = solution(xMat,yMat,0.01)
solMat2 = solution(xMat,yMat,0.5)
diff = abs(solMat2 - solMat1)

pcolor(x1,x1,diff,shading='flat')
colorbar()
axis('equal')
title('Diffusion Difference from T = 0.01 to T = 0.05')
xlabel('x position')
ylabel('y position')
hardcopy('TimeFrom0_01To0_5.png')
hold('off')

# Where T = 0.01
xMat, yMat = ndgrid(x1,x1)
solMat = solution(xMat,yMat,0.01)

pcolor(x1,x1,ExpT0_01,shading='flat')
colorbar()
axis('equal')
title('Diffusion Explicit Scheme 2D T = 0.01')
xlabel('x position')
ylabel('y position')
hardcopy('Explicit2DT0_01.png')
hold('off')

diff = abs(ExpT0_01 - solMat)
pcolor(x1,x1,diff,shading='flat')
colorbar()
axis('equal')
title('Diffusion Error Explicit Scheme 2D T = 0.01')
xlabel('x position')
ylabel('y position')
hardcopy('ErrorExplicit2DT0_01.png')
hold('off')

# Where T = 0.1
xMat, yMat = ndgrid(x2,x2)
solMat = solution(xMat,yMat,0.1)

pcolor(x2,x2,ExpT0_1,shading='flat')
colorbar()
axis('equal')
title('Diffusion Explicit Scheme 2D T = 0.1')
xlabel('x position')
ylabel('y position')
hardcopy('Explicit2DT0_1.png')
hold('off')

diff = abs(ExpT0_1 - solMat)
pcolor(x2,x2,diff,shading='flat')
colorbar()
axis('equal')
title('Diffusion Error Explicit Scheme 2D T = 0.1')
xlabel('x position')
ylabel('y position')
hardcopy('ErrorExplicit2DT0_1.png')
hold('off')

# Where T = 0.5
xMat, yMat = ndgrid(x3,x3)
solMat = solution(xMat,yMat,0.5)

pcolor(x3,x3,ExpT0_5,shading='flat')
colorbar()
axis('equal')
title('Diffusion Explicit Scheme 2D T = 0.5')
xlabel('x position')
ylabel('y position')
hardcopy('Explicit2DT0_5.png')
hold('off')

diff = abs(ExpT0_5 - solMat)
pcolor(x3,x3,diff,shading='flat')
colorbar()
axis('equal')
title('Diffusion Error Explicit Scheme 2D T = 0.5')
xlabel('x position')
ylabel('y position')
hardcopy('ErrorExplicit2DT0_5.png')
hold('off')
