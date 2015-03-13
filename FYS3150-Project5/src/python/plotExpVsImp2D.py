from scitools.std import *
from StringIO import StringIO

ExpT0_01 = loadtxt("../cpp/2DT0_01Exp.dat")
ExpT0_1 = loadtxt("../cpp/2DT0_1Exp.dat")
ExpT0_5 = loadtxt("../cpp/2DT0_5Exp.dat")
ImpT0_01 = loadtxt("../cpp/2DT0_01Imp.dat")
ImpT0_1 = loadtxt("../cpp/2DT0_1Imp.dat")
ImpT0_5 = loadtxt("../cpp/2DT0_5Imp.dat")


x1 = linspace(0,1,len(ExpT0_01))
x2 = linspace(0,1,len(ExpT0_1))
x3 = linspace(0,1,len(ExpT0_5))
x4 = linspace(0,1,len(ExpT0_01))
x5 = linspace(0,1,len(ExpT0_1))
x6 = linspace(0,1,len(ExpT0_5))

def solution(x,y,t):
	return (1-y)*exp(x+t)


# Where T = 0.01
xMat, yMat = ndgrid(x1,x1)
solMat = solution(xMat,yMat,0.01)

diffExp = abs(ExpT0_01 - solMat)
diffImp = abs(ImpT0_01 - solMat)
diff = diffExp - diffImp
pcolor(x1,x1,diff,shading='flat')
colorbar()
axis('equal')
title('Error Explicit vs Implicit 2D T = 0.01')
xlabel('x position')
ylabel('y position')
hardcopy('ErrorExpVsImp2DT0_01.png')
hold('off')

# Where T = 0.1
xMat, yMat = ndgrid(x2,x2)
solMat = solution(xMat,yMat,0.1)

diffExp = abs(ExpT0_1 - solMat)
diffImp = abs(ImpT0_1 - solMat)
diff = diffExp - diffImp
pcolor(x2,x2,diff,shading='flat')
colorbar()
axis('equal')
title('Error Explicit vs Implicit 2D T = 0.1')
xlabel('x position')
ylabel('y position')
hardcopy('ErrorExpVsImp2DT0_1.png')
hold('off')

# Where T = 0.5
xMat, yMat = ndgrid(x3,x3)
solMat = solution(xMat,yMat,0.5)

diffExp = abs(ExpT0_5 - solMat)
diffImp = abs(ImpT0_5 - solMat)
diff = diffExp - diffImp
pcolor(x3,x3,diff,shading='flat')
colorbar()
axis('equal')
title('Error Explicit vs Implicit 2D T = 0.5')
xlabel('x position')
ylabel('y position')
hardcopy('ErrorExpVsImp2DT0_5.png')
hold('off')
