from __future__ import division, print_function
import ctypes
from visual import *
import math

simDLL = ctypes.CDLL("./BracketingLambdaFcnOfS.dll")
simDLL.getBx.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double)
simDLL.getBx.restype = ctypes.c_double
simDLL.getBy.argtypes = (ctypes.c_double, ctypes.c_double, ctypes.c_double)
simDLL.getBy.restype = ctypes.c_double
simDLL.getSAtLambda.argtypes = (ctypes.c_double, ctypes.c_double)
simDLL.getSAtLambda.restype = ctypes.c_double

def drawLine(windObj, p, laxis, headwd=0.005, headln=0.001, shaftwd=0.01,
    col=color.red):
    px = p[0]
    py = p[1]
    pz = p[2]
    axis_x = laxis[0]
    axis_y = laxis[1]
    axis_z = laxis[2]
    windObj.select()
    arrow(pos=(px, py, pz), axis=(axis_x, axis_y, axis_z), color=col, 
        headwidth=headwd, headlength=headln, shaftwidth=shaftwd)

def setupDraw():
	windObj = display(title='Geomagnetic Field Line', autocenter=0, width=1920, 
        height=1080, center=(0.309016994,0.951056516,0), exit=0, range=(15,15,15))
	sphObj = sphere( pos=(0,0,0), color=color.green, opacity=0.35, radius=1 )
	
	#for i in [[10,0,0], [0,10,0], [0,0,10]]:
		#drawLine(windObj, [0,0,0], i, shaftwd=0.05, headwd=0.1, headln=0.15,
            #col=color.cyan)
	xaxpt=[0,1,2,3,4,5,6,7,8,9,10]
	yaxpt=[0,1,2,3,4,5,6,7,8,9,10]
	zaxpt=[0,1,2,3,4,5,6,7,8,9,10]
	xlbl = label(pos=(10,1,0), text='x')
	ylbl =  label(pos=(1,10,0), text='y')
	zlbl =  label(pos=(0,1,10), text='z')
	for i in range(0,10,1):
		points(pos=(xaxpt[i],0,0), size=5, color=color.cyan)
		points(pos=(0,yaxpt[i],0), size=5, color=color.cyan)
		points(pos=(0,0,zaxpt[i]), size=5, color=color.cyan)
	
	return windObj, sphObj

def drawBlines(windObj, p, s_max, L, ILATdegrees, pupbound=[None,None,None], plobound=[None,None,None], 
    numiter=None, linelength=None, multlng=None):
    """Draw B field lines starting at po and ending at ####."""
    loopind = 0
    BoundBool = True
    pts=[]
    pts.append((p[0],p[1],p[2]))
    s = 0
    while BoundBool:
        #print("Hi")
        if numiter is not None:
            loopind += 1
            BoundBool = BoundBool and loopind < numiter
        
        #print("Hi 2")
        
        Bxyz = [simDLL.getBx(s, L, ILATdegrees), simDLL.getBy(s, L, ILATdegrees), 0.0]
        if loopind % 100000 == 0:
            print(loopind, s, math.atan2(p[1], p[0]) * 180.0 / 3.14159265359, Bxyz, p)
        if multlng is not None:
            Bxyz[0] *= multlng; Bxyz[1] *= multlng; Bxyz[2] *= multlng
        #print("Hi 3")
        s += math.sqrt(Bxyz[0]**2 + Bxyz[1]**2) * 6.371e6
        #print("s: ", s)
        p[0] += Bxyz[0]; p[1] += Bxyz[1]; p[2] += Bxyz[2]
        pts.append((p[0],p[1],p[2]))
    windObj.select()
    curve(pos=pts, color=color.yellow)

def drawFieldLine(ILATdegrees, lambdaBins):
    windObj, earthPic = setupDraw()
    windObj.select()
    
    dlambda = ILATdegrees / (lambdaBins)
    lambda_deg = ILATdegrees
    L = 6.371e6/math.cos(ILATdegrees * 3.14159265359 / 180.0)**2
    s_max = simDLL.getSAtLambda(ILATdegrees, L)

    for iii in range(int(lambdaBins + 1)):
        lambda_rad = lambda_deg * 3.14159265359 / 180.0
        xy = [ L/6.371e6 * math.cos(lambda_rad)**3, L/6.371e6 * math.cos(lambda_rad)**2 * math.sin(lambda_rad) ]
        points(pos=(xy[0], xy[1], 0), size=5, color=color.red)
        s = s_max - simDLL.getSAtLambda(math.atan2(xy[1], xy[0]) * 180.0 / 3.14159265359, L)
        Bxyz = [simDLL.getBx(s, L, ILATdegrees), simDLL.getBy(s, L, ILATdegrees), 0.0]
        print(lambda_deg, [Bxyz[0],Bxyz[1]], math.sqrt(Bxyz[0]**2 + Bxyz[1]**2), s)
        lambda_deg -= dlambda
    
    #xyz = [ L / 6.371e6 * math.cos(ILATdegrees * 3.14159265359 / 180)**3, L / 6.371e6 * math.cos(ILATdegrees * 3.14159265359 / 180)**2 * math.sin(ILATdegrees * 3.14159265359 / 180), 0 ]
    #drawBlines(windObj, xyz, s_max, L, ILATdegrees, numiter=1000000, multlng=10)
    print("Done")


drawFieldLine(72.0, 20.0)
while (True):
	rate(30)
