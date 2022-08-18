#!/usr/bin/python
#
# By Xiaotong Li
#
# Hogbom CLEAN
#
import numpy
import scipy.io
from mat4py import loadmat

"""
@INPUT:
    psf   : PSF (dirty beam)
    dirty : Dirty image
    cbeam : CLEAN beam
    
@OUTPUT:
    maximum  : Intensity and position of the peaks
    residual : Residual image
    modelnew : Model image
    skymodel : Restored image (the final output we want!)
"""
def normalminmax(inmodel,m,n,maximum):
    ma = inmodel.max()
    outmodel = numpy.zeros((m,n))
    for i in range(m):
        for j in range(n):
            outmodel[i][j] = inmodel[i][j]/ma*maximum
    return outmodel

    
def HogCLEAN(cycle_previous, cycle, xlen, ylen, p, residual, modelnew, sky, maximum):
    for cycle_number in range(cycle):
        for i in range(xlen):
            for j in range(ylen):
                if i==0 and j==0:
                    maximum[cycle_previous + cycle_number][0] = residual[i][j]
                    maximum[cycle_previous + cycle_number][1] = i
                    maximum[cycle_previous + cycle_number][2] = j
                else:
                    if(residual[i][j]>maximum[cycle_previous + cycle_number][0]):
                        maximum[cycle_previous + cycle_number][0] = residual[i][j]
                        maximum[cycle_previous + cycle_number][1] = i
                        maximum[cycle_previous + cycle_number][2] = j
        modelnew[int(maximum[cycle_previous + cycle_number][1])][int(maximum[cycle_previous + cycle_number][2])] = modelnew[int(maximum[cycle_previous + cycle_number][1])][int(maximum[cycle_previous + cycle_number][2])] + 0.1*maximum[cycle_previous + cycle_number][0]
        psfw = numpy.zeros((xlen,ylen))
        fw = numpy.zeros((xlen,ylen))
        for k1 in range(xlen):
            for k2 in range(ylen):
                psfw[k1][k2] = p[int(xlen-maximum[cycle_previous + cycle_number][1]+k1)][int(ylen-maximum[cycle_previous + cycle_number][2]+k2)]
                fw[k1][k2] = f[int(xlen-maximum[cycle_previous + cycle_number][1]+k1)][int(ylen-maximum[cycle_previous + cycle_number][2]+k2)]
        
        model = normalminmax(psfw,xlen,ylen,maximum[cycle_previous + cycle_number][0])
        sk = normalminmax(fw,xlen,ylen,maximum[cycle_previous + cycle_number][0])
        
        for i1 in range(xlen):
            for j1 in range(ylen):
                residual[i1][j1] = residual[i1][j1] - 0.1*model[i1][j1]
        for i2 in range(xlen):
            for j2 in range(ylen):
                sky[i2][j2] = sky[i2][j2] + 0.1*sk[i2][j2]
        
        if maximum[cycle_previous + cycle_number][0] < maximum[0][0]/50:
            print('Stop after',cycle_previous + cycle_number,'cycles.')
            break
    return residual, modelnew, sky, maximum, cycle_number

if __name__ == "__main__":
          
    Reference_Number = '_HogCLEAN_'
    
    print('Start Hogbom CLEAN ...')
    
    Psf1 = loadmat('psf.mat')#2048
    p = Psf1['psf']
    p = numpy.array(p)
	
    d = loadmat('dirty.mat')#1024
    dirty = d['dirty']
    dirty = numpy.array(dirty)
    
    CBeam1 = loadmat('cbeam.mat')#2048
    f = CBeam1['cbeam']
    f = numpy.array(f)
	
    xlen = len(dirty)
    print('xlen',xlen)
    ylen = len(dirty[0])
    print('ylen',ylen)
    
    residual = dirty
    modelnew = numpy.zeros((xlen,ylen))
    sky = numpy.zeros((xlen,ylen))
    cycle = 2000
    maximum = numpy.zeros((cycle,3))
    cycle_previous = 0
    
    residual, modelnew, sky, maximum, cycle_number = HogCLEAN(cycle_previous, cycle, xlen, ylen, p, residual, modelnew, sky, maximum)
    
    skymodel = sky + residual
    
    cycle_previous = cycle_previous + cycle_number
    
    scipy.io.savemat('maximum'+Reference_Number+str(cycle_previous)+'.mat',mdict = {'maximum':maximum})
    scipy.io.savemat('residual'+Reference_Number+str(cycle_previous)+'.mat',mdict={'residual':residual})
    scipy.io.savemat('modelnew'+Reference_Number+str(cycle_previous)+'.mat',mdict={'modelnew':modelnew})
    scipy.io.savemat('skymodel'+Reference_Number+str(cycle_previous)+'.mat',mdict={'skymodel':skymodel})
    
    print('Hogbom CLEAN finished.')
