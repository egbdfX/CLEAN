#!/usr/bin/python
#
# Multi-Scale CLEAN (MS CLEAN)
#
import numpy
import scipy.io
from mat4py import loadmat

"""
@INPUT:
    res_scalestack      : Scaled residual images
    psf_scalescalestack : Scaled and crossed PSFs
    cbeam               : CLEAN beam
    
@OUTPUT:
    maximum  : Intensity and position of the peaks
    residual : Residual images
    modelnew : Model image
    skymodel : Restored image (the final output we want!)
"""

def normalmax(inmodel,m,n,maximum):
    ma = inmodel.max()
    outmodel = numpy.zeros((m,n))
    for i in range(m):
        for j in range(n):
            outmodel[i][j] = inmodel[i][j]/ma*maximum
    return outmodel
    
def MSCLEAN(cyc_pre, cycle, n_scale, xlen, ylen, residual, modelnew, psfun, psfcross, maximum):
    sky = numpy.zeros((xlen,ylen))
    for cyc_num in range(cycle):
        
        # Find the maximum peak in n_scale residuals (h: order number of scale, i&j: x&y)
        for h in range(n_scale):
            for i in range(xlen):
                for j in range(ylen):
                    if h==0 and i==0 and j==0:
                        maximum[cyc_pre + cyc_num][0] = residual[h][i][j]
                        maximum[cyc_pre + cyc_num][1] = i
                        maximum[cyc_pre + cyc_num][2] = j
                        maximum[cyc_pre + cyc_num][3] = h
                    else:
                        if(residual[h][i][j]>maximum[cyc_pre + cyc_num][0]):
                            maximum[cyc_pre + cyc_num][0] = residual[h][i][j]
                            maximum[cyc_pre + cyc_num][1] = i
                            maximum[cyc_pre + cyc_num][2] = j
                            maximum[cyc_pre + cyc_num][3] = h
        
        # Calculate model and the corresponding residual
        modelnew[int(maximum[cyc_pre + cyc_num][1])][int(maximum[cyc_pre + cyc_num][2])] = modelnew[int(maximum[cyc_pre + cyc_num][1])][int(maximum[cyc_pre + cyc_num][2])] + 0.1*maximum[cyc_pre + cyc_num][0]
        
        psfw = numpy.zeros((xlen,ylen))
        fw = numpy.zeros((xlen,ylen))
        for k1 in range(xlen):
            for k2 in range(ylen):
                psfw[k1][k2] = psfun[int(maximum[cyc_pre + cyc_num][3])][int(xlen-maximum[cyc_pre + cyc_num][1]+k1)][int(ylen-maximum[cyc_pre + cyc_num][2]+k2)]
                fw[k1][k2] = f[int(xlen-maximum[cyc_pre + cyc_num][1]+k1)][int(ylen-maximum[cyc_pre + cyc_num][2]+k2)]
        
        model = normalmax(psfw,xlen,ylen,maximum[cyc_pre + cyc_num][0])
        sk = normalmax(fw,xlen,ylen,maximum[cyc_pre + cyc_num][0])
        for i1 in range(xlen):
            for j1 in range(ylen):
                residual[int(maximum[cyc_pre + cyc_num][3])][i1][j1] = residual[int(maximum[cyc_pre + cyc_num][3])][i1][j1] - 0.1*model[i1][j1]
        for i2 in range(xlen):
            for j2 in range(ylen):
                sky[i2][j2] = sky[i2][j2] + 0.1*sk[i2][j2]
        
        # Calculate other residuals
        psfcro = numpy.zeros((n_scale,xlen,ylen))
        for h1 in range(n_scale):
            if h1 != maximum[cyc_pre + cyc_num][3]:
                for i in range(xlen):
                    for j in range(ylen):
                        psfcro[h1][i][j] = psfcross[int(maximum[cyc_pre + cyc_num][3])][h1][int(xlen-maximum[cyc_pre + cyc_num][1]+i)][int(ylen-maximum[cyc_pre + cyc_num][2]+j)]
        
        modelcro = numpy.zeros((n_scale,xlen,ylen))
        for h1 in range(n_scale):
            if h1 != maximum[cyc_pre + cyc_num][3]:
                modelcro[h1] = normalmax(psfcro[h1],xlen,ylen,maximum[cyc_pre + cyc_num][0])
        for h1 in range(n_scale):
            if h1 != maximum[cyc_pre + cyc_num][3]:
                for i in range(xlen):
                    for j in range(ylen):
                        residual[h1][i][j] = residual[h1][i][j] - 0.1*modelcro[h1][i][j]
                        
        if maximum[cyc_pre + cyc_num][0] < maximum[0][0]/50:
            print('Stop after',cyc_pre + cyc_num,'cycles.')
            break				
    return residual, modelnew, sky, maximum, cyc_num

if __name__ == "__main__":
          
    Reference_Number = '_MsCLEAN_'
    
    print('Start MS CLEAN ...')
    
    d = scipy.io.loadmat('res_scalestack.mat')
    residual = d['res_scalestack']
    residual = numpy.array(residual)
    
    Psf1 = scipy.io.loadmat('psf_scalescalestack.mat')
    psfcross = Psf1['psf_scalescalestack']
    psfcross = numpy.array(psfcross)
    
    CBeam1 = loadmat('cbeam.mat')
    f = CBeam1['cbeam']
    f = numpy.array(f)
    
    n_scale = 6
    xlen = len(residual[0])
    print('xlen',xlen)
    ylen = len(residual[0][0])
    print('ylen',ylen)
    
    psfun = numpy.zeros((n_scale,2*xlen,2*ylen))
    
    for i in range(n_scale):
        psfun[i] = psfcross[i][i]
    
    cyc_pre = 0
    cycle = 2000
    
    modelnew = numpy.zeros((xlen,ylen))
    maximum = numpy.zeros((cycle,4))
    
    residual, modelnew, sky, maximum, cyc_num = MSCLEAN(cyc_pre, cycle, n_scale, xlen, ylen, residual, modelnew, psfun, psfcross, maximum)
    
    skymodel = sky + residual[int(maximum[0][3])]
    
    cyc_pre = cyc_pre + cyc_num
    
    scipy.io.savemat('residual'+Reference_Number+str(cyc_pre)+'.mat',mdict={'residual':residual})
    scipy.io.savemat('modelnew'+Reference_Number+str(cyc_pre)+'.mat',mdict={'modelnew':modelnew})
    scipy.io.savemat('maximum'+Reference_Number+str(cyc_pre)+'.mat',mdict = {'maximum':maximum})
    scipy.io.savemat('skymodel'+Reference_Number+str(cyc_pre)+'.mat',mdict={'skymodel':skymodel})
    
    print('MS CLEAN finished.')
