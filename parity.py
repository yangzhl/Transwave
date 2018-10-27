#!/usr/bin/python
# coding: utf-8

import subprocess

## create GCOEFF.txt file 
file=open("GCOEFF.txt",'w')
file.close()

def parity(soc,kpoint,nbandmax):
    ## convert to int 
    nbandmax=int(nbandmax)
    ## 如果具有inversion symmetry 存在两重简并 step=2
    step=1
    if (soc == 't'):
        step=2
    ## 读取k点处的波函数信息，输出至GCOEFF.txt文件
    for nband in range(1,nbandmax+1,step):
        ## subprocess.call paramater must be string type
        subprocess.call(["Transwave","-soc",soc,"-k",kpoint,"-b",str(nband)])

    coefile=open("GCOEFF.txt",'r') 
    lines=coefile.readlines()
    coefile.close()
    ## 布里渊区高对称点k处 不同能带nband的宇称
    parlist=[]
    ## loop over each band
    for i in range(0,len(lines)):
        line=lines[i]
        if 'nplane' in line:
            shift=[0,0,0]
            igtmp=[0,0,0]
            coefftmp=[0,0]
            kposition=lines[i+2].split()
            kx=float(kposition[0])
            ky=float(kposition[1])
            kz=float(kposition[2])
            ## k+G  和k+G' G'=G+shift-vector 
            if (kx >= 1.0e-5):
                shift[0]=-1
            elif (kx <= 1.0e-5 and kx >= -1.0e-5):
                shift[0]=0
            elif (kx <= -1.0e-5):
                shift[0]=1
            if (ky >= 1.0e-5):
                shift[1]=-1
            elif (ky <= 1.0e-5 and ky >=-1.0e-5):
                shift[1]=0
            elif (ky <= -1.0e-5):
                shift[1]=1
            if (kz >= 1.0e-5 ):
                shift[2]=-1
            elif (kz <= 1.0e-5 and kz >=-1.0e-5):
                shift[2]=0
            elif (kz <=-1.0e-5):
                shift[2]=1
            #print 'shift',shift
            ## 首先找到合适的平面波 系数
            for j in range(4,100):
                coeff=lines[i+j].split()
                ig1=int(coeff[0])
                ig2=int(coeff[1])
                ig3=int(coeff[2])
                ## real part
                coeffr=float(coeff[4])
                ## image part
                coeffi=float(coeff[6])
                if(ig1==0 and ig2==0 and ig3==0):
                    continue
                elif( abs(coeffr)<=1.0e-4 or abs(coeffi)<=1.0e-4 ):
                    continue
                else:
                    ## igtmp 为G',coefftmp为k+G的系数 
                    print 'ig',ig1,ig2,ig3
                    igtmp[0]=-(ig1)+shift[0]
                    igtmp[1]=-(ig2)+shift[1]
                    igtmp[2]=-(ig3)+shift[2]
                    print 'igtmp',igtmp[0],igtmp[1],igtmp[2]
                    coefftmp[0]=coeffr
                    coefftmp[1]=coeffi
                    break
            #print 'coefftmp',coefftmp
            for m in range(1,500):
                coeff=lines[i+j+m].split()
                ig1=int(coeff[0])
                ig2=int(coeff[1])
                ig3=int(coeff[2])
                coeffr=float(coeff[4])
                coeffi=float(coeff[6])

                if(ig1==igtmp[0] and ig2==igtmp[1] and ig3==igtmp[2]):
                    ##找到了k+G',对比于k+G的系数
                    difcoeffr=abs(abs(coeffr)-abs(coefftmp[0]))
                    difcoeffi=abs(abs(coeffi)-abs(coefftmp[1]))
                    print 'diff',difcoeffr,difcoeffi
                    if (abs(abs(coeffr)-abs(coefftmp[0]))<=1.0e-1 and abs(abs(coeffi)-abs(coefftmp[1]))<=1.0e-1):
                        sign1=sign(coeffr)
                        sign2=sign(coefftmp[0])
                        sign3=sign(coeffi)
                        sign4=sign(coefftmp[1])
                        print sign1,sign2,sign3,sign4
                        ## 同号parity=1
                        if (sign1*sign2==1 and sign3*sign4==1):
                            parlist.append(1)
                        elif(sign1*sign2==-1 and sign3*sign4==-1):
                            parlist.append(-1)
                        else:
                            print "program error"
                        break
                else:
                    continue
    print 'the parity list is '
    print '('
    kparity=1
    for item in parlist:
        kparity=kparity*item
        if (item == 1 ):
            p='+'
        elif (item == -1):
            p='-'
        print p+','
    print ')'
    print 'the parity in high symmetry k point','(',kx,ky,kz,')','is', kparity
    #return kparity

##判断浮点数的正负号
def sign(x):
    x=float(x)
    if (x>=1.0e-10):
        signal=1
    elif (x<=1.0e-10 and x>=-1.0e-10):
        signal=0
    else:
        signal=-1
    return signal

def Z2():
    print ''
    

if __name__ == "__main__":
    ## soc kpoint nbandmax are string type
    soc=raw_input("weather the system contains soc, t or f ? \n")
    kpoint=raw_input("which kpoint would you like to calculate parity ? \n")
    nbandmax=raw_input("how many energy bands bellow fermi energy ? \n")
    parity(soc,kpoint,nbandmax)
    

