import CTHGeomCalc

import numpy as np

geomSet=0
run=[]
event=[]
time=[]
drawCdcHits=True
tStart=-10
tEnd=5
fPath="../data/cth48/test.merged."
evList="../data/res_old/res_cth48_thr650.txt"
eThr=0.65 #MeV

def setGeom(gSet):
    global geomSet
    geomSet = gSet

def setList(dataPath,listName):
    global fPath
    global run
    global event
    global time
    fPath=dataPath
    run,event,time=np.loadtxt(listName,unpack=True)
    #print('Dump list:', len(run))
    #for i in range(len(run)):
    #    print(int(run[i]),', ',int(event[i]),', ',int(time[i]))

def setThre(thr):
    global eThr
    eThr=thr

def main():
    if   geomSet==0: ## default 48
        CTHGeomCalc.SetNumSeg(48)
        CTHGeomCalc.SetTiltAngle(20)
        CTHGeomCalc.SetParameters(410,450,55,55,5,10) #Position and size parameters
    elif geomSet==1: ## 64, 5mmScint+10mmCheren
        CTHGeomCalc.SetNumSeg(64)
        CTHGeomCalc.SetTiltAngle(20)
        #CTHGeomCalc.SetParameters(420,460,44,49,7.5,7.5) #Position and size parameters
        CTHGeomCalc.SetParameters(419.5,459.4,42.6,46.8,5,10) #Position and size parameters
    elif geomSet==2: ## 64, 5mmScint+5mmCheren
        CTHGeomCalc.SetNumSeg(64)
        CTHGeomCalc.SetTiltAngle(20)
        CTHGeomCalc.SetParameters(420,460,44,49,5,5) #Position and size parameters
    elif geomSet==3: ## 64, 5mmScint+5mmCheren
        CTHGeomCalc.SetNumSeg(64)
        CTHGeomCalc.SetTiltAngle(17)
        CTHGeomCalc.SetParameters(420,460,42,46,5,5) #Position and size parameters
    CTHGeomCalc.Init(drawCdcHits)
    nTrig=[0,0]
    for i in range(len(time)):
        tRange=[time[i]+tStart,time[i]+tEnd]
        fName=fPath+str(int(run[i])).zfill(3)+".root"
        print("Read :",fName, " , event=",int(event[i]))
        print("   From",tRange[0]," to ",tRange[1]," ns")
        hits=CTHGeomCalc.GetHitsFromTestData(fName,int(event[i]),tRange[0],tRange[1],eThr)
        CTHGeomCalc.AddHits(hits)
        CTHGeomCalc.DrawGeom(True,False)
        drawUS, drawDS = CTHGeomCalc.GetDrawStatus()
        if drawUS==True: nTrig[0]=nTrig[0]+1
        if drawDS==True: nTrig[1]=nTrig[1]+1
    print('Number of trigger found in US/DS: ',nTrig[0],"/",nTrig[1])
