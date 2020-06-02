import uproot

import numpy as np

import matplotlib.pyplot as plt

import math

from matplotlib.collections import LineCollection

kBirks=0.13 # [mm/MeV] Assuming density = 1g/cm^3, https://en.wikipedia.org/wiki/Birks%27_law
####################
# All unit in [mm]
####################
# CDC inner wall
CDCInnerR=495.8
# CTH outer support wall
CTHSupportOuterR=CDCInnerR-5.0
CTHSupportThick=3.0
# CTH counter's possible outer most radius
CTHCounterMaxR=CTHSupportOuterR-CTHSupportThick-3.0 # 5mm clearance from the support wall

# CTH inner support
CTHSupportInnerR=387.0

# CTH counter's central R
LayerGapRadius=40.0 # [D=40.0]
CTHCounterCentROuter=455.0 # [D=450.0]
CTHCounterCentRInner=CTHCounterCentROuter-LayerGapRadius

######################
nSeg=64 # [D=48]
nHod=2 # up/down stream
nType=2 # Cheren/scint
nChannels=nHod*nType*nSeg

tiltAngle=20.0 # definition in ICEDUST [D=20.0]
TiltAngle=math.radians(90.0-tiltAngle)

CentAngle=math.radians(360.0/float(nSeg))

HalfWidthOuter=47.5 # [D=55]
HalfWidthInner=42.5 # [D=55]

ThickOuter=10.0 # [D=10]
ThickInner= 5.0 # [D=5]
addHits=False

def SetNumSeg(nseg):
    global nSeg, CentAngle
    nSeg=nseg
    CentAngle=math.radians(360.0/float(nSeg))
    nChannels=nHod*nType*nSeg

def SetTiltAngle(angle):
    global tiltAngle
    global TiltAngle
    tiltAngle=angle
    TiltAngle=math.radians(90.0-tiltAngle)

def SetParameters(rin,rout,hwin,hwout,thin,thout):
    global CTHCounterCentROuter
    global CTHCounterCentRInner
    global HalfWidthOuter
    global HalfWidthInner
    global ThickOuter
    global ThickInner
    CTHCounterCentROuter=rout
    CTHCounterCentRInner=rin
    HalfWidthOuter=hwout
    HalfWidthInner=hwin
    ThickOuter=thout
    ThickInner=thin
    
### outer layer
centPosXOuter=[]
centPosYOuter=[]
### inner layer
centPosXInner=[]
centPosYInner=[]

lc1=[[],[]]
lc2=[[],[]]
lines1=[]
lines2=[]
lines3=[]
#lc3=[[],[]]

drawCdcHits=False
drawUS=False
drawDS=False
def Init(drawCdc):
    print("#of segments: "+str(nSeg))
    print("Tilt angle: "+str(tiltAngle))
    print("Inner/Outer radial pos: "+str(CTHCounterCentRInner)+" / "+str(CTHCounterCentROuter))
    print("Inner/Outer width:      "+str(2*HalfWidthInner)+" / "+str(2*HalfWidthOuter))
    print("Inner/Outer thickness:  "+str(ThickInner)+" / "+str(ThickOuter))
    global centPosXOuter, centPosXInner
    global centPosYOuter, centPosYInner
    global lines1,lines2#,lines3
    global drawCdcHits
    drawCdcHits=drawCdc
    centPosXOuter=[]
    centPosYOuter=[]
    centPosXInner=[]
    centPosYInner=[]
    lines1=[]
    lines2=[]
    lines3=[]
    for i in range(nSeg):
        ## Outer layer
        x = CTHCounterCentROuter*math.cos((i+0.5)*CentAngle)
        y = CTHCounterCentROuter*math.sin((i+0.5)*CentAngle)
        v = np.array([x,y])
        c, s = np.cos(TiltAngle), np.sin(TiltAngle)
        dirV = np.array([HalfWidthOuter*x/CTHCounterCentROuter,HalfWidthOuter*y/CTHCounterCentROuter])
        newV = np.array([(c * dirV[0] - s * dirV[1]), (s * dirV[0] + c * dirV[1])])
        centPosXOuter.append(v[0])
        centPosYOuter.append(v[1])
        lines1.append([(v[0]+newV[0],v[1]+newV[1]), (v[0]-newV[0],v[1]-newV[1])])
        #### Calculate distance between two counters
        ## Inner layer
        x = CTHCounterCentRInner*math.cos((i+0.5)*CentAngle)
        y = CTHCounterCentRInner*math.sin((i+0.5)*CentAngle)
        v = np.array([x,y])
        #c, s = np.cos(TiltAngle), np.sin(TiltAngle)
        dirV = np.array([HalfWidthInner*x/CTHCounterCentRInner,HalfWidthInner*y/CTHCounterCentRInner])
        newV = np.array([(c * dirV[0] - s * dirV[1]), (s * dirV[0] + c * dirV[1])])
        centPosXInner.append(v[0])
        centPosYInner.append(v[1])
        lines2.append([(v[0]+newV[0],v[1]+newV[1]), (v[0]-newV[0],v[1]-newV[1])])
        ## Radial lines
        lines3.append([(0,0),(CTHCounterMaxR*math.cos((i+0.5)*CentAngle),CTHCounterMaxR*math.sin((i+0.5)*CentAngle))])

def ch_to_info(ch):
    ihod=int(ch/(nType*nSeg))
    itype=int((ch-ihod*nType*nSeg)/nSeg) # cheren/scint=0/1
    iseg=ch%nSeg
    return ihod,itype,iseg

def info_to_ch(ihod,itype,iseg):
    if(iseg<0):
        iseg = iseg+nSeg
    if(iseg>=nSeg):
        iseg = iseg-nSeg
    ch=nType*nSeg*ihod + nSeg*itype + iseg
    return ch

hitInfos=[]
lHitCounterUS=[[],[]] #inner/outer
lHitCounterDS=[[],[]] #inner/outer
nHitCounterUS=[[],[]] #inner/outer
nHitCounterDS=[[],[]] #inner/outer
pArrows=[[],[]]
def AddHits(hitInfoVect): ### ihod, ilayer, iseg, xyz, pxyz, edep, time
    global hitInfos
    global lHitCounterUS
    global lHitCounterDS
    global nHitCounterUS
    global nHitCounterDS
    global pArrows
    ### Init
    hitInfos=[]
    ltemp0=[]
    ltemp1=[]
    ltemp2=[]
    ltemp3=[]
    ptemp=[[],[]]
    if (len(hitInfoVect)): addHits=True
    for ihit in range(len(hitInfoVect)):
        #print(hitInfoVect[ihit])
        hitInfos.append(hitInfoVect[ihit])
        r = CTHCounterCentROuter
        hw= HalfWidthOuter
        if (hitInfos[ihit][1]==1):
            r = CTHCounterCentRInner
            hw= HalfWidthInner
        iseg = segId_to_phiIdx(int(hitInfos[ihit][0]),int(hitInfos[ihit][2]))
        x = r*math.cos((iseg+0.5)*CentAngle)
        y = r*math.sin((iseg+0.5)*CentAngle)
        v = np.array([x,y])
        c, s = np.cos(TiltAngle), np.sin(TiltAngle)
        dirV = np.array([hw*x/r,hw*y/r])
        newV = np.array([(c * dirV[0] - s * dirV[1]), (s * dirV[0] + c * dirV[1])])
        k=2*int(hitInfos[ihit][0])+int(hitInfos[ihit][1])
        #print("### ",x,", ",y)
        #print("$$$ ",hitInfos[ihit][3],", ",hitInfos[ihit][4])
        ptemp[int(k/2)].append([(hitInfos[ihit][3],hitInfos[ihit][4]),
                                (hitInfos[ihit][3]+hitInfos[ihit][6]/2.,hitInfos[ihit][4]+hitInfos[ihit][7]/2.)])
        if   (k==0): ltemp0.append([(v[0]+newV[0],v[1]+newV[1]), (v[0]-newV[0],v[1]-newV[1])])
        elif (k==1): ltemp1.append([(v[0]+newV[0],v[1]+newV[1]), (v[0]-newV[0],v[1]-newV[1])])
        elif (k==2): ltemp2.append([(v[0]+newV[0],v[1]+newV[1]), (v[0]-newV[0],v[1]-newV[1])])
        elif (k==3): ltemp3.append([(v[0]+newV[0],v[1]+newV[1]), (v[0]-newV[0],v[1]-newV[1])])
    lHitCounterUS=[LineCollection(ltemp0, lw=ThickOuter, alpha=1.0, color='red'),
                   LineCollection(ltemp1, lw=ThickInner, alpha=1.0, color='red')]
    lHitCounterDS=[LineCollection(ltemp2, lw=ThickOuter, alpha=1.0, color='red'),
                   LineCollection(ltemp3, lw=ThickInner, alpha=1.0, color='red')]
    nHitCounterUS=[len(ltemp0),len(ltemp1)]
    nHitCounterDS=[len(ltemp2),len(ltemp3)]
    print(nHitCounterUS[0],",",nHitCounterUS[1],",",nHitCounterDS[0],",",nHitCounterDS[1])
    pArrows=[LineCollection(ptemp[0], lw=2, color='blue'),
             LineCollection(ptemp[1], lw=2, color='blue')]

def segId_to_phiIdx(hod,seg):
    phiIdx=-1
    if hod==1:
        phiIdx=(nSeg-seg+int(nSeg/2)-1)%nSeg
    else:
        phiIdx=seg
    return phiIdx

cdcHits=[[],[]]
def GetHitsFromTestData(fileName,evId,tstart,tend,ethr):
    debug=False
    events = uproot.open(fileName)["mc"]
    hitInfos=[]
    hitPhi=[[],[]]
    hitSeg=[[],[]]
    hitXY=[[],[]] # dont care up/down
    global cdcHits
    cdcHits=[[],[]]
    for data in events.iterate(["eventId","cth_*","cdc_x","cdc_y","cdc_z","cdc_t"],entrysteps=1,namedecode="utf-8"):
        eventId=int(data["eventId"])
        if (evId>-1) and (eventId!=evId): continue
        cth_edep=data["cth_edep"]
        cth_stepL=data["cth_stepL"]
        cth_pdgId=data["cth_pdgId"]
        cth_time=data["cth_t"]
        cth_x,cth_y,cth_z=data["cth_x"],data["cth_y"],data["cth_z"]
        cth_px,cth_py,cth_pz=data["cth_px"],data["cth_py"],data["cth_pz"]
        cth_beta=data["cth_beta"]
        cth_segId=data["cth_segId"]
        cth_hodId=data["cth_hodId"] # 0/1 Up/downstream
        cth_cryId=data["cth_cryType"] # 0: Cherenkov, 1: Scintillator, 10: CherenLG, 11: ScintLG
        cth_edep_corr=(data["cth_edep"]/(1.0+kBirks*data["cth_edep"]/data["cth_stepL"]))
        nhits=len(cth_edep[0])
        for ihit in range(nhits):
            if cth_time[0][ihit]<=tstart or cth_time[0][ihit]>tend:continue
            if cth_edep_corr[0][ihit]<ethr: continue
            if cth_cryId[0][ihit]>1: continue
            #if abs(cth_pdgId[0][ihit])!=11: continue
            #if math.sqrt(cth_px[0][ihit]*cth_px[0][ihit]+cth_py[0][ihit]*cth_py[0][ihit]+cth_pz[0][ihit]*cth_pz[0][ihit])<5:continue
            vtemp=np.array([int(cth_hodId[0][ihit]),int(cth_cryId[0][ihit]%10),int(cth_segId[0][ihit]),
                            7650.0-cth_z[0][ihit],cth_y[0][ihit],cth_x[0][ihit],-cth_pz[0][ihit],cth_py[0][ihit],cth_px[0][ihit],
                            cth_edep_corr[0][ihit],cth_time[0][ihit]])
            #print(cth_time[0][ihit],", ",cth_edep[0][ihit]*1e3)
            hitInfos.append(vtemp)
            if (debug==True): 
                hitPhi[int(cth_hodId[0][ihit])].append(math.atan2(cth_y[0][ihit],7650.0-cth_z[0][ihit]))
                hitSeg[int(cth_hodId[0][ihit])].append(segId_to_phiIdx(int(cth_hodId[0][ihit]),int(cth_segId[0][ihit])))
                hitXY[0].append(7650.0-cth_z[0][ihit])
                hitXY[1].append(cth_y[0][ihit])
        if drawCdcHits==True:
            cdc_x=data["cdc_x"]
            cdc_y=data["cdc_y"]
            cdc_z=data["cdc_z"]
            cdc_t=data["cdc_t"]
            nhits=len(cdc_t[0])
            for ihit in range(nhits):
                if cdc_t[0][ihit]<=tstart or cdc_t[0][ihit]>tend:continue
                cdcHits[0].append(cdc_x[0][ihit])
                cdcHits[1].append(cdc_y[0][ihit])
                #print(cdc_x[0][ihit],", ",cdc_y[0][ihit],", ",cdc_z[0][ihit])
   
    if (debug==True): 
        plt.figure(figsize=(10, 10))
        plt.xlim(-1,65)
        plt.ylim(-1.2*math.pi,1.2*math.pi)
        plt.scatter(hitSeg[0],hitPhi[0],color='blue',s=8)
        plt.scatter(hitSeg[1],hitPhi[1],color='red',s=8)
        plt.show()
        plt.figure(figsize=(10, 10))
        plt.scatter(hitXY[0],hitXY[1],color='blue',s=8)
        plt.show()
    return hitInfos

def GetDrawStatus():
    return drawUS, drawDS

def DrawGeom(drawHits):
    global drawUS, drawDS
    drawUS=False
    drawDS=False
    #### All counters
    lc1 = [LineCollection(lines1, lw=ThickOuter, alpha=0.6),
           LineCollection(lines1, lw=ThickOuter, alpha=0.6)]
    lc2 = [LineCollection(lines2, lw=ThickInner, alpha=0.6, color='green'),
           LineCollection(lines2, lw=ThickInner, alpha=0.6, color='green')]
    #lc3 = [LineCollection(lines3, lw=0.5, color='black', ls='--'),
    #       LineCollection(lines3, lw=0.5, color='black', ls='--')]
    if (drawHits==False or (nHitCounterUS[0]>1 and nHitCounterUS[1]>1)):
        drawUS=True
        print("##")
        plt.figure(figsize=(16.0, 16.0))
        plt.title('Upstream',size=22)
        plt.xlim(-CTHCounterMaxR-50,CTHCounterMaxR+50)
        plt.ylim(-CTHCounterMaxR-50,CTHCounterMaxR+50)
        #plt.scatter(centPosXOuter,centPosYOuter,color='blue',linewidth=1.0)
        #plt.scatter(centPosXInner,centPosYInner,color='blue',linewidth=1.0)
        #circle1=plt.Circle((0,0),CTHCounterCentROuter,color='black',fill=False,ls='--')
        #circle2=plt.Circle((0,0),CTHCounterCentRInner,color='black',fill=False,ls='--')
        circle3=plt.Circle((0,0),CTHCounterMaxR,color='black',fill=False,lw=2.0)
        circle4=plt.Circle((0,0),CTHSupportInnerR,color='black',fill=False,lw=2.0) 
        #plt.gcf().gca().add_artist(circle1)
        #plt.gcf().gca().add_artist(circle2)
        plt.gcf().gca().add_artist(circle3)
        plt.gcf().gca().add_artist(circle4)
        plt.gca().add_collection(lc1[0])
        plt.gca().add_collection(lc2[0])
        #plt.gca().add_collection(lc3[0])
        if(drawHits==True):
            plt.gca().add_collection(lHitCounterUS[0])
            plt.gca().add_collection(lHitCounterUS[1])
            plt.gca().add_collection(pArrows[0])
        if(drawCdcHits==True):
            plt.scatter(cdcHits[0],cdcHits[1],s=10,color='red')
        plt.show()

    if (drawHits==False or (nHitCounterDS[0]>1 and nHitCounterDS[1]>1)):
        drawDS=True
        print("##")
        plt.figure(figsize=(16.0, 16.0))
        plt.title('Downstream',size=22)
        plt.xlim(-CTHCounterMaxR-50,CTHCounterMaxR+50)
        plt.ylim(-CTHCounterMaxR-50,CTHCounterMaxR+50)
        #plt.scatter(centPosXOuter,centPosYOuter,color='blue',linewidth=1.0)
        #plt.scatter(centPosXInner,centPosYInner,color='blue',linewidth=1.0)
        #circle1=plt.Circle((0,0),CTHCounterCentROuter,color='black',fill=False,ls='--')
        #circle2=plt.Circle((0,0),CTHCounterCentRInner,color='black',fill=False,ls='--')
        circle3=plt.Circle((0,0),CTHCounterMaxR,color='black',fill=False,lw=2.0)
        circle4=plt.Circle((0,0),CTHSupportInnerR,color='black',fill=False,lw=2.0) 
        #plt.gcf().gca().add_artist(circle1)
        #plt.gcf().gca().add_artist(circle2)
        plt.gcf().gca().add_artist(circle3)
        plt.gcf().gca().add_artist(circle4)
        plt.gca().add_collection(lc1[1])
        plt.gca().add_collection(lc2[1])
        #plt.gca().add_collection(lc3[1])
        if(drawHits==True):
            plt.gca().add_collection(lHitCounterDS[0])
            plt.gca().add_collection(lHitCounterDS[1])
            plt.gca().add_collection(pArrows[1])
        if(drawCdcHits==True):
            plt.scatter(cdcHits[0],cdcHits[1],s=10,color='red')
        plt.show()
