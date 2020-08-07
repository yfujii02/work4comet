import uproot

import numpy as np

import math

import matplotlib.pyplot as plt

#from scipy import stats

#### Constant parameters
ntmax=2400 # [ns]
nSeg=48 # altarnatively 64
nHod=2 # up/down stream
nType=2 # Cheren/scint
nChannels=nHod*nType*nSeg

trgSetStr=["Scint+Cheren w/  Birks law",
           "Scint-only   w/  Birks law",
           "Scint+Cheren w/o B",
           "Scint-only   w/o B"]

kBirks=0.125 # [mm/MeV] Assuming density = 1g/cm^3, https://en.wikipedia.org/wiki/Birks%27_law
### Paramters to mimic the waveform shaping (~9ns FWHM), G=2 (arbitrary)
shapeParams=np.array([1.0,1.6,1.85,2.0,1.95,1.75,1.5,1.25,1.0,0.75,0.45,0.3,0.15,0.075,0.03,0.01])
tExtend=len(shapeParams)
bunchCycle=1170
twindow1=[500,bunchCycle]
twindow2=[twindow1[0]+bunchCycle,twindow1[1]+bunchCycle]
specialDebug=False

#### Trigger settings
trgSetting=0 # 0: Scint+Cheren w/  Birks law
             # 1: Scint-only   w/  Birks law
             # 2: Scint+Cheren w/o B
eneThr=0.6 #[MeV]
#nphThr=10   # number of Cherenkov photons
nphThr=2   # number of Cherenkov photons
tSkip=400

refIndex=1.491
scintiThr=eneThr*2*50 #[MeV] x shaping x gain
cherenThr=nphThr*2 #Number of photon threshold x shaping
threshold=[0.0,0.0]

def trigger_setting(trgset,eth,pth,thickRatio):
    global trgSetting
    global eneThr
    global nphThr
    global scintiThr
    global cherenThr
    global threshold
    trgSetting=trgset
    eneThr=eth
    nphThr=pth
    scintiThr=eneThr*2*50 #[MeV] x shaping x gain
    cherenThr=nphThr*2 #Number of photon threshold x shaping
    if trgSetting%2==0:
        threshold=[cherenThr,scintiThr]
    else:
        threshold=[thickRatio*scintiThr,scintiThr]

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

def checkNumHits(wf):
    nhits=np.zeros(nChannels)
    for ibin in range(1,ntmax):
        for ich in range(nChannels):
            ihod,ityp,iseg=ch_to_info(ich)
            jtyp=(ityp+1)%2 # for different layer
            if wf[ich][ibin]>threshold[ityp] and wf[ich][ibin-1]<threshold[ityp]:
                nhits[ich]=nhits[ich]+1
                
    return nhits

time=np.arange(0,ntmax+tExtend,1.0)

def process_trigger(events,debugLvl,trgSetting):
    print("")
    print("Start of a file")
    
    nevents=events.numentries
    nacc=0
    nSingleHits=np.zeros(nChannels)
    for data in events.iterate(["eventId","cth_*"],entrysteps=1,namedecode="utf-8"):

        eventId=int(data["eventId"])
        if (eventId==4):continue ### Skip for new MC
        print("Event=",eventId)
        wf=np.zeros((nChannels,ntmax+tExtend),dtype=np.float32)
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
        ##### Birks' law
        cth_edep_corr=(data["cth_edep"]/(1.0+kBirks*data["cth_edep"]/data["cth_stepL"]))
        #nhits=int(data["cth_nHits"])
        nhits=len(cth_edep[0])
        
        for ihit in range(nhits):
            if cth_cryId[0][ihit]==11:continue ## Skip if this is a scintillator's light guide
            if cth_edep[0][ihit]<2e-6:continue ## Ignore too small energy deposition less than 2eV
            tbin=np.int32(cth_time[0][ihit])
            if tbin<0:continue
            ch=info_to_ch(cth_hodId[0][ihit],(cth_cryId[0][ihit]%10),cth_segId[0][ihit])
            output=0.0
            if (trgSetting%2==1 and cth_cryId[0][ihit]>1):continue
            
            if (trgSetting%2==1):### Scint-only
                if trgSetting<2:output=50*cth_edep_corr[0][ihit]
                else           :output=50*cth_edep[0][ihit]
            else:### Scint+Cheren
                if ((cth_cryId[0][ihit]%10)==0 and (refIndex*cth_beta[0][ihit])>1 ):
                    npho=cth_stepL[0][ihit]*47.5*0.1 # From Leo p.37 [n/mm],0.1=overall PDE (arbitrary)
                    npho=npho*(1.0-1.0/((refIndex*cth_beta[0][ihit])**2))
                    output=float(np.random.poisson(npho,1))
                elif cth_cryId[0][ihit]==1:
                    if trgSetting<2:output=50*cth_edep_corr[0][ihit]
                    else           :output=50*cth_edep[0][ihit]
            if tbin<ntmax and output>0.0:
                for ish in range(tExtend):wf[ch][tbin+ish]=wf[ch][tbin+ish]+shapeParams[ish]*output
        
        nhitstmp = checkNumHits(wf)
        nSingleHits = nSingleHits + nhitstmp
    print("End of a file")
    return nSingleHits

trgSets=[0]
debugLvl=0
# Turn interactive plotting off
#plt.ioff()

######### Main loop
def main_loop(case,startRun,endRun,dataDir):
    thickRatio=1.0
    thresholds=[0.05]
    global nSeg
    global nChannels
    global trgSets
    if case==0:
        dataDir=dataDir+"cth48/"
        thickRatio=1.8 ### Actually 2.0 but lowered to be conservative
        nSeg=48
        nChannels=nHod*nType*nSeg
    elif case==1:
        dataDir=dataDir+"cth64_SC/"
        thickRatio=1.8 ### Actually 2.0 but lowered to be conservative
        nSeg=64
        nChannels=nHod*nType*nSeg
    elif case==2:
        dataDir=dataDir+"cth64_SS/"
        thickRatio=1.0
        nSeg=64
        nChannels=nHod*nType*nSeg
        trgSets=[1]
    elif case==3:
        dataDir=dataDir+"cth64_SS_thicker/"
        ##thickRatio=1.4  ### Actually 1.5 but to be conservative
        thickRatio=1.0  ### New "thicker" geometry has same thickness in both layers
        nSeg=64
        nChannels=nHod*nType*nSeg
        trgSets=[1]
    
    fileList=[]
    t10n=np.arange(0,ntmax,10.0)
    
    #### Input files
    for i in range(startRun,endRun+1,1):
        num=str(i)
        #print(dataDir+"test.merged."+num.zfill(3)+".root")
        fileList.append(dataDir+"test.merged."+num.zfill(3)+".root")
    
    ####### Main scan loop
    for ith in range(len(thresholds)):
        for iset in range(len(trgSets)):
            trigger_setting(trgSets[iset],thresholds[ith],10,thickRatio)
            nSingleHits=np.zeros(nChannels) 
            print("********************************")
    
            print(trgSetStr[trgSetting])
            print("Energy threshold= ",eneThr," [MeV]")
            if trgSetting%2==0:
                print("#of Cherenkov photon threshold= ",nphThr)
            else:
                print("Enerty threshold for outer layer= ",thickRatio*eneThr," [MeV]")
    
            maxN=0
            for ifile in range(len(fileList)):
                events = uproot.open(fileList[ifile])["mc"]
                nhits  = process_trigger(events,debugLvl,trgSetting)
                nSingleHits = nSingleHits + nhits
                print(ifile+1,"/",len(fileList))

            print(len(nSingleHits))
            print(nSingleHits)
            hitRate=nSingleHits*1000/float(ntmax)/float(len(fileList)*4)
            print(hitRate)
            print("##################################")
            print("# Hit rate...   ")
            x=np.linspace(0,nChannels,nChannels)
            plt.plot(x,hitRate)
            plt.ylabel("Rate [MHz]")
            plt.xlabel("Channels")
            print("##################################")
            print("")
            plt.show()

