import uproot

import numpy as np

import math

import matplotlib.pyplot as plt

from scipy import stats

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
#shapeParams=np.array([1.0,1.25,1.5,1.75,2.0,1.75,1.5,1.25,1.0,0.75,0.5,0.25,0.1,0.05,0.02])
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
             # 3: Scint-only   w/o B

eneThr=0.6 #[MeV]
nphThr=10   # number of Cherenkov photons
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

def findClusterInSameLayer(wf,ihod,ityp,iseg,ibin,dseg,dbin,used):
    nCoinHit=1 # This function is called only when you already have one hit
    usedTmp=np.zeros(2*dseg+1)
    usedTmp[dseg]=1 # this channel itself
    for jseg in range(iseg-dseg,iseg+dseg):
        if jseg==iseg:continue
        jch=info_to_ch(ihod,ityp,jseg)
        for k in range(dbin):
            if ibin+k>=ntmax:break
            if used[jch][ibin+k]==1:continue
            if wf[jch][ibin+k]>threshold[ityp]:
                nCoinHit=nCoinHit+1
                usedTmp[jseg-iseg+dseg]=1
    return nCoinHit,usedTmp

def triggerSimulation(wf,nCoinThr1,nCoinThr2,dBin,dSeg1,dSeg2,skip):
    used=np.zeros((nChannels,ntmax))
    trigs=np.zeros(ntmax)
    nTrig=0
    nTrigAccepted=0
    for ibin in range(1,ntmax):
        if ((ibin%bunchCycle<tSkip-dBin) and skip):continue
        for ich in range(nChannels):
            if used[ich][ibin]==1: continue
            ihod,ityp,iseg=ch_to_info(ich)
            jtyp=(ityp+1)%2 # for different layer
            if wf[ich][ibin]>threshold[ityp] and wf[ich][ibin-1]<threshold[ityp]:
                nCoinHit1=1
                nCoinHit2=0
                usedTmp1=np.zeros(2*dSeg1+1)
                usedTmp2=np.zeros(2*dSeg1+1)
                nCoinHit1,usedTmp1=findClusterInSameLayer(wf,ihod,ityp,iseg,ibin,dSeg1,dBin,used)
                if nCoinHit1<nCoinThr1:continue
                centSeg2=iseg
                for jseg in range(iseg-dSeg2,iseg+dSeg2):# Check the different layer now
                    jch=info_to_ch(ihod,jtyp,jseg)
                    foundClust=False
                    for k in range(dBin):
                        if ibin+k>=ntmax:break
                        if used[jch][ibin+k]==1: continue
                        if wf[jch][ibin+k]>threshold[jtyp] and wf[jch][ibin+k-1]<threshold[jtyp]:
                            ######### Find cluster in the same layer...
                            nCoinHit2,usedTmp2=findClusterInSameLayer(wf,ihod,jtyp,jseg,ibin+k,dSeg1,dBin,used)
                            if nCoinHit2>=nCoinThr2:
                                centSeg2=jseg
                                foundClust=True
                                break
                        if foundClust==True:break
                
                if nCoinHit1>=nCoinThr1 and nCoinHit2>=nCoinThr2:
                    nTrig=nTrig+1
                    trigs[ibin]=trigs[ibin]+1
                    tacc=False
                    if ((ibin>=twindow1[0] and ibin<=twindow1[1]) or
                        (ibin>=twindow2[0] and ibin<=twindow2[1])):
                        tacc=True
                        print(" Trigger found within the time window!!!")
                        nTrigAccepted=nTrigAccepted+1
                    
                    if tacc==True :
                        print("  (ch,time,height)=(",ich,",",ibin,",",wf[ich][ibin],")")
                        print("   & contributing channels are:")
                    
                    # 10ns deadtime
                    for jseg in range(iseg-dSeg1,iseg+dSeg1):
                        jch=info_to_ch(ihod,ityp,jseg)
                        if usedTmp1[jseg-iseg+dSeg1]==1:
                            for j in range(10):
                                used[jch][ibin+j]=1
                                if(ibin+j==ntmax-1):break
                            if tacc==True : print("  ",jch)
                    for jseg in range(centSeg2-dSeg1,centSeg2+dSeg1):
                        jch=info_to_ch(ihod,jtyp,jseg)
                        if usedTmp2[jseg-centSeg2+dSeg1]==1:
                            for j in range(10):
                                used[jch][ibin+j]=1
                                if(ibin+j==ntmax-1):break
                            if tacc==True : print("  ",jch)
    return nTrig,nTrigAccepted,used,trigs

time=np.arange(0,ntmax+tExtend,1.0)

def process_trigger(events,debug,trgSetting,skipPrompt):
    print("")
    print("Start of a file")
    
    nevents=events.numentries
    nacc=0
    trigN_vs_time=np.zeros(ntmax)
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
            if skipPrompt==True and math.fmod(cth_time[0][ihit],float(bunchCycle))<tSkip: continue
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
        
        ###### Trigger
        triggers=np.zeros(ntmax)
        #nTrig,nTrigAccepted,trgch,triggers = triggerSimulation(wf,2,2,5,1,2,skipPrompt) #WF,ncoin1,ncoin2,dTime,dSeg(same layer),dSeg(different layer)
        nTrig,nTrigAccepted,trgch,triggers = triggerSimulation(wf,2,2,5,1,1,skipPrompt) #WF,ncoin1,ncoin2,dTime,dSeg(same layer),dSeg(different layer)
        nacc = nacc + nTrigAccepted
        print(" # of triggers=",nTrig," (accepted=",nTrigAccepted,")")
        
        for ibin in range(ntmax):
            trigN_vs_time[ibin]=trigN_vs_time[ibin]+triggers[ibin]
        ###### Draw waveforms
        if debug==True:###
            trigN=np.zeros(nChannels)
            trigT=[]
            trigA=[]
            for ich in range(nChannels):
                n=0
                tmpt=[]
                tmpa=[]
                for ibin in range(1,ntmax):
                    #### Get the earliest bin that contributes the trigger
                    if (trgch[ich][ibin]==1 and trgch[ich][ibin-1]==0):
                        tmpt.append(ibin)
                        tmpa.append(min(wf[ich][ibin],499))
                        if( (ibin>=twindow1[0] and ibin<=twindow1[1]) or
                            (ibin>=twindow2[0] and ibin<=twindow2[1])):
                            n=n+1
                trigT.append(np.array(tmpt))
                trigA.append(np.array(tmpa))
                trigN[ich]=n # Number of accepted trigger within the measurement time window
            trigT=np.array(trigT)
            trigA=np.array(trigA)
            #### Upstream
            #if (nTrigAccepted>0): print("  US Cherenkov <----> Scintillator")
            for ch in range(nSeg):
                if (trigN[ch]==0) and (trigN[ch+nSeg]==0):continue
                plt.figure(figsize=(15.0, 3.6))
                print("   CTH channel=",ch," & ",ch+nSeg)
                ### Cherenkov layer
                plt.subplot(1,2,1)
                plt.plot(time,wf[ch])
                plt.grid(True)
                plt.xlim(0, ntmax)
                plt.ylim(-1, 500)
                if trigN[ch]>0:
                    plt.scatter(trigT[ch],trigA[ch],color='red',linewidth=2.0)
                ### Scint layer
                plt.subplot(1,2,2)
                plt.plot(time,wf[ch+nSeg])
                plt.grid(True)
                plt.xlim(0, ntmax)
                plt.ylim(-1, 500)
                if trigN[ch+nSeg]>0:
                    plt.scatter(trigT[ch+nSeg],trigA[ch+nSeg],color='red',linewidth=2.0)
                plt.show()
            #### Downstream
            #if (nTrigAccepted>0): print("  DS Cherenkov <----> Scintillator")
            for ch in range(2*nSeg,3*nSeg):
                if (trigN[ch]==0) and (trigN[ch+nSeg]==0):continue
                plt.figure(figsize=(15.0, 3.6))
                print("   CTH channel=",ch," & ",ch+nSeg)
                ### Cherenkov layer
                plt.subplot(1,2,1)
                plt.plot(time,wf[ch])
                plt.grid(True)
                plt.xlim(0, ntmax)
                plt.ylim(-1, 500)
                if trigN[ch]>0:
                    plt.scatter(trigT[ch],trigA[ch],color='red',linewidth=2.0)
                ### Scint layer
                plt.subplot(1,2,2)
                plt.plot(time,wf[ch+nSeg])
                plt.grid(True)
                plt.xlim(0, ntmax)
                plt.ylim(-1, 500)
                if trigN[ch+nSeg]>0:
                    plt.scatter(trigT[ch+nSeg],trigA[ch+nSeg],color='red',linewidth=2.0)
                plt.show()
    print("End of a file")
    return nevents,trigN_vs_time,nacc

#thresholds=[0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85]
trgSets=[0,1]
debug=False
skipPrompt=True # Skip prompt (tSkip=400ns) to speed up the simulation
# Turn interactive plotting off
plt.ioff()

######### Main loop
def main_loop(case,startRun,endRun,thr,dataDir):
    thickRatio=1.0
    thresholds=[thr]
    global nSeg
    global nChannels
    global trgSets
    if case==0:
        dataDir=dataDir+"cth48/"
        thickRatio=1.75 ### Actually 2.0 but to be conservative
        nSeg=48
        nChannels=nHod*nType*nSeg
    elif case==1:
        dataDir=dataDir+"cth64_SC/"
        thickRatio=1.75 ### Actually 2.0 but to be conservative
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
        thickRatio=1.4  ### Actually 1.5 but to be conservative
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
    
            print("********************************")
            trigN_vs_time_All=np.zeros(int(ntmax/10))
    
            print(trgSetStr[trgSetting])
            print("Energy threshold= ",eneThr," [MeV]")
            if trgSetting%2==0:
                print("#of Cherenkov photon threshold= ",nphThr)
            else:
                print("Enerty threshold for outer layer= ",thickRatio*eneThr," [MeV]")
    
            naccTot=0
            maxN=0
            for ifile in range(len(fileList)):
                trigN_vs_time=np.zeros(ntmax)
                events = uproot.open(fileList[ifile])["mc"]
                nevents,trigN_vs_time,nacc = process_trigger(events,debug,trgSetting,skipPrompt)
                for i in range(int(ntmax/10)):
                    for j in range(10):
                        trigN_vs_time_All[i]=trigN_vs_time_All[i]+trigN_vs_time[i*10+j]
                    if maxN<trigN_vs_time_All[i]:
                        maxN=trigN_vs_time_All[i]
                print(ifile+1,"/",len(fileList))
                naccTot=naccTot+nacc
            
            print("##################################")
            print("# Number of accepted trigger is   ")
            print("#      ",naccTot)
            #print("#   out of ",len(fileList)*2*5," bunches")
            print("#   out of ",len(fileList)*2*4," bunches") # for new files
            print("##################################")
    
            print("#of Triggers vs timing")
            fig=plt.figure(figsize=(15.0, 3.6))
            plt.plot(t10n,trigN_vs_time_All,color='red',linewidth=2.0)
            plt.grid(True)
            plt.axvspan(twindow1[0], twindow1[1], facecolor='blue', alpha=0.4)
            plt.axvspan(twindow2[0], twindow2[1], facecolor='blue', alpha=0.4)
            plt.title("#of triggers vs timing")
            plt.xlabel('time [ns]')
            plt.ylabel('#of triggers / 10ns')
            plt.xlim(-1, ntmax+1)
            plt.ylim(-1, maxN+1)
            #plt.show()
            plt.savefig('figs2/trgDistDataset'+str(case)+'Run'+str(startRun)+'-'+str(endRun)
                        +'th'+str(int(1000*thresholds[ith]))+'keV_TrgSet'+str(trgSets[iset])+'.png')
            plt.close(fig)
            print("********************************")
            print("")

