import uproot

import numpy as np

import math

import matplotlib.pyplot as plt

nBunch=1

pNames=['e-','e+','mu-','mu+','pi-','pi+','n','p','others']
pdgCodes=[11,-11,13,-13,-211,211,22,2112,2212,'others']
def makeBarVals(data):
    values = np.zeros(len(pdgCodes))
    for i in range(len(data)):
        if   data[i]==pdgCodes[0]: values[0]=values[0]+1.0
        elif data[i]==pdgCodes[1]: values[1]=values[1]+1.0
        elif data[i]==pdgCodes[2]: values[2]=values[2]+1.0
        elif data[i]==pdgCodes[3]: values[3]=values[3]+1.0
        elif data[i]==pdgCodes[4]: values[4]=values[4]+1.0
        elif data[i]==pdgCodes[5]: values[5]=values[5]+1.0
        elif data[i]==pdgCodes[6]: values[6]=values[6]+1.0
        elif data[i]==pdgCodes[7]: values[7]=values[7]+1.0
        elif data[i]==pdgCodes[8]: values[8]=values[8]+1.0
        else :                     values[9]=values[9]+1.0

    for j in range(len(values)):
        values[j] = values[j]/float(nBunch)
    
    return values

def setNumBunch(num):
    global nBunch
    nBunch=num

offsetX=7650.0

def main(fname):
    ttree = uproot.open(fname)["trjTree"]
    print(ttree.keys())
    ptree = uproot.open(fname)["parentTree"]
    print(ptree.keys())
    
    trj_pId, trj_nhits, trj_pos, trj_mom = ttree.arrays(["pId", "numHits", "pos", "mom"], outputtype=tuple)
    prt_pId, prt_pos, prt_mom, prt_stop, prt_lvl = ptree.arrays(["pId", "pos", "mom", "posF", "pLevel"], outputtype=tuple)
    print("#of trajectories: ",len(trj_pId))

    prt_pos_ent = prt_pos[prt_pos.x<3500]
    prt_mom_ent = prt_mom[prt_pos.x<3500]
    prt_pId_ent = prt_pId[prt_pos.x<3500]
    
    tposz=np.array(trj_pos.x)
    tposy=np.array(trj_pos.y)
    tposx=np.array(offsetX-trj_pos.z)
    tpost=np.array(trj_pos.t)
    tposr=np.hypot(tposx,tposy)
    
    pposz=np.array(prt_pos[prt_lvl==1].x)
    pposy=np.array(prt_pos[prt_lvl==1].y)
    pposx=np.array(offsetX-prt_pos[prt_lvl==1].z)
    ppost=np.array(prt_pos[prt_lvl==1].t)
    pposr=np.hypot(pposx,pposy)

    pposz_ent=np.array(prt_pos_ent.x)
    pposy_ent=np.array(prt_pos_ent.y)
    pposx_ent=np.array(offsetX-prt_pos_ent.z)
    ppost_ent=np.array(prt_pos_ent.t)

    gposz=np.array(prt_pos[prt_lvl==2].x)
    gposy=np.array(prt_pos[prt_lvl==2].y)
    gposx=np.array(offsetX-prt_pos[prt_lvl==2].z)
    gpost=np.array(prt_pos[prt_lvl==2].t)
    gposr=np.hypot(gposx,gposy)
    
    x_bins=[np.linspace( 3100, 9100,120),
            np.linspace( 3100, 9100,120)]
    y_bins=[np.linspace(-1500, 1500,100),
            np.linspace(-1500, 1500,100)]
    z_bins=np.linspace(3100,9100, 600)

    plt.figure(figsize=(12,8))
    plt.scatter(tposz,tposx,s=8,label='Trajectories')
    plt.scatter(pposz,pposx,s=8,label='Parents'     )
    plt.scatter(gposz,gposx,s=8,label='Grandparents')
    plt.title('Z-X 2D for origin of triggered trajectories')
    plt.xlim( 3100, 9100)
    plt.ylim(-1500, 1500)
    plt.legend(loc='upper right')
    plt.show()

    plt.figure(figsize=(12,8))
    plt.title('Z origin')
    plt.hist(tposz,bins=z_bins,label='Trajectories')
    plt.hist(pposz,bins=z_bins,label='Parents'     )
    plt.hist(gposz,bins=z_bins,label='Grandparents')
    plt.xlim(3100,9100)
    plt.legend(loc='upper right')
    plt.show()

    plt.figure(figsize=(8,8))
    plt.scatter(tposx,tposy,s=8,label='Trajectories')
    plt.scatter(pposx,pposy,s=8,label='Parents'     )
    plt.scatter(gposx,gposy,s=8,label='Grandparents')
    plt.title('X-Y 2D for origin of triggered trajectories')
    plt.xlim(-1500, 1500)
    plt.ylim(-1500, 1500)
    plt.legend(loc='upper right')
    plt.show()

    plt.figure(figsize=(8,6))
    r_bins=np.linspace(0,1500,300)
    plt.hist(tposr,bins=r_bins,label='Trajectories',density=True)
    plt.hist(pposr,bins=r_bins,label='Trajectories',density=True,alpha=0.7)
    plt.xlim(0,1500)
    plt.xlabel('R [mm]')
    plt.yscale('log')
    plt.show()

    tposx2=tposx[tposz>5000] 
    tposy2=tposy[tposz>5000] 
    tposz2=tposz[tposz>5000] 
    pposx2=pposx[pposz>5000] 
    pposy2=pposy[pposz>5000] 
    pposz2=pposz[pposz>5000] 
    gposx2=gposx[gposz>5000] 
    gposy2=gposy[gposz>5000] 
    gposz2=gposz[gposz>5000] 
    tposx3=tposx2[tposz2<8000] 
    tposy3=tposy2[tposz2<8000] 
    tposz3=tposz2[tposz2<8000] 
    tposr3=np.hypot(tposx3,tposy3)
    pposx3=pposx2[pposz2<8000] 
    pposy3=pposy2[pposz2<8000] 
    pposz3=pposz2[pposz2<8000] 
    pposr3=np.hypot(pposx3,pposy3)
    gposx3=gposx2[gposz2<8000] 
    gposy3=gposy2[gposz2<8000] 
    gposz3=gposz2[gposz2<8000] 

    ### zoomed plot around the detector..
    plt.figure(figsize=(12,8))
    plt.scatter(tposz3,tposx3,s=8,label='Trajectories')
    plt.scatter(pposz3,pposx3,s=8,label='Parents'     )
    plt.scatter(gposz3,gposx3,s=8,label='Grandparents')
    plt.title('Z-X 2D for origin of triggered trajectories')
    plt.xlim( 5000, 8000)
    plt.ylim( -600,  600)
    plt.legend(loc='upper right')
    plt.show()

    ### zoomed plot around the detector..
    plt.figure(figsize=(12,8))
    plt.title('Z origin')
    plt.hist(tposz3,bins=z_bins,label='Trajectories')
    plt.hist(pposz3,bins=z_bins,label='Parents'     )
    plt.hist(gposz3,bins=z_bins,label='Grandparents')
    plt.xlim(5000,8000)
    plt.legend(loc='upper right')
    plt.show()
    
    plt.figure(figsize=(8,8))
    plt.scatter(tposx3,tposy3,s=8,label='Trajectories')
    plt.scatter(pposx3,pposy3,s=8,label='Parents'     )
    plt.scatter(gposx3,gposy3,s=8,label='Grandparents')
    plt.title('X-Y 2D for origin of triggered trajectories')
    plt.xlim(-600, 600)
    plt.ylim(-600, 600)
    plt.legend(loc='upper right')
    plt.show()

    plt.figure(figsize=(8,6))
    r_bins=np.linspace(0,600,120)
    plt.hist(tposr3,bins=r_bins,label='Trajectories',density=True)
    plt.hist(pposr3,bins=r_bins,label='Trajectories',density=True,alpha=0.7)
    plt.xlim(0,600)
    plt.xlabel('R [mm]')
    plt.yscale('log')
    plt.show()
    
    plt.figure(figsize=(12,6))
    plt.title('Initial time distributions of trajectories and their parents')
    tbins=np.linspace(0,2400,120)
    plt.hist(tpost,bins=tbins,label='Trajectories')
    plt.hist(ppost,bins=tbins,label='Parents'     )
    plt.hist(ppost_ent,bins=tbins,label='Parents@Ent')
    plt.xlim(0,2400)
    plt.legend(loc='upper right')
    plt.xlabel('Time [ns]')
    plt.ylabel('Entries/(20 ns)')
    plt.show()

    pbins=np.linspace(0,200,200)
    plt.figure(figsize=(12,8))
    plt.title('Initial momentum distributions of trajectories and their parents')
    plt.hist(trj_mom.p,bins=pbins,label='Trajectories')
    plt.hist(prt_mom[prt_lvl==1].p,bins=pbins,label='Parents'     )
    plt.hist(prt_mom_ent.p,bins=pbins,label='Parents@Ent')
    plt.xlim(0,200)
    plt.legend(loc='upper right')
    plt.ylabel('Entries/MeV')
    plt.xlabel('Momentum [MeV/c]')
    plt.show()

    cbins=np.linspace(-1.1,1.1,100)
    plt.figure(figsize=(12,8))
    plt.title('cos')
    plt.hist(trj_mom.x/trj_mom.p,bins=cbins,label='Trajectories')
    plt.hist(prt_mom[prt_lvl==1].x/prt_mom[prt_lvl==1].p,bins=cbins,label='Parents'     )
    plt.hist(prt_mom_ent.x/prt_mom_ent.p,bins=cbins,label='Parents@Ent')
    plt.xlim(-1.1,1.1)
    plt.legend(loc='upper right')
    plt.show()

    # Estimate the 2D histogram
    nbins = 80
    H, xedges, yedges = np.histogram2d(pposx_ent,pposy_ent,bins=nbins)
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    # Mask zeros
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
    # Plot 2D histogram using pcolor
    fig2 = plt.figure(figsize=(9,8))
    plt.title('XY at the entrance')
    plt.pcolormesh(xedges,yedges,Hmasked)
    plt.xlim(-250,250)
    plt.xlabel('X [mm]')
    plt.ylabel('Y [mm]')
    plt.ylim(-250,250)
    plt.colorbar()
    plt.show()

    plt.figure(figsize=(12,4))
    left = np.arange(len(pdgCodes))
    width= 0.25
    label='Triggered trajectories'
    values1=makeBarVals(trj_pId)
    plt.bar(left,        values1,width=width,color='blue', label=label)
    label='Parents of trg. trjs.'
    values2=makeBarVals(prt_pId[prt_lvl==1])
    plt.bar(left+width,  values2,width=width,color='red',  label=label)
    label='Grandparents of trg. trjs.'
    values3=makeBarVals(prt_pId[prt_lvl==2])
    plt.bar(left+2*width,values3,width=width,color='green',label=label)
    plt.xticks(left+width/2, pdgCodes)
    plt.ylim(0.8/float(nBunch),len(trj_pId)/(float(nBunch)))
    plt.grid(True)
    plt.yscale('log')
    plt.title('PID of each trajectories associated to the trigger')
    plt.legend(loc='upper right')
    plt.ylabel('Number of entries')
    plt.xlabel('PDG code')
    plt.show()

    #print('Parents PID')
    #print(prt_pId[prt_lvl==1])
    #print('Grandparents PID')
    #print(prt_pId[prt_lvl==2])
