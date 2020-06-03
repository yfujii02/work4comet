import uproot

import numpy as np

import math

import matplotlib.pyplot as plt

pNames=['e-','e+','mu-','mu+','pi-','pi+','n','p','others']
pdgCodes=[11,-11,13,-13,-211,211,22,2112,2212,'others']
def makeBarVals(data):
    values = np.zeros(len(pdgCodes))
    for i in range(len(data)):
        if   data[i]==pdgCodes[0]: values[0]=values[0]+1
        elif data[i]==pdgCodes[1]: values[1]=values[1]+1
        elif data[i]==pdgCodes[2]: values[2]=values[2]+1
        elif data[i]==pdgCodes[3]: values[3]=values[3]+1
        elif data[i]==pdgCodes[4]: values[4]=values[4]+1
        elif data[i]==pdgCodes[5]: values[5]=values[5]+1
        elif data[i]==pdgCodes[6]: values[6]=values[6]+1
        elif data[i]==pdgCodes[7]: values[7]=values[7]+1
        elif data[i]==pdgCodes[8]: values[8]=values[8]+1
        else :                     values[9]=values[9]+1
    
    return values

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
    
    tposx=np.array(trj_pos.x)
    tposy=np.array(trj_pos.y)
    tposz=np.array(trj_pos.z-7650.0)
    tpost=np.array(trj_pos.t)
    
    pposx=np.array(prt_pos[prt_lvl==1].x)
    pposy=np.array(prt_pos[prt_lvl==1].y)
    pposz=np.array(prt_pos[prt_lvl==1].z-7650.0)
    ppost=np.array(prt_pos[prt_lvl==1].t)

    pposx2=np.array(prt_pos_ent.x)
    pposy2=np.array(prt_pos_ent.y)
    pposz2=np.array(prt_pos_ent.z-7650.0)
    ppost2=np.array(prt_pos_ent.t)

    gposx=np.array(prt_pos[prt_lvl==2].x)
    gposy=np.array(prt_pos[prt_lvl==2].y)
    gposz=np.array(prt_pos[prt_lvl==2].z-7650.0)
    gpost=np.array(prt_pos[prt_lvl==2].t)
    
    plt.figure(figsize=(12,8))
    x_bins=[np.linspace( 3100, 9100,120),
            np.linspace( 3100, 9100,120)]
    y_bins=[np.linspace(-1500, 1500,100),
            np.linspace(-1500, 1500,100)]
    #plt.hist2d(tposx,tposz,bins=[x_bins[0],y_bins[0]])
    plt.scatter(tposx,tposz,s=8,label='Trajectories')
    plt.scatter(pposx,pposz,s=8,label='Parents'     )
    plt.scatter(gposx,gposz,s=8,label='Grandparents')
    plt.title('Z-X 2D for origin of triggered trajectories')
    plt.xlim( 3100, 9100)
    plt.ylim(-1500, 1500)
    plt.legend(loc='upper right')
    plt.show()

    zbins=np.linspace(3100,9100, 600)
    plt.figure(figsize=(12,8))
    plt.title('Z origin')
    plt.hist(tposx,bins=zbins,label='Trajectories')
    plt.hist(pposx,bins=zbins,label='Parents'     )
    plt.hist(gposx,bins=zbins,label='Grandparents')
    plt.xlim(3100,9100)
    plt.legend(loc='upper right')
    plt.show()
    
    plt.figure(figsize=(8,8))
    plt.scatter(tposz,tposy,s=8,label='Trajectories')
    plt.scatter(pposz,pposy,s=8,label='Parents'     )
    plt.scatter(gposz,gposy,s=8,label='Grandparents')
    plt.title('X-Y 2D for origin of triggered trajectories')
    plt.xlim(-1500, 1500)
    plt.ylim(-1500, 1500)
    plt.legend(loc='upper right')
    plt.show()
    
    plt.figure(figsize=(12,6))
    plt.title('Initial time distributions of trajectories and their parents')
    tbins=np.linspace(0,2400,120)
    plt.hist(tpost,bins=tbins,label='Trajectories')
    plt.hist(ppost,bins=tbins,label='Parents'     )
    plt.hist(ppost2,bins=tbins,label='Parents@Ent')
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
    nbins = 50
    H, xedges, yedges = np.histogram2d(pposz2,pposy2,bins=nbins)
    # H needs to be rotated and flipped
    H = np.rot90(H)
    H = np.flipud(H)
    # Mask zeros
    Hmasked = np.ma.masked_where(H==0,H) # Mask pixels with a value of zero
    # Plot 2D histogram using pcolor
    fig2 = plt.figure(figsize=(9,8))
    plt.title('XY at the entrance')
    #plt.hist2d(pposz2,pposy2,bins=[np.linspace(-500,500,50),np.linspace(-500,500,50)],cmap=plt.cm.jet,vmin=0.5,vmax=4)
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
    plt.ylim(0.5,len(trj_pId))
    plt.grid(True)
    plt.yscale('log')
    plt.title('PID of each trajectories associated to the trigger')
    plt.legend(loc='upper right')
    plt.ylabel('Number of entries')
    plt.xlabel('PDG code')
    plt.show()

    print('Parents PID')
    print(prt_pId[prt_lvl==1])
    print('Grandparents PID')
    print(prt_pId[prt_lvl==2])
