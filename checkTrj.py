import uproot

import numpy as np

import math

import matplotlib.pyplot as plt

def main(fname):
    ttree = uproot.open(fname)["trjTree"]
    print(ttree.keys())
    ptree = uproot.open(fname)["parentTree"]
    print(ptree.keys())
    
    trj_pId, trj_nhits, trj_pos, trj_mom = ttree.arrays(["pId", "numHits", "pos", "mom"], outputtype=tuple)
    prt_pId, prt_pos, prt_mom, prt_stop = ptree.arrays(["pId", "pos", "mom", "posF"], outputtype=tuple)
    print("#of trajectories: ",len(trj_pId))

    prt_pos_ent = prt_pos[prt_pos.x<3500]
    prt_mom_ent = prt_mom[prt_pos.x<3500]
    prt_pId_ent = prt_pId[prt_pos.x<3500]
    
    tposx=np.array(trj_pos.x)
    tposy=np.array(trj_pos.y)
    tposz=np.array(trj_pos.z-7650.0)
    tpost=np.array(trj_pos.t)
    
    pposx=np.array(prt_pos.x)
    pposy=np.array(prt_pos.y)
    pposz=np.array(prt_pos.z-7650.0)
    ppost=np.array(prt_pos.t)

    pposx2=np.array(prt_pos_ent.x)
    pposy2=np.array(prt_pos_ent.y)
    pposz2=np.array(prt_pos_ent.z-7650.0)
    ppost2=np.array(prt_pos_ent.t)
    
    plt.figure(figsize=(12,8))
    x_bins=[np.linspace( 3100, 9100,120),
            np.linspace( 3100, 9100,120)]
    y_bins=[np.linspace(-1500, 1500,100),
            np.linspace(-1500, 1500,100)]
    #plt.hist2d(tposx,tposz,bins=[x_bins[0],y_bins[0]])
    plt.scatter(tposx,tposz,s=8,label='Trajectories')
    plt.scatter(pposx,pposz,s=8,label='Parents'     )
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
    plt.xlim(3100,9100)
    plt.legend(loc='upper right')
    plt.show()
    
    plt.figure(figsize=(8,8))
    plt.scatter(tposz,tposy,s=8,label='Trajectories')
    plt.scatter(pposz,pposy,s=8,label='Parents'     )
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
    plt.hist(prt_mom.p,bins=pbins,label='Parents'     )
    plt.hist(prt_mom_ent.p,bins=pbins,label='Parents@Ent')
    plt.xlim(0,200)
    plt.legend(loc='upper right')
    plt.ylabel('Entries/MeV')
    plt.xlabel('Momentum [MeV/c]')
    plt.show()
#
#    plt.figure(figsize=(12,8))
#    plt.title('Pz')
#    plt.hist(trj_mom.x,bins=pbins,label='Trajectories')
#    plt.hist(prt_mom.x,bins=pbins,label='Parents'     )
#    plt.hist(prt_mom_ent.x,bins=pbins,label='Parents@Ent')
#    plt.xlim(0,200)
#    plt.legend(loc='upper right')
#    plt.show()

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
    plt.hist(trj_pId,bins=np.linspace(-24.5,24.5,50),label='Trajectories')
    plt.hist(prt_pId,bins=np.linspace(-24.5,24.5,50),label='Parents'     )
    plt.hist(prt_pId_ent,bins=np.linspace(-24.5,24.5,50),label='Parents@Ent')
    plt.title('PID of each trajectories associated to the trigger')
    plt.legend(loc='upper left')
    plt.xlabel('PDG code')
    plt.yscale('log')
    plt.show()

    #print(prt_pId)
