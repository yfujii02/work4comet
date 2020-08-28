import sys

import CTHTriggerStudy

args = sys.argv
startRun=int(args[1])
numRun  =int(args[2])
iset    =int(args[3])
thrkev  =int(args[4])
dataDir =args[5]
scintOnly=int(args[6])
if scintOnly>0:scintOnly=1
print("Start Run: " + args[1])
print("End   Run: " + args[2])
strSetting=['CTH 48 Segments','CTH 64 Segments Scint+Cheren(10mm Acrylic)',
            'CTH 64 Segments Scint+Scint 5mmT','CTH 64 Segments Scint+Scint 7.5mmT',
            'CHT 64 Segments Scint+Scint 5mmT, tilt angle=17deg',
            'CTH 64 Segments Scint+Scint 5mmT / thinner lead shielding']
print(" Run description :")
print("    " + strSetting[iset])
print("Threshold: " + args[4] + " [keV]")
print(" Use Cherenkov if (0) --> " + args[6])
thr=0.001*float(thrkev)
CTHTriggerStudy.main_loop(iset,startRun,numRun,thr,dataDir,scintOnly)
