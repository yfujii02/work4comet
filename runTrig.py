import sys

import CTHTriggerStudy

args = sys.argv
startRun=int(args[1])
endRun  =int(args[2])
iset    =int(args[3])
thrkev  =int(args[4])
print("Start Run: " + args[1])
print("End   Run: " + args[2])
strSetting=['CTH 48 Segments','CTH 64 Segments Scint+Cheren(10mm Acrylic)',
            'CTH 64 Segments Scint+Scint 5mmT','CTH 64 Segments Scint 5mmT + Scint 7.5mmT']
print(" Run description :")
print("    " + strSetting[iset])
print("Threshold: " + args[4] + " [keV]")
thr=0.001*float(thrkev)
dataDir="/Users/yfujii/work/COMET/data/"
CTHTriggerStudy.main_loop(iset,startRun,endRun,thr,dataDir)
