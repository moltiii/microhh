from pylab import *

data100_2nd = loadtxt('conservation100_2nd/conservation.out', skiprows=1)
data200_2nd = loadtxt('conservation200_2nd/conservation.out', skiprows=1)
data400_2nd = loadtxt('conservation400_2nd/conservation.out', skiprows=1)

data100_4th = loadtxt('conservation100_4th/conservation.out', skiprows=1)
data200_4th = loadtxt('conservation200_4th/conservation.out', skiprows=1)
data400_4th = loadtxt('conservation400_4th/conservation.out', skiprows=1)

time100_2nd = data100_2nd[:,1]
mom100_2nd  = data100_2nd[:,7] / data100_2nd[1,7]
tke100_2nd  = data100_2nd[:,8] / data100_2nd[1,8]
mass100_2nd = data100_2nd[:,9] / data100_2nd[1,9]

time200_2nd = data200_2nd[:,1]
mom200_2nd  = data200_2nd[:,7] / data200_2nd[1,7]
tke200_2nd  = data200_2nd[:,8] / data200_2nd[1,8]
mass200_2nd = data200_2nd[:,9] / data200_2nd[1,9]

time400_2nd = data400_2nd[:,1]
mom400_2nd  = data400_2nd[:,7] / data400_2nd[1,7]
tke400_2nd  = data400_2nd[:,8] / data400_2nd[1,8]
mass400_2nd = data400_2nd[:,9] / data400_2nd[1,9]

time100_4th = data100_4th[:,1]
mom100_4th  = data100_4th[:,7] / data100_4th[1,7]
tke100_4th  = data100_4th[:,8] / data100_4th[1,8]
mass100_4th = data100_4th[:,9] / data100_4th[1,9]

time200_4th = data200_4th[:,1]
mom200_4th  = data200_4th[:,7] / data200_4th[1,7]
tke200_4th  = data200_4th[:,8] / data200_4th[1,8]
mass200_4th = data200_4th[:,9] / data200_4th[1,9]

time400_4th = data400_4th[:,1]
mom400_4th  = data400_4th[:,7] / data400_4th[1,7]
tke400_4th  = data400_4th[:,8] / data400_4th[1,8]
mass400_4th = data400_4th[:,9] / data400_4th[1,9]

"""
figure()
subplot(131)
plot(time100_2nd[1:200], mom100_2nd [1:200] - mom100_2nd [1], 'b-' , label = 'momdiff100 ')
plot(time200_2nd[1:400], mom200_2nd [1:400] - mom200_2nd [1], 'b--', label = 'momdiff200 ')
plot(time400_2nd[1:800], mom400_2nd [1:800] - mom400_2nd [1], 'b:' , label = 'momdiff400 ')
legend(loc=0, frameon=False)
subplot(132)
plot(time100_2nd[1:200], tke100_2nd [1:200] - tke100_2nd [1], 'b-' , label = 'tkediff100 ')
plot(time200_2nd[1:400], tke200_2nd [1:400] - tke200_2nd [1], 'b--', label = 'tkediff200 ')
plot(time400_2nd[1:800], tke400_2nd [1:800] - tke400_2nd [1], 'b:' , label = 'tkediff400 ')
legend(loc=0, frameon=False)
subplot(133)
plot(time100_2nd[1:200], mass100_2nd[1:200] - mass100_2nd[1], 'b-' , label = 'massdiff100')
plot(time200_2nd[1:400], mass200_2nd[1:400] - mass200_2nd[1], 'b--', label = 'massdiff200')
plot(time400_2nd[1:800], mass400_2nd[1:800] - mass400_2nd[1], 'b:' , label = 'massdiff400')
legend(loc=0, frameon=False)

figure()
subplot(131)
plot(time100_4th[1:200], mom100_4th [1:200] - mom100_4th [1], 'r-' , label = 'momdiff100 ')
plot(time200_4th[1:400], mom200_4th [1:400] - mom200_4th [1], 'r--', label = 'momdiff200 ')
plot(time400_4th[1:800], mom400_4th [1:800] - mom400_4th [1], 'r:' , label = 'momdiff400 ')
legend(loc=0, frameon=False)
subplot(132)
plot(time100_4th[1:200], tke100_4th [1:200] - tke100_4th [1], 'r-' , label = 'tkediff100 ')
plot(time200_4th[1:400], tke200_4th [1:400] - tke200_4th [1], 'r--', label = 'tkediff200 ')
plot(time400_4th[1:800], tke400_4th [1:800] - tke400_4th [1], 'r:' , label = 'tkediff400 ')
legend(loc=0, frameon=False)
subplot(133)
plot(time100_4th[1:200], mass100_4th[1:200] - mass100_4th[1], 'r-' , label = 'massdiff100')
plot(time200_4th[1:400], mass200_4th[1:400] - mass200_4th[1], 'r--', label = 'massdiff200')
plot(time400_4th[1:800], mass400_4th[1:800] - mass400_4th[1], 'r:' , label = 'massdiff400')
legend(loc=0, frameon=False)
"""

figure()
subplot(131)
plot(time100_2nd[1: 500], mom100_2nd [1: 500] - mom100_2nd [1], 'b-' , label = 'cfl08_2nd' )
plot(time100_4th[1: 500], mom100_4th [1: 500] - mom100_4th [1], 'r-' , label = 'cfl08_4th' )
plot(time200_2nd[1:1000], mom200_2nd [1:1000] - mom200_2nd [1], 'b--', label = 'cfl04_2nd' )
plot(time200_4th[1:1000], mom200_4th [1:1000] - mom200_4th [1], 'r--', label = 'cfl04_4th' )
plot(time400_2nd[1:2000], mom400_2nd [1:2000] - mom400_2nd [1], 'b:' , label = 'cfl02_2nd' )
plot(time400_4th[1:2000], mom400_4th [1:2000] - mom400_4th [1], 'r:' , label = 'cfl02_4th' )
legend(loc=0, frameon=False)
title('momentum conservation')
subplot(132)
plot(time100_2nd[1: 500], tke100_2nd [1: 500] - tke100_2nd [1], 'b-' , label = 'cfl08_2nd' )
plot(time100_4th[1: 500], tke100_4th [1: 500] - tke100_4th [1], 'r-' , label = 'cfl08_4th' )
plot(time200_2nd[1:1000], tke200_2nd [1:1000] - tke200_2nd [1], 'b--', label = 'cfl04_2nd' )
plot(time200_4th[1:1000], tke200_4th [1:1000] - tke200_4th [1], 'r--', label = 'cfl04_4th' )
plot(time400_2nd[1:2000], tke400_2nd [1:2000] - tke400_2nd [1], 'b:' , label = 'cfl02_2nd' )
plot(time400_4th[1:2000], tke400_4th [1:2000] - tke400_4th [1], 'r:' , label = 'cfl02_4th' )
legend(loc=0, frameon=False)
title('energy conservation')
subplot(133)
plot(time100_2nd[1: 500], mass100_2nd[1: 500] - mass100_2nd[1], 'b-' , label = 'cfl08_2nd')
plot(time100_4th[1: 500], mass100_4th[1: 500] - mass100_4th[1], 'r-' , label = 'cfl08_4th')
plot(time200_2nd[1:1000], mass200_2nd[1:1000] - mass200_2nd[1], 'b--', label = 'cfl04_2nd')
plot(time200_4th[1:1000], mass200_4th[1:1000] - mass200_4th[1], 'r--', label = 'cfl04_4th')
plot(time400_2nd[1:2000], mass400_2nd[1:2000] - mass400_2nd[1], 'b:' , label = 'cfl02_2nd')
plot(time400_4th[1:2000], mass400_4th[1:2000] - mass400_4th[1], 'r:' , label = 'cfl02_4th')
legend(loc=0, frameon=False)
title('mass conservation')

"""
timesteps = array([100, 200, 400])
momerror  = array([mom100[-1]-mom100[1], mom200[-1]-mom200[1], mom400[-1]-mom400[1]])
tkeerror  = array([tke100[-1]-tke100[1], tke200[-1]-tke200[1], tke400[-1]-tke400[1]])
masserror = array([mass100[-1]-mass100[1], mass200[-1]-mass200[1], mass400[-1]-mass400[1]])

figure()
subplot(131)
loglog(timesteps, abs(momerror))
subplot(132)
loglog(timesteps, abs(tkeerror))
subplot(133)
loglog(timesteps, abs(masserror))
"""
