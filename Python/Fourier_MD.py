import numpy as np
from scipy.fftpack import fft,fftfreq,fftshift
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from math import factorial
print ' '

# Indicate which plots should show up

sig_plot=0                                      # Raw Signal
diff_sig_plot=1                                 # Positive signal minus negative Signal
diff_sig_plot_minus_decay=1                     # Positive signal minus negative Signal minus decay
Transform_plot=1                                # Fourier Transform
Gauss_plot=0                                    # Gaussian Fit of the Peak Frequency
temp_profs=1                                    # Temperature profiles


## Conditions

MAXSAWv=4000                                    # Maximum speed guess (m/s), used to limit range in which the peak frequency exists
removetransspike=0                              # Point (in GHz) where below, the transform is set to 0
widthGauss=5                                    # Distance (in GHz) in either direction used to search for peak
temp_index=25                                   # Time index when temperature fit starts


########## Define the smoothing function #############
def savitzky_golay(y, window_size, order, deriv=0, rate=1):

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')



###################### Functions to fit #######################

# Basic decay curve
def func(x, a, b, c):
    return a * np.exp(-b * (x/1e-9)) + c

# Gaussian fit for peak frequency
def funcgauss(x, a, b, c):
    return a * np.exp(-((x-b)/c)**2)

# Fit for acoustic damping parameter
def functimecon(x, a, b, c):
    return a * np.exp(-x/(b*1e-9)) + c

###################### Begin Data Reading For Acoustic Fit #######################

dt=0.2e-12
f=open('COM.txt','r')
f.readline()
wavelengthCOM=float(f.readline().split()[19])/10
freqcut=MAXSAWv/wavelengthCOM
f.seek(0)
g=open('DATA.txt','w')
g.write('#Parameters of MD run \n')
lengthCOM=len(f.readlines())-1
y=np.linspace(0.0,0.0,lengthCOM)
NewY=np.linspace(0.0,0.0,lengthCOM)
COM1z=np.linspace(0.0,0.0,lengthCOM) # z dimension for signal
COM5z=np.linspace(0.0,0.0,lengthCOM) # z dimension for signal

COM1x=np.linspace(0.0,0.0,lengthCOM) # x dimension for temperature
COM2x=np.linspace(0.0,0.0,lengthCOM) # x dimension for temperature
COM3x=np.linspace(0.0,0.0,lengthCOM) # x dimension for temperature
COM4x=np.linspace(0.0,0.0,lengthCOM) # x dimension for temperature
COM5x=np.linspace(0.0,0.0,lengthCOM) # x dimension for temperature
t=np.linspace(0.0,dt*lengthCOM,lengthCOM)

f.seek(0)
f.readline()
for i in xrange(lengthCOM):             # Read in the COM.txt data
    line=f.readline().split()
    drift=float(line[18])
    # Signal
    COM1z[i]=float(line[3])-drift
    COM5z[i]=float(line[15])-drift
    # Temperature
    COM1x[i]=float(line[1])
    COM2x[i]=float(line[4])
    COM3x[i]=float(line[7])
    COM4x[i]=float(line[10])
    COM5x[i]=float(line[13])
y=COM1z-COM5z
COM2x=(COM2x-COM1x)*10**(-10)
COM3x=(COM3x-COM1x)*10**(-10)
COM4x=(COM4x-COM1x)*10**(-10)
COM5x=(COM5x-COM1x)*10**(-10)
COM1x=(COM1x-COM1x)*10**(-10)


################################## End Data Reading ###################################

y=(y-y[0])*10**(-10)
COM1z=COM1z-COM1z[0]
COM5z=COM5z-COM5z[0]

try:        ######## Fit the exponential decay of POS-NEG ##########
    popt, pcov = curve_fit(func, t[0:lengthCOM], y[0:lengthCOM])
except Exception:
    print 'COULD NOT FIT EXPONENTIAL DECAY'

######  Make the signal data set of (POS-NEG) - decay curve with extra 0's so transform will be more refined
for i in xrange(lengthCOM):
    NewY[i]=y[i]-func(t[i], *popt)

############### FFT #####################
yf=fft(NewY)
tf = np.linspace(0.0, 1.0/(2.0*dt), lengthCOM//2)
tf=tf/1e9
dtf=tf[1]-tf[0]
yplot= (2.0/lengthCOM) * np.abs(yf[0:lengthCOM//2])

#######  Remove the damped "0" parts of data sets
NewY=NewY[0:lengthCOM]
y=y[0:lengthCOM]
t=t[0:lengthCOM]

####### Find first resonance peak frequency
maxcheck=freqcut/dtf
forindy=max(yplot)
avg_yplot=np.linspace(0.0,0.0,len(yplot))
avg_yplot[:]=forindy/2
numsteps=int(widthGauss/dtf)

#######  Find Lowest peak frequency to eliminate resonances
findmax=1
for i in xrange(len(yplot)):
    if yplot[i]==forindy:
        index=i
peakfreq_guess=tf[index]


######  Print the peak frequency guess
g.write( 'Peak Frequency Guess (GHz):\n ' + str(peakfreq_guess) + '\n')

#######  Use a gaussian fit to find peak frequency based on peak frequency guess
try:
    popt2, pcov2 = curve_fit(funcgauss, tf[index-numsteps:index+numsteps]-tf[index], yplot[index-numsteps:index+numsteps],maxfev=10000)
    g.write( 'Peak Frequency Fit (GHz):\n ' + str(popt2[1]+tf[index]) + '\n')
except Exception:
    g.write( 'Peak Frequency Fit (GHz):\n ' + 'COULD NOT FIT' + '\n')

freq=popt2[1]+tf[index]

t1=t[2000]
t2=t[2000]+1e-9/(freq)
T=int(round((t2-t1)/dt))
y1=max(NewY)
y2=min(NewY)

########  Print the wavelengthCOM and SAW speed
g.write( 'WavelengthCOM (nm):\n ' + str(wavelengthCOM) + '\n')
g.write( 'SAW speed (m/s):\n ' + str(wavelengthCOM*(popt2[1]+tf[index])) + '\n')

########  Initialize values to find the acoustic damping constant                                                                                                             # This needs to be fixed so that "iter" corresponds to the fit index and not the max index
numiter=(lengthCOM/(T))
localmaxpts=np.linspace(0.0,0.0,numiter)
shortx=np.linspace(0.0,0.0,numiter)

########  Iterate to find maximum in each frequency bin and its x-position
for i in xrange(int(numiter)):
    localmaxpts[i]=max(NewY[(i)*T:(i+1)*T])
    for j in xrange(T):
        if NewY[(i)*T+j]==localmaxpts[i]:
            shortx[i]=t[i*T+j]

#######  Fit and Print the acoustic damping constant values
try:
    popt3, pcov3 = curve_fit(functimecon, shortx, localmaxpts)
    g.write( 'Acoustic Damping Constant (ns):\n ' + str(popt3[1]) + '\n' )
except Exception:
    g.write( 'Acoustic Damping Constant (ns):\n ' + 'COULD NOT FIT' + '\n' )


#######  Fit and Print the acoustic damping constant values


# Fit for thermal diffusivity
if temp_profs:
    h=open('TEMP.txt','r')
    lengthTEMP=len(h.readlines())-1

    temptime=np.linspace(0.0,0.0,lengthTEMP)
    temp1=np.linspace(0.0,0.0,lengthTEMP)
    temp2=np.linspace(0.0,0.0,lengthTEMP)
    temp3=np.linspace(0.0,0.0,lengthTEMP)
    temp4=np.linspace(0.0,0.0,lengthTEMP)
    temp5=np.linspace(0.0,0.0,lengthTEMP)
    h.seek(0)
    h.readline()
    for i in xrange(lengthTEMP):            # Read in the TEMP.txt data
        line=h.readline().split()
        temptime[i]=float(line[0])
        temp1[i]=float(line[1])
        temp2[i]=float(line[2])
        temp3[i]=float(line[3])
        temp4[i]=float(line[4])
        temp5[i]=float(line[5])

    temp1=savitzky_golay(temp1,11,4)
    temp2=savitzky_golay(temp2,11,4)
    temp3=savitzky_golay(temp3,11,4)
    temp4=savitzky_golay(temp4,11,4)
    temp5=savitzky_golay(temp5,11,4)

    temptime=(temptime-temptime[0])*10**(-12)
    qo = (2*np.pi) / (wavelengthCOM*10**(-9))
    start_time = temp_index * (temptime[1]-temptime[0])

    def funcThermal1(x, a, k, c, b):
        return ( a / np.sqrt(x + start_time)) * (1 + b * np.exp(- (qo**2) * k * (x + start_time) * 1e-6) * np.cos(qo * np.mean(COM1x))) + c

    def funcThermal2(x, a, k, c, b):
        return ( a / np.sqrt(x + start_time)) * (1 + b * np.exp(- (qo**2) * k * (x + start_time) * 1e-6) * np.cos(qo * np.mean(COM2x))) + c

    def funcThermal3(x, a, k, c, b):
        return ( a / np.sqrt(x + start_time)) * (1 + b * np.exp(- (qo**2) * k * (x + start_time) * 1e-6) * np.cos(qo * np.mean(COM3x))) + c

    def funcThermal4(x, a, k, c, b):
        return ( a / np.sqrt(x + start_time)) * (1 + b * np.exp(- (qo**2) * k * (x + start_time) * 1e-6) * np.cos(qo * np.mean(COM4x))) + c

    def funcThermal5(x, a, k, c, b):
        return ( a / np.sqrt(x + start_time)) * (1 + b * np.exp(- (qo**2) * k * (x + start_time) * 1e-6) * np.cos(qo * np.mean(COM5x))) + c

    poptCOM1, pcovCOM1 = curve_fit(funcThermal1, temptime[temp_index:], temp1[temp_index:],maxfev=10000)
    poptCOM2, pcovCOM2 = curve_fit(funcThermal2, temptime[temp_index:], temp2[temp_index:],maxfev=10000)
    poptCOM3, pcovCOM3 = curve_fit(funcThermal3, temptime[temp_index:], temp3[temp_index:],maxfev=10000)
    poptCOM4, pcovCOM4 = curve_fit(funcThermal4, temptime[temp_index:], temp4[temp_index:],maxfev=10000)
    poptCOM5, pcovCOM5 = curve_fit(funcThermal5, temptime[temp_index:], temp5[temp_index:],maxfev=10000)


    ks=np.array([poptCOM1[1],poptCOM2[1],poptCOM3[1],poptCOM4[1],poptCOM5[1]])
    kavg=np.mean(ks)
    kstd=np.std(ks)

    g.write( 'Thermal Diffusivity (m^2/s):\n ' + str(kavg*1e-6) + '  (' + str(kstd*1e-6) + ')\n' )

    h.close()

################################## Plots ################################

if sig_plot:            # Raw Signal
    plt.figure(1)
    plt.xlabel('time (s)')
    plt.yticks([])
    plt.ylabel('Amplitude (arb.)')
    plt.title('Signal')
    plt.plot(t, COM1z)
    #plt.xlim(0,0.5)
    plt.grid()
    plt.show()

if diff_sig_plot:       # Positive signal minus negative Signal
    plt.figure(2)
    plt.xlabel('time (ns)')
    plt.yticks([])
    plt.ylabel('Amplitude (arb.)')
    plt.title('Positive Signal minus Negative Signal')
    plt.plot(t/1e-9, y)
    plt.plot(t/1e-9,func(t, *popt))
    plt.xlim(0,1.5)
    plt.grid()
    plt.savefig('Signal.png')
    # plt.show()

if diff_sig_plot_minus_decay:       # Positive signal minus negative Signal minus decay
    t1=t1/1e-9
    t2=t2/1e-9
    plt.figure(3)
    plt.xlabel('time (ns)')
    plt.yticks([])
    plt.ylabel('Amplitude (arb.)')
    plt.title('(Pos - Neg - Decay) and Acoustic Damping Parameter Fit')
    plt.plot(t/1e-9, NewY)
    plt.plot(shortx/1e-9,localmaxpts,'o')
    plt.plot(shortx/1e-9,functimecon(shortx,*popt3))
    plt.plot((t1,t1),(y1,y2),'k-')
    plt.plot((t2,t2),(y1,y2),'k-')
    plt.xlim(0,1.5)
    ax1 = plt.axes()  # standard axes
    ax2 = plt.axes([0.55, 0.2, 0.3, 0.3])
    plt.yticks([])
    plt.plot(t/1e-9, NewY)
    plt.plot(shortx/1e-9,localmaxpts,'o')
    plt.plot(shortx/1e-9,functimecon(shortx,*popt3))
    plt.plot((t1,t1),(y1,y2),'k-')
    plt.plot((t2,t2),(y1,y2),'k-')
    plt.xlim(0.38,0.435)
    plt.grid()
    plt.savefig('AcousDampFit.png')
    # plt.show()

if Transform_plot:      # Fourier Transform
    plt.figure(4)
    plt.xlabel('Frequency (GHz)')
    plt.yticks([])
    plt.ylabel('Amplitude (arb.)')
    plt.title('Fourier Transform of Signal')
    plt.plot(tf,yplot)
    plt.xlim(0,300)
    plt.grid()
    ax1 = plt.axes()  # standard axes
    ax2 = plt.axes([0.45, 0.45, 0.4, 0.4])
    plt.xlabel('time (ns)')
    plt.yticks([])
    # plt.ylabel('Amplitude (arb.)')
    plt.plot(t/1e-9, y)
    plt.plot(t/1e-9,func(t, *popt))
    plt.xlim(0,1.5)
    plt.grid()
    # plt.show()
    plt.savefig('Transform.png')

if Gauss_plot:              # Gaussian Fit of the Peak Frequency
    plt.figure(5)
    plt.xlabel('Frequency (GHz)')
    plt.yticks([])
    plt.ylabel('Amplitude (arb.)')
    plt.title('Peak Frequency and Gaussian Fit')
    plt.plot(tf[index-numsteps:index+numsteps],yplot[index-numsteps:index+numsteps])
    plt.plot(tf[index-numsteps:index+numsteps],funcgauss(tf[index-numsteps:index+numsteps]-tf[index],*popt2))
    plt.grid()
    plt.show()

if temp_profs:            # Raw Signal
    plt.figure(6)
    plt.xlabel('time (ns)')
    plt.ylabel('Temperature (K)')
    plt.title('Surface Temperature Profiles')
    plt.plot(temptime[0:lengthTEMP/2]/1e-9,temp1[0:lengthTEMP/2],'k-',label='Region 1')
    plt.plot(temptime[0:lengthTEMP/2]/1e-9,temp3[0:lengthTEMP/2],'b-',label='Region 3')
    plt.plot(temptime[0:lengthTEMP/2]/1e-9,temp5[0:lengthTEMP/2],'g-',label='Region 5')
    plt.plot(temptime[temp_index:lengthTEMP/2]/1e-9,funcThermal1(temptime[temp_index:lengthTEMP/2],*poptCOM1),'r--',linewidth=2)
    plt.plot(temptime[temp_index:lengthTEMP/2]/1e-9,funcThermal3(temptime[temp_index:lengthTEMP/2],*poptCOM3),'r--',linewidth=2)
    plt.plot(temptime[temp_index:lengthTEMP/2]/1e-9,funcThermal5(temptime[temp_index:lengthTEMP/2],*poptCOM5),'r--',linewidth=2)
    plt.legend()
    plt.xlim(0,0.5)
    plt.grid()
    # plt.show()
    plt.savefig('TempProf.png')


f.close()
g.close()
