from numpy import zeros, arange, array, str, exp, log, logspace
from matplotlib.pyplot import plot, show, xlabel, ylabel, title, hold, legend, yscale, xscale, ylim, xlim, grid
from scipy import special


kerg = 1.380658e-16         # Boltzmann's constant [erg K]
h = 6.62607e-27             # Planck's constant [erg s]
c = 2.99792e10              # Speed of light [cm/s]
wav = arange(1000,20801,200)
b = zeros(wav.shape)

B = 2.
tau = arange(0.01,10.01, 0.01)
intensity = zeros(tau.shape)

for I0 in range(4,-1,-1):
    intensity[:] = I0 * exp(-tau[:]) + B*(1-exp(-tau[:]))
    plot(tau, intensity, linewidth=1.0, label = 'intensity I(0) = ' + str(I0))
xlabel(r'Optical depth $\tau$', size=14)
ylabel('Intensity', size=14)
legend(fontsize=12)
grid()
show()


def planck(temp, wav):
    kT = kerg * temp
    B = (2.* h * c**2. / wav**5.) * (1./(exp(h*c / (wav*kT)) - 1))  # [erg/cm^2/s/cm/steradian]
    return B


def voigt(gamma,x):
    z = (x+1*gamma)
    V = special.wofz(z).real
    return V

u = arange(-10, 10.1, 0.1)
a = array([0.001, 0.01, 0.1, 1])
vau = zeros((a.shape[0], u.shape[0]))

for i in range(4):
    vau[i, :] = voigt(a[i], u[:])

for i in range(4):
    vau[i, :] = voigt(a[i], u[:])
    plot(u[:], vau[i, :], linewidth=1.0, label = 'a = ' + str(a[i]))
ylim(0,1)
xlim(-10,10)
legend(fontsize=12)
xlabel('u', size=14)
ylabel('Voigt function', size=14)
grid()
show()

for i in range(4):
    vau[i, :] = voigt(a[i], u[:])
    plot(u[:], vau[i, :], linewidth=1.0, label = 'a = ' + str(a[i]))
yscale('log')
legend(fontsize=12, loc = 8)
xlabel('u', size=14)
ylabel('Logarithmic voigt function', size=14)
grid()
show()

Ts = 5700.      # solar surface temperature
Tl = 4200.      # solar T-min temperature = 'reversing layer'
a = 0.1         # damping parameter
wav = 5000.0e-8 # wavelength in cm
tau0 = 1.       # reversing layer thickness at line center
u = arange(-10,10.1,0.1)
intensity = zeros(u.shape)

# Schuster-Schwarzchild line profile
for i in range(201):
    tau = tau0 * voigt(a, u[i])
    intensity[i] = planck(Ts,wav) * exp(-tau) + planck(Tl, wav)*(1.- exp(-tau))
    plot(u, intensity, linewidth=1.0)
xlabel('u', size=14)
ylabel('Intensity', size=14)
grid()
show()

logtau0 = arange(-2,2.1,0.5)

for itau in range(9):
    for i in range(201):
        tau = 10.**(logtau0[itau]) * voigt(a, u[i])
        intensity[i] = planck(Ts,wav) * exp(-tau) + planck(Tl,wav)*(1.-exp(-tau))
    plot(u,intensity, linewidth=1.0, label = r'$\log{(\tau_0)} = $' + str(logtau0[itau]))
legend(loc=3, fontsize=10)
xlabel('u', size=14)
ylabel('Intensity', size=14)
grid()
show()

for iwav in range(1,4):
    wav = (iwav**2+1.)*1.0e-5 # wav = 2000, 5000, 10000 angstrom
    for itau in range(8):
        for i in range(201):
            tau = 10.**(logtau0[itau]) * voigt(a,u[i])
            intensity[i] = planck(Ts,wav) * exp(-tau) + planck(Tl,wav)*(1.-exp(-tau))
        intensity = intensity / intensity[0]
        plot(u, intensity[:], linewidth=1.0, label = r'$\log{(\tau_0)} = $' + str(logtau0[itau]))
legend(loc=3, fontsize=6)
xlabel('u', size=14)
ylabel('I/I(0)', size=14)
grid()
show()


# 3.4 The equivalent width of spectral lines

def profile(a, tau0, u):
    Ts = 5700.
    Tl = 4200.
    wav = 5000.0e-8
    intensity = zeros(u.size)
    usize = u.size
    for i in range(usize):
        tau = tau0 * voigt(a, u[i])
        intensity[i] = planck(Ts,wav)*exp(-tau) + planck(Tl,wav)*(1.-exp(-tau))

    return intensity    # Schuster-Schwarzschild profile

# checking the profile
u = arange(-200,200.4,0.4)
a = 0.1
tau0 = 1.0e2
intensity = profile(a, tau0, u)
plot(u, intensity, linewidth=1.0)           # absorption line
xlabel('u', size=14)
ylabel('Profile', size=14)
grid()
show()

# relative
reldepth = (intensity[0]-intensity)/intensity[0]
plot(u, reldepth, linewidth=1.0)           # emission line
xlabel('u', size=14)
ylabel('Relative depth', size=14)
grid()
show()

eqw = sum(reldepth)*0.4
print(eqw)


# 3.5 The curve of growth

tau0 = logspace(-2, 4, 61)
eqw = zeros(tau0.size)
for i in range(61):
    intensity = profile(a,tau0[i],u)
    reldepth = (intensity[0] - intensity) / intensity[0]
    eqw[i] = sum(reldepth)*0.4

plot(tau0,eqw, linewidth=1.0)
xlabel(r'$\tau_0$', size=14)
ylabel(r'Equivalent width $W_{\lambda}$', size=14)
grid()
xscale('log')
yscale('log')
show()
