from numpy import zeros, arange, array, str, exp, log, logspace
from matplotlib.pyplot import plot, show, xlabel, ylabel, title, hold, legend, yscale, xscale, ylim, xlim
from scipy import special


kerg = 1.380658e-16         # Boltzmann's constant [erg K]
h = 6.62607e-27             # Planck's constant [erg s]
c = 2.99792e10              # Speed of light [cm/s]
wav = arange(1000,20801,200)
b = zeros(wav.shape)


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
    plot(u[:], vau[i, :], label = 'a = ' + str(a[i]))
ylim(0,1)
xlim(-10,10)
legend(fontsize=8)
xlabel('u', size=14)
ylabel('Voigt function', size=12)
show()

for i in range(4):
    vau[i, :] = voigt(a[i], u[:])
    plot(u[:], vau[i, :], label = 'a = ' + str(a[i]))
yscale('log')
legend(fontsize=8, loc = 8)
xlabel('u', size=14)
ylabel('Logarithmic voigt function', size=12)
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
    plot(u, intensity)
xlabel('u', size=14)
ylabel('Intensity', size=12)
show()

logtau0 = arange(-2,2.1,0.5)

for itau in range(9):
    for i in range(201):
        tau = 10.**(logtau0[itau]) * voigt(a, u[i])
        intensity[i] = planck(Ts,wav) * exp(-tau) + planck(Tl,wav)*(1.-exp(-tau))
    plot(u,intensity, label = r'$\log{(\tau_0)} = $' + str(logtau0[itau]))
legend(loc=3, fontsize=8)
xlabel('u', size=14)
ylabel('Intensity', size=12)
show()

for iwav in range(1,4):
    wav = (iwav**2+1.)*1.0e-5 # wav = 2000, 5000, 10000 angstrom
    for itau in range(8):
        for i in range(201):
            tau = 10.**(logtau0[itau]) * voigt(a,u[i])
            intensity[i] = planck(Ts,wav) * exp(-tau) + planck(Tl,wav)*(1.-exp(-tau))
        intensity = intensity / intensity[0]
        plot(u,intensity[:], linewidth=1., label = r'$\log{(\tau_0)} = $' + str(logtau0[itau]))
legend(loc=3, fontsize=5)
xlabel('u', size=14)
ylabel('I/I(0)', size=12)
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
plot(u, intensity)           # absorption line
xlabel('u')
ylabel('Profile')
show()

# relative
reldepth = (intensity[0]-intensity)/intensity[0]
plot(u, reldepth)           # emission line
xlabel('u')
ylabel('Relative depth')
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

plot(tau0,eqw)
xlabel(r'$\tau_0$', size=12)
ylabel(r'Equivalent width $W_{\lambda}$', size=12)
xscale('log')
yscale('log')
show()
