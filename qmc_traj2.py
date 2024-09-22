import numpy as np
import matplotlib.pyplot as plt

# Define the parameters
Gamma = 1#1
omega_eg = 3#3
omega_r = 0.5
n = 10

Omega = 21 # Classical Rabi freq
Delta = 5


# Define the system settings
tmax = 1 #0.5
dt = 0.001
steps = int(tmax/dt) + 1


# psi_0 = np.array([1/np.sqrt(2), 1/np.sqrt(2)])
psi_0 = np.array([0,1]) # excited state



# Define the time array that stores the time of the system
t = np.linspace(0, tmax, steps)


# An array containing the calcluated density matrix at each timestep, with an initial density matrix rhon_0 at the start
rhon = np.zeros((steps, len(psi_0), len(psi_0)), dtype="complex_")
rhon[0]= np.outer(psi_0, psi_0)



# Define the interacting Hamiltonian part for the driving Rabi freq part
# def H_int(Omegas, Deltas):
#      levels = len(Omegas) + 1
#      H = np.zeros((levels, levels))
#      for i in range(levels):
#           if np.logical_and()




# Define the operators used in the Lindblad equation
sigma_z = np.array([[1, 0], [0, -1]])
sigma_p = np.array([[0, 1], [0, 0]])
sigma_m = np.array([[0, 0], [1, 0]])
sigma_pm = np.matmul(sigma_p, sigma_m)
sigma_mp = np.matmul(sigma_m, sigma_p)
H_r = 0.5*np.array([[0, Omega], [Omega, -2*Delta]]) # photon anti bunching case
H = 0.5*omega_eg*sigma_z #+ H_r # H = H_0 + H_r


# Here we define som functions/operators needed to calculate the density matrix
# commutator = lambda A, B : np.matmul(A, B) - np.matmul(B, A)
def commutator(A, B):
    return np.matmul(A, B) - np.matmul(B, A)

def norm(v):
    return np.sqrt(np.conjugate(np.transpose(v))@v)


# Calculate the lindblad eq at each timestep, adding the contribution for each timestep as the time evolves (finite difference method)
for i in range(steps-1):
    rhon[i+1] = (rhon[i] + dt*(-1j*commutator(H, rhon[i]) + 
                Gamma*(n+1)*(np.matmul(sigma_m, np.matmul(rhon[i], sigma_p)) - 0.5*(np.matmul(sigma_p, np.matmul(sigma_m, rhon[i])) + np.matmul(rhon[i], np.matmul(sigma_p, sigma_m)))) + 
                Gamma*n*(np.matmul(sigma_p, np.matmul(rhon[i], sigma_m)) - 0.5*(np.matmul(sigma_m, np.matmul(sigma_p, rhon[i])) + np.matmul(rhon[i], np.matmul(sigma_m, sigma_p))))))



# def main():
    
#     plt.figure()
#     plt.plot(t, rhon[:, 0, 0], label=r"$\rho_{ee}$") # rho_ee 
#     plt.plot(t, rhon[:, 1, 1], label=r"$\rho_{gg}$") # rho_gg
#     # plt.plot(t, np.imag(rhon[:, 0, 1]), label="Im("r"$\rho_{eg}$"")")#rho_eg
#     # plt.plot(t, np.imag(rhon[:, 1, 0]), label="Im("r"$\rho_{ge}$"")")#rho_ge
#     # plt.plot(t, np.real(rhon[:, 0, 1]), label="Re("r"$\rho_{eg}$"")")#rho_eg
#     # plt.plot(t, np.real(rhon[:, 1, 0]), label="Re("r"$\rho_{ge}$"")")#rho_ge
#     plt.ylabel(r"$\rho(t)$")
#     plt.xlabel(r"$t$")
#     plt.legend()
#     plt.savefig("rho_Omega" + str(Omega) + "_Delta" + str(Delta) + ".pdf", format="pdf")
#     plt.show()


# if __name__ == "__main__":
# 	main()

plt.figure()
plt.plot(t, rhon[:, 0, 0], label="Master equation solution") # rho_ee 
# plt.plot(t, rhon[:, 1, 1], label=r"$\rho_{gg}$") # rho_gg
# plt.plot(t, np.imag(rhon[:, 0, 1]), label="Im("r"$\rho_{eg}$"")")#rho_eg
# plt.plot(t, np.imag(rhon[:, 1, 0]), label="Im("r"$\rho_{ge}$"")")#rho_ge
# plt.plot(t, np.real(rhon[:, 0, 1]), label="Re("r"$\rho_{eg}$"")")#rho_eg
# plt.plot(t, np.real(rhon[:, 1, 0]), label="Re("r"$\rho_{ge}$"")")#rho_ge
# plt.ylabel(r"$\rho(t)$")
# plt.xlabel(r"$t$")
# plt.legend()
# plt.savefig("rho_Omega" + str(Omega) + "_Delta" + str(Delta) + ".pdf", format="pdf")
# plt.show()










###########################################################################################################################################


def norm(vec):
    return np.sqrt(np.conjugate(np.transpose(vec))@vec)

#------params------
G=1#Gamma
w=3#omega
n=10#n
psi0=np.array([1/np.sqrt(2), 1/np.sqrt(2)])#psi(0)
# tmax=0.5#calculate time evolution until here
dt=0.0001#time step size
steps=int(tmax/dt)+1#+1 so that the point tmax is also included

rho = np.zeros((steps, len(psi0), len(psi0)), dtype = 'complex_')#an array with space for a density matrix at each time that we calculate
rho[0] = np.outer(psi0, psi0)#starting value for rho is the pure preparated state

sz = np.array([[1, 0], [0, -1]])#pauli matrix sigma_z
sm = np.array([[0, 0], [1, 0]])#simga_minus
sp = np.array([[0, 1], [0, 0]])#simga_plus
H = 1/2*w*sz#Hamiltonian




#problem 4b)
tmax=1 #0.5#calculate time evolution until here
dt=0.001#time step size
steps=int(tmax/dt)+1#+1 so that the point tmax is also included
M=2000#number of samples

# psi0=np.array([1/np.sqrt(2),1/np.sqrt(2)])#psi(0)
psi0 = np.array([0,1]) # excited state
I = np.eye(len(psi0))#identity
Lm = np.sqrt(G*(n+1))*sm#jump operator L-
Lmt = np.conjugate(np.transpose(Lm))
Lp = np.sqrt(G*n)*sp#jump operator L+
Lpt = np.conjugate(np.transpose(Lp))
J = H - 1j/2*(Lmt@Lm+Lpt@Lp)#non-hermitian operator for time evolution

psi = np.zeros((M, steps, len(psi0)), dtype='complex_')
psi[:,0] = psi0
rho = np.zeros((steps, len(psi0), len(psi0)), dtype='complex_')

for m in range(M):
    for i in range(steps-1):
        brapsi = np.conjugate(np.transpose(psi[m,i]))
        Pm = brapsi@Lmt@Lm@psi[m,i]*dt
        Pp = brapsi@Lpt@Lp@psi[m,i]*dt
        P0 = 1-Pm-Pp
        r = np.random.uniform(0, 1)
        if r<P0:
            psi[m,i+1] = (I-1j*dt*J)@psi[m,i]
        elif (r>P0 and r<P0+Pm):
            psi[m,i+1] = Lm@psi[m,i]
        else:
            psi[m,i+1] = Lp@psi[m,i]
        psi[m,i+1] /= norm(psi[m,i+1])

ts=np.linspace(0, tmax, steps)#all the values t takes
# for m in range(M):
#     plt.plot(ts, psi[m,:,0])







for i in range(steps):
    rhoi = np.zeros((len(psi0), len(psi0)), dtype='complex_')
    for m in range(M):
        rhoi += np.outer(psi[m,i], np.conjugate(psi[m,i]))
    rho[i] = rhoi/M
    #print(rho[i])
#plt.plot(ts, psi[4,:,0])
labels = [r"$\rho_{ee}$", r"$\rho_{gg}$", r"$Im(\rho_{eg})$", r"$Im(\rho_{ge})$", r"$Re(\rho_{eg})$", r"$Re(\rho_{ge})$"]
#one entry (for one specific time) of rho looks like this: [[rho_ee, rho_eg], [rho_ge, rho_gg]] -> e.g. rho_ee(t) = rho[t,0,0]
#to get the values of the rho entries at all calculated times, use the placeholder ':' instead of a specific value for t.
plt.figure()
plt.plot(t, rho[:, 0, 0], label=r"$\rho_{ee}$") # rho_ee
plt.plot(t, rho[:, 1, 1], label=r"$\rho_{gg}$") # rho_gg
plt.plot(t, np.imag(rho[:, 0, 1]), label="Im("r"$\rho_{eg}$"")")#rho_eg
plt.plot(t, np.imag(rho[:, 1, 0]), label="Im("r"$\rho_{ge}$"")")#rho_ge
plt.plot(t, np.real(rho[:, 0, 1]), label="Re("r"$\rho_{eg}$"")")#rho_eg
plt.plot(t, np.real(rho[:, 1, 0]), label="Re("r"$\rho_{ge}$"")")#rho_ge
plt.ylabel(r"$\rho(t)$")
# plt.ylabel(r"$\rho_{ee}$")
plt.xlabel(r"$t$")
plt.legend()
plt.show()