from curses import savetty
import numpy as np
import matplotlib.pyplot as plt


Gamma = 3#1
omega_eg = 3#3
omega_r = 0.5
n = 10
tmax = 2 # 0.5
dt = 0.001
steps = int(tmax/dt) + 1
M= 2000 # number of trajectories

Omega = 21*Gamma # Classical Rabi freq
Delta = 0


# Define the initial state of the system, the wave equation and the density matrix
# psi_0 = 1/np.sqrt(2)*np.array([1, 1])
psi_0 = np.array([0,1]) # excited state


# Define the operators used in the Lindblad equation
sigma_z = np.array([[1, 0], [0, -1]])
sigma_p = np.array([[0, 1], [0, 0]])
sigma_m = np.array([[0, 0], [1, 0]])
sigma_pm = np.matmul(sigma_p, sigma_m)
sigma_mp = np.matmul(sigma_m, sigma_p)
I = np.eye(len(psi_0)) # identity

H_r = 0.5*np.array([[0, Omega], [Omega, -2*Delta]]) # photon anti bunching case
H = 0.5*omega_eg*sigma_z + H_r # H = H_0 + H_r


# Define operators used for spontaneous emission without thermal cavity
L_m = np.sqrt(Gamma)*sigma_m#jump operator L-
L_mdag = np.conjugate(np.transpose(L_m))
# L_p = np.sqrt(Gamma*n)*sigma_p # L_+ jump operator
# L_pdag = np.conjugate(np.transpose(L_p))
J = H - 0.5*1j*(L_mdag@L_m)


# Initalize system
psi = np.zeros((M, steps, len(psi_0)), dtype="complex_") # containing M number of psi wave functions
psi[:,0] = psi_0 #inintal state
rho = np.zeros((steps, len(psi_0), len(psi_0)), dtype='complex_') # the rho is claculated after the algorithm is done and we obtain the same strucutre as in a)
t = np.linspace(0, tmax, steps)


# Store the probabilities of a quantum jump towards the ground state have happened
# Pms_list = []
# P0s_list = []
Pes_list = [] #prob of finding it in e state, contains M trajectrories

# implenting the algorithm
for m in range(M):
    # Pm_list = np.zeros(len(t))
    # P0_list = np.zeros(len(t))
    Pe_list = np.zeros(len(t))
    for i in range(steps-1):
        bra_psi = np.conjugate(np.transpose(psi[m,i]))
        P_m = bra_psi@L_mdag@L_m@psi[m,i]*dt
        # P_p = bra_psi@L_pdag@L_p@psi[m,i]*dt
        P_0 = 1-P_m
        r = np.random.uniform(0, 1)
        if r<P_0:
            psi[m,i+1] = (I-1j*dt*J)@psi[m,i]
        elif (r>P_0 and r<P_0+P_m): # pick smallest index
            psi[m,i+1] = L_m@psi[m,i]
        # else:
        #     psi[m,i+1] = L_p@psi[m,i]
        psi[m,i+1] /= np.linalg.norm(psi[m,i+1])
        
        Pe_list[i] = np.linalg.norm(np.conjugate(np.transpose(np.array([1,0]))) @ psi[m,i+1])**2

        # Pm_list[i] = P_m
        # P0_list[i] = P_0

    # Pms_list.append(Pm_list)
    # P0s_list.append(P0_list)
    Pes_list.append(Pe_list)



def mean_mat(mat):
    dim_axis_0 = mat.__len__()
    mean = [0 for i in range(dim_axis_0)]
    for vector in mat:
        for i, value in enumerate(vector):
            mean[i] += (value / dim_axis_0)
    return mean



# # Obtaining a list containing the average probability at each timestep
# average_Pe= []
# for i in range(len(t)):
#     for sublist in Pes_list:
#         average_Pe.append(sublist[i])

# for i in range(len(t)):
#     average_Pe[i] /= len(Pes_list)


average_Pe = [sum(row)/len(row) for row in zip(*Pes_list)]
        


spec_traj = 2

plt.figure()
plt.plot(Omega*t/(2*np.pi), Pes_list[spec_traj])
plt.xlabel("$\Omega_r t / 2 \pi$")
plt.ylabel(r"$P_e^{j}(t)$".format(j=spec_traj))
# plt.savefig(r"spontaneous_emission_j{j}.pdf".format(j=spec_traj))
plt.show()

plt.figure()
plt.plot(Omega*t[:-1]/(2*np.pi), average_Pe[:-1])
plt.xlabel("$\Omega_r t / 2 \pi$")
plt.ylabel(r"$p_e(t)$")
plt.savefig(r"driven_average_prob")
plt.show()




# plt.figure()
# for m in range(M):
#     plt.plot(t, psi[m,:,0])
# plt.xlabel("$t$")
# plt.ylabel("Re("r"$\rho_{ge}$"")")
# plt.show()


# Now for each timestep, calculate a new density matrix rho, and use the psi vector obatined from the algorithm to update the matrix

for i in range(steps):
    rho_i = np.zeros((len(psi_0), len(psi_0)), dtype='complex_')
    for m in range(M):
        rho_i += np.outer(psi[m,i], np.conjugate(psi[m,i]))
    rho[i] = rho_i/M


# Finally we plot the rho matrix just as before

# plt.figure()
# plt.plot(t, rho[:, 0, 0], label=r"$\rho_{ee}$") # rho_ee
# plt.plot(t, rho[:, 1, 1], label=r"$\rho_{gg}$") # rho_gg
# plt.plot(t, np.imag(rho[:, 0, 1]), label="Im("r"$\rho_{eg}$"")")#rho_eg
# plt.plot(t, np.imag(rho[:, 1, 0]), label="Im("r"$\rho_{ge}$"")")#rho_ge
# plt.plot(t, np.real(rho[:, 0, 1]), label="Re("r"$\rho_{eg}$"")")#rho_eg
# plt.plot(t, np.real(rho[:, 1, 0]), label="Re("r"$\rho_{ge}$"")")#rho_ge
# plt.xlabel(r"$t$")
# plt.ylabel(r"$\rho(t)$")
# plt.savefig(r"rho(t)_traj{j}.pdf".format(j=M))
# plt.legend(loc=4, fontsize=11)

# plt.show()
