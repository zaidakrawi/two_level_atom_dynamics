# The Dynamics of a Two-Level Atom using Quantum Monte Carlo


This project investigates the dynamics of a simple two-level atom using a formalism of quantum mechanics that introduces the notion of random evolutions of the state vectors of a system. Over many realizations this method is shown to be successful in simulating the systems that will be investigated. The algorithm is first shown to simulate the simple case of a relaxing atom in a thermal cavity well. It is then applied to a system consisting of driven two-level atom undergoing spontaneous emission. Here, the algorithm also succeeds in providing insights into the anti-bunching phenomena occurring


## Introduction

In our study of physical systems, quantum mechanics has proven to be the most fundamental theory and
is at the very core of our present understanding of the laws of physics. To utilise this physical theory,
quantum mechanical systems must be studied as open systems, this is because in practice the systems
will almost always be subjected to a coupling to an uncontrolled environment. In many applications of
quantum physics, it is not possible to completely isolate the effects of the environment, neither is it feasible
to completely describe or control the degrees of freedom of the environment. Thus, one is forced to seek for
a simpler description of an open system’s dynamics – and one quickly realizes that this description is of a
probabilistic form.
In addition to the complexity of the systems, such a description is also warranted given the probabilistic
nature of quantum mechanics. There is thus a relation between the evolution of a system and the notion of
stochastic processes that will be of use later. A quantum mechanical master equation has been shown to
describe the dynamics of open quantum systems. Derived by G¨oran Lindblad, former professor at KTH, the
Lindblad equation is a form of the master equation utilised to study different open systems.
Even though the information of the system is not restricted – achieved by for instance treating the
environment as noise – the number of parameters required to be kept track of increases extremely fast as
the system increases in complexity. This is where computation power is required as analytical investigations
quickly become too tedious. A numerical scheme is needed to efficiently simulate these complex systems
with, hopefully, reasonable accuracy.
This project investigates the dynamics of a simple two-level atom using a formalism of quantum mechanics
that introduces the notion of random evolutions of the state vectors of a system. Over many realizations
this method is shown to be successful in simulating the systems that will be investigated. The project starts
by introducing how the dynamics of the system is obtained analytically, through the Lindblad equation, and
then introduces the stochastic algorithm which will be implemented to simulate the concerned systems. The
results of the implementation will then be presented and discussed.



## The Dynamics of the System

We would like to describe the evolution of a system coupled to an environment in thermal equilibrium. Focusing only on the relevant degrees of freedom is possible if we divide the whole system into two subsystems, referred to as simply the system $S$ and the environment $E$. The Hilbert space of the full system is then a tensor product of the two subsystems and we describe the state of the system by a density matrix $\rho$. This allows us to trace out the degrees of freedom of the environment from the full density matrix $\rho$ to obtain the reduced density matrix $\rho_S$ describing the state of the system.
