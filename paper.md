---
title: 'DiracSVT: Numerical solver of the Dirac equation with scalar,
vector and tensor potentials'
tags:
  - Dirac equation
  - Potentials
  - Numerical solver
authors:
  - name: Alexander Wallén Kiessling
    orcid: 0009-0008-5433-1258
    affiliation: "1" 
  - name: Daniel Karlsson
    affiliation: "1"
  - name: Yuxin Zhao
    affiliation: "2"
  - name: Chong Qi
    orcid: 0000-0002-1406-5695
    affiliation: "1"
    
affiliations:
 - name: Department of Physics, KTH Royal Institute of Technology, Sweden
   index: 1
 - name: Department of Physics, University of Florida, USA
   index: 2
date: 01 November 2023 
bibliography: paper.bib

---

# Summary
The Dirac equation is the fundamental relativistic wave equation in quantum mechanics that describes the behavior of spin-1/2 massive particles. It is the generalization of the Schrodinger equation to account for relativistic effects. The solving of Dirac equations with various potentials have played a very important role in many areas of fundamental physics including atomic, nuclear, hadron and molecular systems. The software DiracSVT focuses on numerically solving the Dirac equation with Scalar, Vector as well as Tensor potentials in the spherical coordinate space. The shooting method is utilized with a Runge-Kutta 4 integration scheme. A Universal Vector+Scalar potential parameterized in a Woods-Saxon form is provided, which reproduces well the known single-particle states around all doubly-magic nuclei and enables us to study the effect of Tensor potential on the shell evolution of exotic nuclei. The code is prepared in three different languages, namely Python, C++ and Matlab. The code can be easily extended to study other systems including atomic, hadron and molecular physics.

# Statement of need
DiracSVT solves the Dirac equation in the spherical coordinate space considering all Scalar, Vector as well as Tensor potentials. There is no such open source code available as far as we know. Previous works on solving the Dirac equation focus in
particular on the vector potential and have utilized various integration schemes. We consider in this code in particular the application to nuclear physics.  One of the most important aspects of contemporary nuclear physics is the study of the shell evolution of exotic nuclei, which refers to the potentially dramatic changes in the shell structure as one approaches the driplines with excessive protons or neutrons (Otsuka & others, 2020). The evolution of the shell structure, which is crucial for our understanding of nuclear stability as well as the origin of heavy elements, can be induced by the isospin dependence of the spin-orbital in various nonrelativistic mean field model approaches. The Dirac equation provides not only  a more fundamental understanding of the spin-orbit interaction and nuclear shell structure in terms of vector and scalar potentials but also provides a unique possibility to study the effect of tensor potential. The shell structure is often studied with various density functional approaches with effective nucleon-nucleon interactions. We believe the present code provides a much simpler way to study those effects and possibly at a much higher precison for the single-particle spectroscopy. It also enables us to study systematically the evolution of the pseudospin symmetry in atomic nuclei.
We hope it can be a useful tool in studying the structure of exotic nuclei as well as other
quantum systems where the spin-orbit and tensor effects can be important.

# Mathematics

We limit ourselves to spherical symmetry. The Dirac Hamiltonian \(H\) with a scalar potential \(S\), a vector potential \(V\) and a tensor potential \(U\) can be given by (see (Akcay, 2009) and (Hassanabadi & others, 2012))

$$\begin{equation}
   H=\vec{\alpha} \cdot \vec{p}+\beta(m+S)+V-i \beta \vec{\alpha} \cdot \hat{r} U.
\end{equation}$$

The corresponding Dirac equation for the radial wave function can be expressed as

$$\begin{eqnarray}\left(\frac{d}{d r}+\frac{\kappa}{r}-U(r)\right) g_{\kappa}(r)=(E+m-\Delta(r)) f_{\kappa}(r) \\ \left(\frac{d}{d r}-\frac{\kappa}{r}+U(r)\right) f_{\kappa}(r)=-(E-m-\Sigma(r)) g_{\kappa}(r)\end{eqnarray}$$

where $f$ and $g$ are the two components of the radial wave function and $\kappa=-(l+1)$ for $j=l+\frac{1}{2}$ and $\kappa=l$ for $j=l-\frac{1}{2}$.
As a matter of convenience, we have defined our potentials as $\Sigma = V + S$ and $\Delta = V - S$ as sum of the vector and scalar potentials. Note that the binding energy B is taken from the total energy $E = B + m.$

The potentials can be parameterized in a standard Woods-Saxon shape (see (Kennedy, 2002) and (Xu & Qi, 2013))

$$\begin{eqnarray}
 \Sigma = \frac{\Sigma_{0,n(p)}}{1+e^{\frac{r - R_{\sigma}}{a_{\sigma}}}} + V_{\mathrm{Coulomb}}\\
 \Delta = \frac{\Delta_{0,n(p)}}{1+e^{\frac{r - R_{\delta}}{a_{\delta}}}} + V_{\mathrm{Coulomb}},\end{eqnarray}$$
 \end{equation}$$

where $R$ and $a$ are the radius and diffuseness parameters.

The quantities $\Sigma_{0,n(p)},$ and $\Delta_{0,n(p)}$ are defined as follows for three different scenarios (for details see (Xu & Qi, 2013), where $n$ and $p$ stand for neutron and proton states respectively.
The coulomb barrier is defined as:

$$\begin{equation}
V_\mathrm{Coulomb} = 
\begin{cases}
	\alpha \frac{Z}{r} \qquad \qquad  r > r_{\sigma} \qquad \qquad \\
	\alpha \frac{Z (3r_{\sigma}^2 - r^2)}{2r_{\sigma}^3} \qquad r \leq r_{\sigma},
\end{cases}
\end{equation}$$

where $\alpha$ is the fine structure constant.

In addition we have the tensor potential defined similarly in a Woods-Saxon shape: 

$$\begin{equation}
  U(r) = \frac{U_0}{1+e^{\frac{\Delta^{\sigma}(r)}{a_{\sigma}}}}  
\end{equation}$$


The boundary condition for bound states is defined separately for cases with small and large radii. For small $r$, we have

$$\begin{eqnarray} 
f = 
\begin{cases}
	-a_0\left(\frac{-B+\Sigma}{k}\right) \epsilon^{l+2} \qquad k < 0 \\
	a_0\epsilon^{l}  \qquad \qquad \qquad \quad  k > 0 
\end{cases}\\
g = 
\begin{cases}
	a_0 \epsilon^{l+1} \qquad \qquad \qquad  k < 0   \\
	a_0\left(\frac{2m + B -\Delta}{k}\right) \epsilon^{l+1} \quad  k > 0.
\end{cases}
\end{eqnarray}$$

Both wave functions approach zero as the radius goes to infinity.

To handle unbound resonance states with positive energy, we implement for  large $r$ similar boundary conditions to Eq. (26) as in (Alonso & De Vincenzo & Mondino, 1997):

$$\psi = 
\begin{pmatrix}
	[C j_lkr) + D y_l(r)]Y_{j,l,j_z} \\
	\frac{\sqrt{B^2 + 2Bm}}{B+2m} [C j_{l'}(r) + D y_{l'}(r)]Y_{j,l',j_z}
\end{pmatrix}$$

Where $j$ and $y$ indicate the spherical Bessel functions, $l' = j\pm \frac{1}{2}$ and $j_z = -j, -j+1, \dots j.$ Here, $Y_{j,l,jz}$ are the spinor spherical harmonics as detailed in (Varshalovich, 1988). In the 1D case, one of the two components of the spinor spherical harmonics vanish and the solution is reduced to the two components f and g. As such the condition becomes

$$\begin{align} 
&f = C j_l(r) + D y_l(r),\\
&g = \frac{\sqrt{B^2 + 2Bm}}{B+2m} [C j_{l'}(r) + D y_{l'}(r)]. 
\end{align}$$

We attempt to conjoin an inner solution and an outer solution together, in
what is known as the shooting method (Silbar & Goldman, 2010). In our case, we employ shooting in two directions, both outward and
inward. We make use of the Runge-Kutta 4 integration scheme as an ordinary
differential equation solver. The outwards shooting solution starts from the origin. The inwards shooting originates from an asymptotic distance with the boundary conditions as described above.

With the code and parameterization of the potential as prodived, one is ready to do a systematic calculation on the shell evolution over the whole nuclear chart.

# References
H. Akcay, Dirac equation with scalar and vector quadratic potentials and
coulomb-like tensor potential, Physics Letters A 373 (2009) 616–620.

V. Alonso, S. De Vincenzo, L. Mondino, On the boundary conditions
for the dirac equation, European Journal of Physics 18 (5) (1997) 315.

H. Hassanabadi, E. Maghsoodi, S. Zarrinkamar, H. Rahimov, Dirac
equation under scalar, vector, and tensor cornell interactions, Advances
in High Energy Physics 2012 (2012) 1–17.

P. Kennedy, The woods–saxon potential in the dirac equation, Journal
of Physics A: Mathematical and General 35 (3) (2002) 689.

T. Otsuka, A. Gade, O. Sorlin, T. Suzuki, Y. Utsuno, Evolution of
shell structure in exotic nuclei, Rev. Mod. Phys. 92 (2020) 015002.
doi:10.1103/RevModPhys.92.015002.

R. R. Silbar, T. Goldman, Solving the radial dirac equations: a numerical odyssey, European journal of physics 32 (1) (2010) 217.

D. A. Varshalovich, A. N. Moskalev, V. K. Khersonskii, Quantum theory
of angular momentum, World Scientific, 1988.

Z. Xu, C. Qi, Shell evolution and its indication on the isospin
dependence of the spin–orbit splitting, Physics Letters B 724 (4) (2013)
247–252. doi:https://doi.org/10.1016/j.physletb.2013.06.018.


