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
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Daniel Karlsson
    affiliation: 1
  - name: Yuxin Zhao
    affiliation: 2
  - name: Chong Qi
    affiliation: 1
    
affiliations:
 - name: Department of Physics, KTH Royal Institute of Technology, Sweden
   index: 1
 - name: Department of Physics, University of Florida, USA
   index: 2
date: 01 November 2023 
bibliography: paper.bib

---

# Summary
The solving of the Dirac equation has played a very important role in many
areas of fundamental physics. In this work we present the Dirac equation
solver, DiracSVT, which solve the Dirac equation with Scalar, Vector as well
as Tensor nuclear potentials in the spherical coordinate space. The shooting
method is utilized with a Runge-Kutta 4 integration scheme. The potentials
are parameterized in a Woods-Saxon form, which reproduce well the known
single-particle states around all doubly-magic nuclei and can be applied to
study the shell evolution of exotic nuclei. The code can be easily extended
to study other systems including atomic, hadron and molecular physics.

The code solves the Dirac equation in the spherical
coordinate space utilizing the shooting method with a Runge-Kutta 4 inte-
gration scheme. The code is prepared in three different languages, namely 
Python, C++ and Matlab. Data is provided to readily study the available single-particle
states around doubly-magic nuclei that include 16O, 40Ca, 48Ca, 56Ni, 100Sn,
132Cs and 208Pb. The same potentials can also reproduce well the shell struc-
ture of many open nuclei over the whole nuclear chart.

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

We limit ourselves to spherical symmetry. The Dirac Hamiltonian \(H\) with a scalar potential \(S\), a vector potential \(V\) and a tensor potential \(U\) can be given by (see, e.g., SOURCE)

$$\begin{equation}
   H=\vec{\alpha} \cdot \vec{p}+\beta(m+S)+V-i \beta \vec{\alpha} \cdot \hat{r} U.
\end{equation}$$

The corresponding Dirac equation for the radial wave function can be expressed as

$$\begin{eqnarray}\left(\frac{d}{d r}+\frac{\kappa}{r}-U(r)\right) g_{\kappa}(r)=(E+m-\Delta(r)) f_{\kappa}(r) \\ \left(\frac{d}{d r}-\frac{\kappa}{r}+U(r)\right) f_{\kappa}(r)=-(E-m-\Sigma(r)) g_{\kappa}(r)\end{eqnarray}$$

where $f$ and $g$ are the two components of the radial wave function and $\kappa=-(l+1)$ for $j=l+\frac{1}{2}$ and $\kappa=l$ for $j=l-\frac{1}{2}$.
As a matter of convenience, we have defined our potentials as $\Sigma = V + S$ and $\Delta = V - S$ as sum of the vector and scalar potentials. The above equations can also be rewritten as

$$\begin{eqnarray}
g'_{\kappa} = (2m + B - \Delta)f_{\kappa} + \left(U-\frac{\kappa}{r}\right)g_{\kappa} \\
f'_{\kappa} = (-B + \Sigma)g_{\kappa} + \left(\frac{\kappa}{r} - U\right)f_{\kappa} 
\end{eqnarray}$$

where the binding energy B is taken from the total energy $E = B + m.$


The potentials can be parameterized in a standard Woods-Saxon shape (see, e.g., \cite{XU2013247,kennedy2002woods})

$$\begin{eqnarray}
 \Sigma = \frac{\Sigma_{0,n(p)}}{1+e^{\frac{r - R_{\sigma}}{a_{\sigma}}}} + V_{\mathrm{Coulomb}}\\
 \Delta = \frac{\Delta_{0,n(p)}}{1+e^{\frac{r - R_{\delta}}{a_{\delta}}}} + V_{\mathrm{Coulomb}},\end{eqnarray}$$
 
where $R$ and $a$ are the radius and diffuseness parameters.

The quantities $\Sigma_{0,n(p)},$ and $\Delta_{0,n(p)}$ are defined as follows for three different scenarios (for details see Ref. \cite{XU2013247}), where $n$ and $p$ stand for neutron and proton states respectively.

$$\begin{eqnarray}
\Sigma_{0,p} = V_0\left(1+\delta\frac{N-Z}{A}\right) \\
\Sigma_{0,n} = V_0\left(1-\delta\frac{N-Z}{A}\right)
\end{eqnarray}$$

Scenario 1:

$$\begin{eqnarray}
\Delta_0 = -\lambda \Sigma_0
\end{eqnarray}$$

Scenario 2:

$$\begin{eqnarray}
\Delta_{0,p} = -\lambda V_0\left(1-\delta\frac{N-Z}{A}\right) \\
\Delta_{0,n} = -\lambda V_0\left(1+\delta\frac{N-Z}{A}\right) 
\end{eqnarray}$$

Scenario 3:

$$\begin{eqnarray}
\Delta_{0,p} = -\lambda V_0\left(1-\delta_{so}\frac{N-Z}{A}\right) \\
\Delta_{0,n} = -\lambda V_0\left(1+\delta_{so}\frac{N-Z}{A}\right) 
\end{eqnarray}$$


The parameters $V_0, \delta, \delta_{so}, \lambda, a_{\sigma} = a_\delta, r_0,$ and $r_{0,ls}$ are fitted to data for the three scenarios separately. In our simplest case the diffuseness parameters  $a_{\sigma}$, $a_\delta$ are taken to be identical.

The coulomb barrier is defined as:

$$\begin{equation}
V_\mathrm{Coulomb} = 
\begin{cases}
	c \frac{Z}{r} \qquad \qquad  r > r_{\sigma} \qquad \qquad \\
	c \frac{Z (3r_{\sigma}^2 - r^2)}{2r_{\sigma}^3} \qquad r \leq r_{\sigma},
\end{cases}
\end{equation}$$

where $c = 0.0072923$.

In addition we have the tensor potential defined similarly in a Woods-Saxon shape: 

$$\begin{equation}
  U(r) = \frac{U_0}{1+e^{\frac{\Delta^{\sigma}(r)}{a_{\sigma}}}}  
\end{equation}$$


The boundary condition for bound states is defined separately for cases with small and large radii. For small $r$, we have

$$\begin{eqnarray} \label{boun1}
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

To handle unbound resonance states with positive energy, we implement for  large $r$ similar boundary conditions to Eq. (26) as in \cite{alonso1997boundary}:

$$\psi = 
\begin{pmatrix}
	[C j_lkr) + D y_l(r)]Y_{j,l,j_z} \\
	\frac{\sqrt{B^2 + 2Bm}}{B+2m} [C j_{l'}(r) + D y_{l'}(r)]Y_{j,l',j_z}
\end{pmatrix}$$

Where $j$ and $y$ indicate the spherical Bessel functions, $l' = j\pm \frac{1}{2}$ and $j_z = -j, -j+1, \dots j.$ Here, $Y_{j,l,jz}$ are the spinor spherical harmonics as detailed in \cite{varshalovich1988quantum}. In the 1D case, one of the two components of the spinor spherical harmonics vanish and the solution is reduced to the two components f and g. As such the condition becomes

$$\begin{align} 
&f = C j_l(r) + D y_l(r),\\
&g = \frac{\sqrt{B^2 + 2Bm}}{B+2m} [C j_{l'}(r) + D y_{l'}(r)]. \label{boun2}
\end{align}$$

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
