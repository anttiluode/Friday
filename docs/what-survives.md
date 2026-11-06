# What survives, what can be read

*A research note on symmetry, hidden modes, and controlled access — 6 November 2026.*

The discussion began with π: why a simple rotation generator closes after 2π, how integer steps nearly close into cycles, and whether rational versus irrational phase relationships might organize computation. Fourier flowers, magnetic graphs, and writable wave media seemed to be different views of one mechanism.

There are real mathematical connections among these subjects. A shared vocabulary is not enough to establish a shared computational advantage. The productive step was to ask what the supplied graph actually does, identify an exact obstruction, and distinguish changes that affect it from changes that cannot.

The resulting lesson is specific: **a system can contain differences that its communication interface neither detects nor excites. Their accessibility is a property of the dynamics, the interface, and the allowed interventions together.**

## 1. A spectral observation with an exact explanation

The supplied `magnetic.py` builds a Hermitian magnetic Laplacian

$$L=D-H,$$

on five six-node rings. On an oriented internal edge from local node j to j+1, the adjacency entry is

$$H_{j,j+1}=\exp(i\pi j/3),\qquad j=0,\ldots,5.$$

The reverse entry is its complex conjugate. Node 0 of each ring connects to node 0 of the next ring. All edge magnitudes are one. There are 30 vertices and 35 undirected edges.

The original comparison changes the five outer phases from `2π(c+1)/5` to `(c+1)` times the golden angle. The reported repeated-eigenvalue count falls from 21 to 12, while the smallest eigenvalue remains 0.267949.

The internal phases sum to 5π, equivalent to π around a closed loop. They are not globally compatible: no assignment of vertex phases satisfies every internal connection simultaneously. “Rational edges” does not imply “an unfrustrated cycle.”

The unchanged minimum is exactly 2−√3. It is part of a three-eigenvalue family that survives every alteration of the outer phases.

## 2. An explicit basis for the hidden modes

For one ring, write

$$
\alpha_r=\frac{(2r+1)\pi}{6},\qquad r=0,1,2.
$$

Define the normalized local vectors

$$
v_r(j)=\frac{1}{\sqrt3}
\exp\!\left(-\frac{i\pi j(j-1)}6\right)
\sin(\alpha_r j),\qquad j=0,\ldots,5.
$$

These vectors vanish at j=0. Substitution into the ring Laplacian gives

$$
L_{\rm ring}v_r=(2-2\cos\alpha_r)v_r.
$$

For the interior nodes, the exponential factor removes the local edge phases, reducing the calculation to the usual sine recurrence. The connecting-node equation also holds: the contributions from nodes 1 and 5 cancel because `sin(5α_r) = sin(α_r)`, with the required opposite phase. The three vectors are orthonormal.

Extend each vector by zero outside its ring. Since its connecting-node amplitude is zero, every outer edge annihilates it. Five rings give an orthonormal 30×15 matrix Q satisfying

$$
LQ=Q\Lambda,\qquad CQ=0,
$$

where C reads the five connecting nodes and Λ contains five copies of each of

$$2-\sqrt3,\quad 2,\quad 2+\sqrt3.$$

This is an analytic construction. The numerical audit verifies its residuals at approximately machine precision.

Two consequences follow immediately.

**The 12 remaining repeated gaps are accounted for.** In the golden configuration, each of these eigenvalues has multiplicity five, contributing 3×4 repeated adjacent gaps. The numerical spectrum has no additional degeneracy at the stated tolerance. The original rational configuration has additional coincidences.

**The spectral minimum cannot move under outer phase changes.** Each isolated ring has minimum 2−√3. An outer magnetic edge with nonnegative weight contributes a positive semidefinite term, because its quadratic energy is a weighted squared phase-adjusted difference. Thus adding such edges cannot lower the minimum. The extended local modes attain the old minimum, so equality remains.

The same argument preserves this 15-dimensional subspace for arbitrary nonnegative outer weights and arbitrary connections supported only on the five connecting nodes, with the internal rings fixed. The runnable audit directly checks the original outer topology, including several arbitrary phase choices; the broader statement follows from the zero-support argument.

## 3. What the golden angle changes

The connected graph has cycle-space dimension

$$E-V+1=35-30+1=6.$$

Five independent cycle fluxes belong to the internal rings. The remaining one is the phase accumulated around the outer ring. With the internal edges fixed, the spectrum depends on the outer phases only through that total flux modulo 2π.

A constant phase change on all six vertices of one flower changes the phases of its outer edges without changing its internal edges. Successive choices move all the outer flux onto one edge. The resulting matrix is related to the original by a diagonal unitary change of basis and therefore has the same eigenvalues. This is standard magnetic gauge equivalence [1].

The audit compares the golden phase pattern with `[0, 0, 0, 0, sum(golden_phases)]`; their spectra agree to rounding precision. A separate rational control, using `2π(c+1)/11`, also removes the additional degeneracies and gives the same observable dimensions as the golden configuration for the tested sensors.

The demonstrated mechanism is a change in loop flux and degeneracy structure. It does not identify an advantage specific to irrationality. Finite matrices also depend continuously on their edge phases: sufficiently close rational fractions of a turn approximate these spectral results. That continuity does not negate asymptotic distinctions between rational and irrational rotations; it limits what this finite experiment establishes.

Gauge-equivalent spectra need not mean identical raw responses with every laboratory phase convention held fixed. Observables and input conventions must be transformed consistently when claiming equivalence of a full experiment. Coordinate sensors retain the observable dimensions used here under these diagonal phase changes.

## 4. Observability is a property of an interface

To discuss observable state, this note adds a dynamical model that the original script did not specify:

$$\dot x=-iLx+Bu,\qquad y=Cx.$$

States and outputs are complex amplitudes. For node actuation and measurement at the same connecting nodes, B=C*.

If a hidden mode starts with amplitude a, its contribution to the output is

$$Ce^{-itL}v_r a=e^{-it\lambda_r}Cv_r a=0.$$

Conversely, `Q*B = 0`: the same ports cannot write those amplitudes through additive input. This reciprocal lack of access is not a shortage of measurement duration.

For a Hermitian L, the continuous-time observable dimension can be computed by restricting C to each eigenspace and summing those ranks. It is essential to treat each degenerate eigenspace as a whole: individual eigenvectors returned by a numerical solver are not uniquely defined within it. This is the PBH observability criterion specialized to this setting [2].

The measured dimensions are:

| Sensors, zero-based indices | Original rational | Golden | Rational 1/11 turns |
|---|---:|---:|---:|
| 0, 6, 12, 18, 24 | 15 | 15 | 15 |
| 1, 7, 13, 19, 25 | 27 | 30 | 30 |
| Both sets | 30 | 30 | 30 |

These are complex state dimensions out of 30. They assume ideal complex readings, not intensity measurements. Numerical rank tolerances approximate an algebraic question; they are not a specified noise floor.

## 5. Why more frequencies cannot always help

At any z outside the spectrum,

$$C(zI-L)^{-1}Q=0.$$

Every derivative of this expression with respect to z also vanishes. Thus a frequency sweep or a packet of frequency derivatives cannot expose these initial-state directions through C. For collocated port inputs, the same modes are also absent from the port-to-port transfer response.

Even time-dependent outer phases do not solve the problem. Every allowed outer operator has the same common invariant hidden subspace:

$$L(t)Q=Q\Lambda,\qquad CQ=0.$$

The hidden component therefore evolves within that subspace and remains unseen. Adaptive choices based on past readings cannot distinguish two states differing only in it: they receive identical readings and hence choose the same interventions.

This is the precise limit of “ask a different question.” It works only if the permitted family of questions reaches different directions.

There are three distinct difficulties:

| Difficulty | Appropriate change |
|---|---|
| A direction is visible but weak | Improve signal, integration, or conditioning |
| A direction is hidden from this probe but visible to another | Change sensor, input, or measurement |
| A direction is hidden under the entire allowed family | Change the coupling, dynamics, or allowed operations |

Finite-difference derivatives computed from already collected samples do not create new information by themselves; they transform those samples and their noise. Direct derivative-sensitive measurements could be useful, but must be compared at a matched measurement budget.

## 6. Opening access, and the difference between rank and strength

Keep the golden outer edges and the five connecting-node sensors. Increase the phase of edge (0,1) within every flower by δ, updating the conjugate entry as well. This changes the internal loop flux and breaks the cancellation that kept the connecting node still.

Define the finite observation stack

$$
\mathcal O_{T}=\frac1{\sqrt T}
\begin{bmatrix}
C\\
Ce^{-i\Delta t L}\\
\vdots\\
Ce^{-i(T-1)\Delta t L}
\end{bmatrix},\qquad T=80,\quad\Delta t=0.37.
$$

The audit computes this stack in an orthonormal eigenbasis, which preserves singular values. All five sensors are read at every time, including zero. No noise is added in this diagnostic.

| δ | Continuous-time observable dimension | Smallest singular value of normalized stack |
|---|---:|---:|
| 0 | 15 | Roundoff zero |
| 0.0001 | 30 | 0.0000072518 |
| 0.01 | 30 | 0.00072518 |
| 0.1 | 30 | 0.0072454 |
| 0.5 | 30 | 0.035361 |

The near-linear scaling in this small-δ range is a measured local trend, not a global theorem about every perturbation. It demonstrates that algebraic access can be restored while practical access remains weak. For a full-column-rank observation matrix, the norm of the least-squares inverse is `1/σ_min`; noise in an unfavorable output direction can therefore be strongly amplified.

This is not a measured decoding error rate. Such a claim would require a payload ensemble, a noise model, a decoder, and a read-time budget. It is also not evidence that arbitrarily fast sampling guarantees observability: sampling can introduce additional aliases. Here the finite stack is checked explicitly.

## 7. Isolation, retention, and writing are different requirements

The phrase “protected memory” needs care.

Under ideal unitary evolution, a dark eigenmode retains its magnitude and accumulates a known phase. Under diffusion,

$$\dot x=-Lx,\qquad x(t)=e^{-\lambda t}v,$$

the very same mode decays. All three hidden eigenvalues are positive. No permanent state is supplied by this graph's dark subspace under that dynamics.

Even in the unitary model, protection here means isolation from the specified outer connections. It does not mean immunity to internal disorder, damping, sensor changes, or nonlinear couplings. The internal phase perturbation used to open access is itself an example of such a symmetry-breaking change.

There is also a writing obstruction: an inaccessible mode cannot be loaded through the same inaccessible ports. A complete memory architecture would need a write operation that couples into it, a retention phase with controlled disturbance, and a read operation that exposes the required payload. Turning the internal phase back to zero does not automatically restore a state disturbed during opening.

What is established here is the distinction between a closed interface and an opened one, plus the strength of the resulting observation map. A write–retain–read protocol remains untested. The next useful evidence would concern that complete cycle at a specified noise and energy budget, rather than another favorable eigenvalue statistic.

## 8. What the π paper contributes

The external manuscript *The irrationality exponent of π is 2* [3] motivates the discussion but is not an assumption needed for any result in this repository.

Its central architecture, as presented in the supplied September 24, 2026 manuscript, separates three steps:

1. A separated-weight interpolation theorem guarantees independent rows and a nonzero determinant.
2. Gaussian-rational entries give an arithmetic lower bound after denominators are cleared.
3. Near-period expansions, repeated Taylor-degree cancellations, and approximation-error factors give an incompatible analytic upper bound under the contradiction hypothesis.

The subtle point is not merely that near-resonances make a matrix nearly singular. The argument first establishes nonvanishing, then accounts for how small a nonzero determinant can be given its denominators. An exactly zero determinant would not produce that contradiction.

Our graph illustrates a related distinction without reproducing the proof. At δ=0 there are exact invisible directions. At small nonzero δ there can be full rank with very weak directions. No arithmetic lower bound protects those weak directions from measurement noise.

A determinant is a product of singular values. It does not, by itself, measure the weakest observable direction. The paper uses determinant smallness to derive a contradiction; an inverse problem needs quantitative control of the directions relevant to recovery.

Nor does an irrationality exponent of two provide a practical finite-precision memory guarantee. The manuscript's eventual approximation bound is not a uniform bound with a known usable constant at every denominator. Exact nonrepetition of a rotation also does not give infinitely many noise-distinguishable phases.

As checked on 6 November 2026, OpenAI's repository describes the work as output from an unnamed internal model. Its formalization notes report the π exponent statement and explicitly put Flint–Hills convergence outside the selected Lean statement [4]. This note does not independently certify that proof, identify the model, or establish specialist review status.

## 9. The connection to the flowers and beaches

The FourierFlower construction discussed alongside this experiment constrains its planar Fourier frequencies to `n_k = 1 + p m_k`. For integer `m_k`, that implies

$$z(u+2\pi/p)=e^{2\pi i/p}z(u).$$

Many harmonic details obey one shared petal symmetry. That is a useful reduction: the structure remains coordinated despite its complexity. This identity concerns the stated planar construction; additional 3D modifications need not preserve it.

In the graph, internal structure also produces a reduction, but now at the measurement interface: distinct states yield identical port readings. Symmetry can both organize a representation and hide distinctions from an observer. Whether that is helpful depends on the task.

The beaches add a crucial difference. Their stored variables may change the propagation operator itself. Identifying a material parameter from echoes is not the same problem as identifying the initial state x of a known L. One needs the response Jacobian with respect to those parameters, and possibly nonlinear or sequential identification methods. Friday's 15-dimensional state result cannot simply be copied into a beach's claimed memory capacity.

The conceptual connection is operational: **which histories remain distinguishable by the future experiments the system is allowed to perform?**

For a fixed linear system, the unobservable subspace describes state differences that all future outputs ignore. With interventions, the relevant equivalence becomes stricter: two states are equivalent only if every permitted future experiment gives the same observations. Changing the permitted operations can separate states that previously belonged to the same class.

That is an interpretation in terms of control and observation, not a new definition of intelligence. It gives a concrete way to avoid counting hidden physical degrees of freedom as usable memory when no allowed task can retrieve or exploit them.

## 10. What deserves to survive this discussion

The most useful conclusion is a requirement on a system, not a special constant:

> Preserve the distinctions the task needs, suppress the distinctions it can ignore, and provide an operation that exposes the needed distinctions at a measurable cost.

Friday provides a small exact example of the access part. It also provides two reasons to be demanding: the original external phase change cannot touch half the state space, and an internal change can restore full rank while leaving the weakest readings very small.

The broader ingredients are established. Novel utility would have to come from a working mechanism that uses them better under an explicit constraint. Nothing here establishes a general advantage for irrational phases, a learning algorithm, or a new law connecting π to cognition.

What is worth keeping is the move from resemblance to mechanism: an explicit hidden subspace, an explicit no-go result for a family of interventions, and an explicit change that opens access.

## Sources and provenance

1. **Magnetic gauge equivalence:** *Magnetic-Sparseness and Schrödinger Operators on Graphs*, Annales Henri Poincaré (2020). [DOI and article](https://doi.org/10.1007/s00023-020-00885-6). The article states unitary equivalence for magnetic potentials with the same cycle fluxes.
2. **Observability and symmetry:** Giuseppe Notarstefano and Gianfranco Parlangeli, *Controllability and Observability of Grid Graphs via Reduction and Symmetries*. [arXiv:1203.0129](https://arxiv.org/abs/1203.0129). Background for eigenspace-based observability and graph symmetry; not the source of the specific five-flower calculation above.
3. **External π manuscript:** OpenAI, *The irrationality exponent of π is 2*, September 24, 2026. [Manuscript and sources](https://github.com/openai/math/tree/main/preprints/The-irrationality-exponent-of-pi-is-2-September-24-2026).
4. **Reported formalization scope and provenance:** [π formalization notes](https://github.com/openai/math/blob/main/lean/docs/017.md); [repository description](https://github.com/openai/math/blob/main/README.md). Status statements above are dated, not promises about later revisions.
5. **Related projects:** [FourierFlower](https://github.com/anttiluode/FourierFlower), [SelfIndexingBeach](https://github.com/anttiluode/SelfIndexingBeach), and [ComputingBeachClaude](https://github.com/anttiluode/ComputingBeachClaude). Connections here refer to the constructions described in the supplied discussion, not a fresh audit of every linked repository.

The original experiment was supplied during Antti Luode's discussion with several assistants. It is preserved in [original/magnetic.py](../original/magnetic.py). The written analysis, closed-form basis, and reproducible audit were prepared with ChatGPT assistance. The original naming and labels are preserved only in the historical script; the audit distinguishes Laplacian eigenvalues from observation singular values.
