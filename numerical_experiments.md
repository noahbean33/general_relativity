# Numerical Demonstrations & Experiments for the GR Lecture Series

Comprehensive list of demonstrations that could be performed using
**EinsteinPy** (symbolic/semi-analytic Python) and/or the
**Einstein Toolkit** (full numerical relativity) to complement each lecture.

---

## Lecture 8 – Parallel Transport & Curvature

### 1. Parallel transport on the round sphere
- **Tool**: EinsteinPy / custom Python
- **Description**: Numerically parallel-transport a vector around a closed loop
  (e.g., a spherical triangle) on S² with the round metric and visualise the
  resulting holonomy (angle deficit). Lecture 8 derives the connection
  coefficients for the round sphere (`Γ¹₂₂ = −sin θ cos θ`,
  `Γ²₁₂ = Γ²₂₁ = cot θ`) — use them to integrate the parallel-transport ODE
  and show the rotation angle equals the solid angle enclosed.

### 2. Autoparallel (geodesic) curves on S²
- **Tool**: EinsteinPy (`Geodesic` module)
- **Description**: Integrate the autoparallel equation from Lecture 8
  (`γ̈ᵍ + Γᵍₘₙ γ̇ᵐ γ̇ⁿ = 0`) on the round sphere. Show that great circles
  are the solutions. Compare equatorial geodesics (θ = π/2) with
  non-equatorial ones and visualise on a 3-D sphere.

---

## Lecture 9 – Newtonian Spacetime is Curved

### 3. Newtonian gravity as spacetime curvature
- **Tool**: Custom Python / SymPy
- **Description**: Implement the Newtonian spacetime connection
  (`Γᵅ₀₀ = −fᵅ` with all other Γ = 0) for a point-mass gravitational field
  `f = −GM/r²r̂`. Compute R^α_{0β0} (tidal tensor) numerically, verify
  R₀₀ = 4πGρ (Poisson equation), and plot tidal forces for orbits around
  Earth.

### 4. Geodesic deviation in Newtonian spacetime
- **Tool**: Custom Python
- **Description**: Integrate two nearby worldlines in the Newtonian spacetime
  connection and show the deviation matches the Newtonian tidal acceleration.
  Visualise the stretching/squeezing of a cloud of test particles (tidal
  deformation).

---

## Lecture 10 – Metric Manifolds

### 5. Christoffel symbols and geodesics for common metrics
- **Tool**: EinsteinPy (`ChristoffelSymbols`, `MetricTensor`)
- **Description**: For the round-sphere metric
  `g = R²dθ² + R²sin²θ dφ²`, compute all Christoffel symbols symbolically,
  verify the geodesic equation, and compute the Riemann tensor, Ricci tensor,
  and scalar curvature `R = 2/R²`.

### 6. Length of curves on the sphere
- **Tool**: EinsteinPy / SciPy
- **Description**: Numerically compute the length functional
  `L[γ] = ∫√(g(v,v)) dλ` for various parametrisations of the equator (as
  done in Lecture 10 with `φ(λ)=2πλ³`) and confirm reparametrisation
  invariance.

### 7. Geodesic equation from Euler-Lagrange equations
- **Tool**: SymPy + numerical integration
- **Description**: Derive the geodesic equation for a general diagonal metric
  via the Euler-Lagrange equations (as done in Lecture 10), then numerically
  integrate for the Schwarzschild, FLRW, and sphere metrics. Compare geodesics
  obtained from the Lagrangian approach vs. the Christoffel-symbol approach.

---

## Lecture 11 – Symmetry

### 8. Killing vector fields and Lie derivatives
- **Tool**: EinsteinPy / SageMath
- **Description**: For the round sphere, compute the Lie derivative
  `L_X g` for the three SO(3) generators
  (`X₁ = −sin φ ∂_θ − cot θ cos φ ∂_φ`, etc.) from Lecture 11 and verify
  `L_X g = 0` numerically at sample points. Repeat for the Schwarzschild
  metric with `K_t = ∂/∂t` and `K_φ = ∂/∂φ`.

### 9. Induced metric via pull-back (ellipsoid)
- **Tool**: EinsteinPy / SymPy
- **Description**: As in Tutorial 11, compute the pull-back of the Euclidean
  metric onto an ellipsoid embedded in R³ for various (a, b, c). Visualise how
  the intrinsic geometry (Gaussian curvature) changes with the aspect ratios.

---

## Lecture 12 – Integration

### 10. Volume element and integration on curved manifolds
- **Tool**: SciPy / custom Python
- **Description**: Numerically integrate `∫√det(g) d^n x` for the round sphere
  (verify 4πR²), for the Schwarzschild spatial slices, and for FLRW spatial
  slices. Demonstrate how the volume element `√g` encodes curvature effects.

---

## Lecture 13 – Relativistic Spacetime

### 11. Proper time and the twin paradox
- **Tool**: Custom Python / EinsteinPy
- **Description**: In Minkowski spacetime with metric `η = diag(1,−1,−1,−1)`,
  numerically compute proper time `τ = ∫√(g(v,v)) dλ` for the two observers
  from Lecture 13 — the stationary observer `γ` and the travelling observer `δ`
  with speed parameter α. Plot `τ_δ = √(1−α²)` vs α and demonstrate time
  dilation.

### 12. Light cones and causal structure
- **Tool**: Custom Python (matplotlib)
- **Description**: For both Minkowski and Schwarzschild spacetimes, plot the
  light-cone structure (null geodesics) at various points. Show how the cones
  tilt as r → 2m in Schwarzschild coordinates.

---

## Lecture 14 – Matter

### 13. Charged particle in electromagnetic field (Lorentz force)
- **Tool**: Custom Python / SciPy
- **Description**: Integrate the equation of motion
  `m(∇_v v)^a = −q F^a_m γ̇^m` from Lecture 14 in a Schwarzschild background
  with an external electromagnetic field. Compare with flat-spacetime Lorentz
  force orbits.

### 14. Energy-momentum tensor of electromagnetic field
- **Tool**: EinsteinPy / SymPy
- **Description**: For a given electromagnetic potential A_a, compute
  `F_ab = 2∂[a A_b]` and the Maxwell energy-momentum tensor
  `T_ab = F_am F_bn g^mn − ¼F_mn F^mn g_ab`. Verify T^a_a = 0 (tracelessness)
  and compute energy density `T(e₀,e₀) = E² + B²` for specific field
  configurations.

### 15. Perfect fluid energy-momentum tensor
- **Tool**: EinsteinPy
- **Description**: Construct the perfect-fluid T^ab = (ρ+p)u^a u^b − p g^ab
  for a cosmological (FLRW) background. Verify `T^ab g_ab = 0` implies
  `p = ρ/3` (radiation). Explore the equation of state for dust (p=0) and
  dark energy (p = −ρ).

---

## Lecture 15 – Einstein Gravity

### 16. Verify known exact solutions of the Einstein equations
- **Tool**: EinsteinPy (`RicciTensor`, `EinsteinTensor`)
- **Description**: As suggested in Lecture 15's tutorial section, compute the
  Einstein tensor G_ab for the following metrics and verify they satisfy
  Einstein's equations `G_ab = 8πG T_ab`:
  - **Schwarzschild** (vacuum, T=0)
  - **Friedmann-Robertson-Walker (FRW)** (perfect fluid)
  - **pp-wave** (vacuum)
  - **Reissner-Nordström** (electromagnetic T_ab)

### 17. Hilbert action variation demo
- **Tool**: SymPy
- **Description**: For a simple 2-D metric (e.g., sphere), symbolically compute
  `δS_Hilbert = δ∫√g R` and show the resulting Euler-Lagrange equation yields
  the Einstein tensor `G^mn = R^mn − ½g^mn R`.

### 18. Cosmological constant and dark energy
- **Tool**: EinsteinPy + SciPy
- **Description**: Modify the Einstein equations with Λ ≠ 0 and compute the
  resulting Friedmann equations. Numerically integrate the scale factor a(t)
  for Λ < 0 (contracting), Λ = 0 (matter-dominated), and Λ > 0 (accelerating
  expansion). Plot the Hubble parameter H(t) in each case.

---

## Lecture 18 – Canonical Formulation of GR (3+1 Decomposition)

### 19. ADM decomposition
- **Tool**: Einstein Toolkit (Cactus/Carpet/McLachlan)
- **Description**: Perform a 3+1 decomposition of the Schwarzschild spacetime.
  Compute the spatial metric γ_ij, lapse N, shift β^i, and extrinsic curvature
  K_ij. This is the foundation for numerical relativity and directly connects
  to Lecture 18's goal of "integrating Einstein's equations by numerical codes."

### 20. Initial value problem for GR
- **Tool**: Einstein Toolkit
- **Description**: Set up constraint-satisfying initial data (Hamiltonian and
  momentum constraints) for:
  - A single Schwarzschild black hole (Brill-Lindquist data)
  - Brill wave initial data
  Evolve forward in time using the BSSN formulation and monitor constraint
  violations.

### 21. Gravitational wave extraction
- **Tool**: Einstein Toolkit (WeylScal4 + Multipole thorns)
- **Description**: Evolve a perturbed black hole and extract the Weyl scalar
  Ψ₄, decompose into spin-weighted spherical harmonics, and plot the
  gravitational waveform h₊, h× vs. time. Connects to the satellite lectures
  on gravitational waves (Bernard Schutz).

---

## Lecture 22 – Black Holes

### 22. Schwarzschild geodesics (radial & orbital)
- **Tool**: EinsteinPy (`Schwarzschild` geodesic classes)
- **Description**: Integrate timelike and null geodesics in the Schwarzschild
  metric. Reproduce:
  - Radial infall (proper time vs. coordinate time)
  - Circular orbits and ISCO (innermost stable circular orbit)
  - Precessing elliptical orbits (Mercury perihelion advance)
  - Photon sphere and light bending

### 23. Eddington-Finkelstein coordinates and horizon penetration
- **Tool**: EinsteinPy / custom Python
- **Description**: As directly motivated by Lecture 22, transform the
  Schwarzschild metric to Eddington-Finkelstein coordinates
  (`t̄ = t + 2m ln|r−2m|`), compute the metric in the new chart (as done with
  SageManifolds in the lecture notes), and plot ingoing/outgoing null geodesics.
  Show that ingoing geodesics smoothly cross r = 2m while outgoing ones cannot
  escape.

### 24. Kruskal-Szekeres maximal extension
- **Tool**: Custom Python
- **Description**: Plot the full Kruskal diagram of Schwarzschild spacetime.
  Overlay radial null geodesics, r = const curves, t = const curves, and the
  singularity r = 0. Visualise the four regions (exterior, black hole, white
  hole, parallel universe).

### 25. Black hole shadow
- **Tool**: EinsteinPy (`Shadow` module)
- **Description**: Ray-trace null geodesics around a Schwarzschild (and
  optionally Kerr) black hole to compute the apparent shadow as seen by a
  distant observer. Compare with the Event Horizon Telescope observations.

---

## Full Numerical Relativity (Einstein Toolkit)

### 26. Binary black hole inspiral and merger
- **Tool**: Einstein Toolkit (McLachlan + Carpet AMR)
- **Description**: The flagship NR experiment. Set up Bowen-York initial data
  for two orbiting black holes, evolve with BSSN or CCZ4, extract
  gravitational waveforms, and compute the final black hole mass and spin.
  Directly realises Lecture 18's purpose item 2: "integrate Einstein's
  equations by numerical codes."

### 27. Gravitational collapse to a black hole
- **Tool**: Einstein Toolkit (GRHydro)
- **Description**: Simulate the collapse of a spherically symmetric perfect
  fluid star (Oppenheimer-Snyder collapse). Monitor apparent horizon formation.
  Connects Lecture 14 (matter/perfect fluid) with Lecture 22 (black holes).

### 28. Neutron star oscillations
- **Tool**: Einstein Toolkit (GRHydro + TOV initial data)
- **Description**: Set up a TOV (Tolman-Oppenheimer-Volkoff) equilibrium
  neutron star, perturb it, and evolve. Extract quasi-normal mode frequencies.
  Demonstrates the full pipeline from matter (Lecture 14) through Einstein's
  equations (Lecture 15) to observable oscillation modes.

### 29. Cosmological (FLRW) evolution
- **Tool**: Einstein Toolkit (FLRWSolver) or custom Python
- **Description**: Evolve the Friedmann equations numerically for different
  matter content (dust, radiation, cosmological constant). Plot a(t), H(t),
  deceleration parameter q(t). Verify the analytic solutions for matter-
  dominated, radiation-dominated, and Λ-dominated universes.

---

## Cross-Cutting / Foundational Demonstrations

### 30. Riemann tensor symmetries and Bianchi identity
- **Tool**: EinsteinPy / SymPy
- **Description**: For several metrics (sphere, Schwarzschild, FRW), compute
  all components of R_abcd and verify:
  - Antisymmetry: R_abcd = −R_bacd = −R_abdc
  - Pair symmetry: R_abcd = R_cdab
  - First Bianchi identity: R_a[bcd] = 0
  - Contracted Bianchi identity: ∇_a G^ab = 0

### 31. Coordinate transformations and chart transitions
- **Tool**: EinsteinPy / SymPy
- **Description**: Implement the Rⁿ coordinate systems from `Rn.sage`
  (Cartesian, spherical, cylindrical for R², R³, R⁴) in pure Python using
  EinsteinPy. Verify transition maps, compute metrics in all coordinate
  systems, and construct orthonormal frames.

### 32. Geodesic deviation equation
- **Tool**: SciPy + EinsteinPy
- **Description**: Numerically integrate the geodesic deviation equation
  `D²ξᵃ/dτ² = −R^a_{bcd} u^b ξ^c u^d` in Schwarzschild spacetime. Visualise
  how a bundle of nearby geodesics diverges/converges, demonstrating the
  geometric meaning of the Riemann tensor from Lecture 8.

### 33. Gravitational lensing
- **Tool**: EinsteinPy + custom ray-tracing
- **Description**: Trace null geodesics in the Schwarzschild metric to compute
  the deflection angle of light passing near a massive object. Compare
  numerical result with the weak-field analytic formula `Δφ = 4GM/c²b`.
  Connects to the satellite lecture by Marcus C. Werner on gravitational
  lensing mentioned in Lecture 15.

### 34. Gravitational redshift
- **Tool**: EinsteinPy / custom Python
- **Description**: Compute the ratio of proper times for two stationary
  observers at different radii in Schwarzschild spacetime. Show
  `ν_received/ν_emitted = √(g_tt(r_emit)/g_tt(r_recv))` and plot as a
  function of radius, demonstrating the divergence at r = 2m.

---

## Summary by Tool

| Tool | Experiments |
|------|------------|
| **EinsteinPy** (symbolic + geodesics) | 1-2, 5-7, 8-9, 11-12, 14-16, 22-25, 30-34 |
| **Einstein Toolkit** (full NR) | 19-21, 26-29 |
| **Custom Python** (SymPy/SciPy/matplotlib) | 3-4, 6-7, 10-13, 17-18, 23-24, 29, 32-34 |

## Summary by Lecture

| Lecture | Experiments |
|---------|------------|
| L8: Parallel Transport & Curvature | 1, 2, 32 |
| L9: Newtonian Spacetime | 3, 4 |
| L10: Metric Manifolds | 5, 6, 7 |
| L11: Symmetry | 8, 9 |
| L12: Integration | 10 |
| L13: Relativistic Spacetime | 11, 12 |
| L14: Matter | 13, 14, 15 |
| L15: Einstein Gravity | 16, 17, 18, 30 |
| L18: Canonical Formulation | 19, 20, 21 |
| L22: Black Holes | 22, 23, 24, 25 |
| Full NR / Cross-cutting | 26, 27, 28, 29, 31, 33, 34 |
