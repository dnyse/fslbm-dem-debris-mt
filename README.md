# Numerical Simulation of Floating Debris Transport and Clogging

**Master's Thesis — Dennys Huber, Friedrich-Alexander-Universität Erlangen-Nürnberg, Lehrstuhl für Informatik 10 (Systemsimulation)**
Supervisor: Prof. Dr.-Ing. Harald Köstler · Advisor: M.Sc. Behzad Safaei

Floating debris that accumulates at bridge piers during floods raises the upstream water level and can cause a bridge to fail. This thesis tests whether a particle-resolved simulation — one that resolves every floater and the flow around it — can reproduce how such a jam forms from the interactions of the floaters with each other and with the flow. It builds on a coupled simulation in the **waLBerla** framework that combines a free-surface lattice Boltzmann method (FSLBM) with the discrete element method (DEM) and a momentum exchange method (MEM), extending it from fully submerged spheres to floating cylinders and cuboids, fixed pillars, and a measurement of the clog in the simulated flow. The headline results compare two synthetic obstruction layouts at the same scale — a **centred pillar** (a single-body geometric blockage) and a **hopper gap** (a narrow opening that drifting rods jam by themselves through multi-body interlocking) — plus a two-pillar layout where the rods form only a porous, non-jamming accumulation. The clogging results are purely numerical methodology demonstrations; they are not yet validated against flume experiments and do not predict the clogging of any particular real bridge.

<p align="center">
  <img src="apps/showcases/Antidunes/docs/thesis_readme_figures/twopillar_jam_3d.png" width="85%" alt="3D render of a debris jam forming at a two-pillar obstruction">
</p>
<p align="center"><em>Figure 1 — Simulated jam of floating rods, cuboids and spheres piling up against a two-pillar obstruction in the periodic open-channel flow (coupled FSLBM-DEM-MEM, Antidunes showcase). Colours denote particle shape/type.</em></p>

<p align="center">
  <img src="apps/showcases/Antidunes/docs/thesis_readme_figures/hopper_gap_jam_side.png" width="85%" alt="Side view of rods self-jamming in a narrow hopper gap, versus clear flow downstream">
</p>
<p align="center"><em>Figure 2 — The "hopper gap" layout: rods drifting in from upstream (left) jam a gap wide enough for a single rod, while the channel downstream of the obstruction (right) stays essentially clear.</em></p>

<p align="center">
  <img src="apps/showcases/Antidunes/docs/thesis_readme_figures/twopillar_backwater_side.png" width="85%" alt="Side profile showing the backwater rise upstream of a jammed pillar compared to the lower water level downstream">
</p>
<p align="center"><em>Figure 3 — Backwater Rise: the free surface sits visibly higher upstream of the jammed obstruction (left) than in the clear channel downstream (right), consistent with the discharge reduction reported in the thesis.</em></p>

## Thesis implementation

The thesis implementation lives in [`apps/showcases/Antidunes/`](apps/showcases/Antidunes/) — a waLBerla showcase that couples FSLBM (two-phase free-surface flow), DEM (granular particle dynamics) and MEM (fluid-particle momentum exchange) to simulate antidune bedforms and, for this thesis, floating-debris clogging at bridge-pier-like obstructions. See that directory's own [`CLAUDE.md`](apps/showcases/Antidunes/CLAUDE.md) for the full architecture notes, key source files and per-timestep execution order.

Build (from the waLBerla repo root; requires C++20, CMake ≥ 3.26, MPI, and Python ≥ 3.10 with `pystencils`/`lbmpy`/`jinja2` for LBM kernel code generation):

```bash
python -m venv .venv && source .venv/bin/activate
pip install "pystencils~=1.3" "lbmpy~=1.3" jinja2

cmake -G Ninja \
  -DWALBERLA_BUILD_WITH_CODEGEN=ON \
  -DWALBERLA_BUILD_SHOWCASES=ON \
  -DWALBERLA_BUILD_WITH_MPI=ON \
  -DPython_ROOT_DIR="$(pwd)/.venv/bin" \
  -B build/ .

cmake --build build/ --target Antidunes
cmake --build build/ --target BedGeneration
```

Full run instructions (MPI process-count constraints, the two-stage bed-generation → main-simulation workflow, experiment configurations E1–E4) are documented in the showcase's `CLAUDE.md`.

## Key References

- Schwarzmeier, C., Rettinger, C., Kemmler, S., Plewinski, J., Núñez-González, F., Köstler, H., Rüde, U., & Vowinckel, B. (2023). *Particle-resolved simulation of antidunes in free-surface flows*. Journal of Fluid Mechanics, 961, R1. https://doi.org/10.1017/jfm.2023.262 — the base FSLBM-DEM coupling that the Antidunes showcase (and this thesis) builds on.
- Kodam, M., Bharadwaj, R., Curtis, J., Hancock, B., & Wassgren, C. (2010). *Cylindrical object contact detection for use in discrete element method simulations. Part I — Contact detection algorithms*. Chemical Engineering Science, 65(22), 5852–5862. https://doi.org/10.1016/j.ces.2010.08.006
- Kodam, M., Bharadwaj, R., Curtis, J., Hancock, B., & Wassgren, C. (2010). *Cylindrical object contact detection for use in discrete element method simulations, Part II — Experimental validation*. Chemical Engineering Science, 65(22), 5863–5871. https://doi.org/10.1016/j.ces.2010.08.007
- De Cicco, P. N., Paris, E., Solari, L., & Ruiz-Villanueva, V. (2020). *Bridge pier shape influence on wood accumulation: Outcomes from flume experiments and numerical modelling*. Journal of Flood Risk Management, 13(2), e12599. https://doi.org/10.1111/jfr3.12599
- Erpicum, S., Poppema, D., Burghardt, L., Benet, L., Wüthrich, D., Klopries, E.-M., & Dewals, B. (2024). *A dataset of floating debris accumulation at bridges after July 2021 flood in Germany and Belgium*. Scientific Data, 11, 1092. https://doi.org/10.1038/s41597-024-03907-8
- Marin, A., & Souzy, M. (2025). *Clogging of Noncohesive Suspension Flows*. Annual Review of Fluid Mechanics, 57, 89–116. https://doi.org/10.1146/annurev-fluid-030124-112742

---

# waLBerla

waLBerla (widely applicable Lattice Boltzmann from Erlangen) is a massively parallel framework for multiphysics simulation applications.
Beyond computational fluid dynamics with the lattice Boltzmann method, the framework now features multiphase and free-surface flows,
rigid body and particle dynamics as well as fluid-structure coupling with moving geometries.
It scales from laptops to current and future supercomputers while maintaining
near-perfect efficiency.

## Getting Started with waLBerla

Refer to our [Setup Guide](https://walberla.pages.i10git.cs.fau.de/walberla/setup-instructions.html) for instructions
on setting up and building waLBerla.

## Documentation and Tutorials

You can find our framework documentation, guides, tutorials, and examples on the following pages:

- Latest Release: [C++ Framework](http://walberla.net/doxygen/index.html), [Python Interface](http://walberla.net/sphinx/index.html)
- Current Development Revision: [C++ Framework](walberla.pages.i10git.cs.fau.de/walberla/)

## Get involved

### Contributing

Please refer to the [contribution guide](CONTRIBUTING.md) for guidance on contributing to waLBerla.

### Support

To get in touch with the waLBerla developers, use our [Issue Tracker](https://i10git.cs.fau.de/walberla/walberla/issues)
or the waLBerla mailing list ([cs10-walberla@fau.de](mailto:cs10-walberla@fau.de)).

## Authors

Many thanks go to waLBerla's [contributors](AUTHORS.txt)

### Please cite us

If you use waLBerla in a publication, please cite the following articles:

Overview:

- M. Bauer et al., _waLBerla: A block-structured high-performance framework for
  multiphysics simulations_. Computers & Mathematics with Applications, 2020.
  https://doi.org/10.1016/j.camwa.2020.01.007.

Grid Refinement:

- F. Schornbaum and U. Rüde, _Massively parallel algorithms for the lattice Boltzmann
  method on nonuniform grids_. SIAM Journal on Scientific Computing, 2016.
  https://doi.org/10.1137/15M1035240

LBM - Particle Coupling:

- C. Rettinger and U. Rüde, _A comparative study of fluid-particle coupling methods for
  fully resolved lattice Boltzmann simulations_. Computers & Fluids, 2017.
  https://doi.org/10.1016/j.compfluid.2017.05.033

Free-surface LBM:

- C. Schwarzmeier et al., _Comparison of free-surface and conservative Allen-Cahn phase-field
  lattice Boltzmann method_. Journal of Computational Physics, 2023.
  https://doi.org/10.1016/j.jcp.2022.111753

Allen-Cahn phase-field LBM

- M. Holzer et al., _Highly efficient lattice Boltzmann multiphase simulations of immiscible
  fluids at high-density ratios on CPUs and GPUs through code generation_. The International Journal of High Performance Computing Applications, 2021.
  https://doi.org/10.1177/10943420211016525

MESA-PD:

- S. Eibl and U. Rüde, _A Modular and Extensible Software Architecture for Particle Dynamics_.
  Proceedings of the 8th International Conference on Discrete Element Methods.
  https://doi.org/10.48550/arXiv.1906.10963

Carbon Nanotubes:

- G. Drozdov et al., _Densification of single-walled carbon nanotube films:
  Mesoscopic distinct element method simulations and experimental validation_.
  Journal of Applied Physics, 2020. https://doi.org/10.1063/5.0025505

## License

waLBerla is licensed under [GPLv3](COPYING.txt).
