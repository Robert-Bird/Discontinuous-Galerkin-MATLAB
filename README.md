# Discontinuous-Galerkin-MATLAB

An MATLAB framework for **hp-adaptive Discontinuous Galerkin (DG)** finite
element analysis of **2D linear elasticity**, using the **Symmetric Interior Penalty
Galerkin (SIPG)** formulation on unstructured triangular meshes.

Given a problem (geometry, boundary conditions, material properties), the code repeatedly
solves the elasticity problem, estimates the error in each element, and refines the mesh -
either by subdividing elements (**h**-refinement) or by raising the local polynomial order
(**p**-refinement) - driven by two user thresholds, `delta1` and `delta2`.

## Features

- SIPG discontinuous Galerkin discretisation of 2D linear elasticity
- Hierarchical shape functions supporting variable polynomial order `p`
- `h`, `p`, and combined `hp` mesh adaptivity with polynomial-order smoothing across faces
- Element-wise error estimator plus `L2`- and `DG`-norm error measures
- Manufactured (analytical) solutions via the Symbolic Math Toolbox for verification and
  convergence studies
- A built-in, dependency-free Delaunay mesher (no external mesher required)
- A regression test harness with stored expected values

## Requirements

- **MATLAB** (R2018a or later recommended) - mesh generation uses only base MATLAB
  (`delaunayTriangulation`), so no meshing toolbox is required
- **Symbolic Math Toolbox** - used by the analytical problem generators to manufacture the
  body force/stress from a chosen displacement field

## Installation

From the repository root, run:

```matlab
path_add
```

This adds all of the framework's subfolders (`adaptivity/`, `stiffness_matrix/`,
`problem_set_up/`, …) to the MATLAB path. That's the only setup step - mesh generation is
built in.

## Mesh generation

The mesh is produced by
[`mesh_generation/triangle_mesh.m`](mesh_generation/triangle_mesh.m), a dependency-free
mesher built on MATLAB's `delaunayTriangulation` (base MATLAB - no toolboxes). It resamples
the boundary and fills the interior with a grid at the target spacing `hfun`, then
triangulates. Element size follows `hfun` (no quality/size refinement); the hp-adaptivity
loop refines further. Called from
[`Mesh_gen_square.m`](mesh_generation/square/Mesh_gen_square.m):

```matlab
[vert,conn,tria] = triangle_mesh(node,edge,hfun);   % hfun = target element size
```

## Quick start

With the path set up, run one of the worked examples or the analytical driver from the
repository root:

```matlab
path_add                % set up the path (do this once per session)
example_problem_1       % a worked example
analytical_problem      % manufactured-solution driver with convergence plots
```

To run the regression test suite:

```matlab
tests                   % or: main_tests
```

### Defining your own problem

Problems are defined in `problem_set_up/`. For a problem with a known (manufactured)
solution, edit
[`problem_set_up/analytical_problem/analytical_problem_generator.m`](problem_set_up/analytical_problem/analytical_problem_generator.m):

- `node` / `edge` - domain geometry (vertices and the connecting boundary edges)
- `BC` - boundary-condition flags per edge (Dirichlet / Neumann / mixed; see the header
  comment in that file for the flag table)
- `E`, `nu` - material properties (Young's modulus, Poisson's ratio)
- `u`, `v` - the symbolic displacement field; the generator differentiates it to produce
  the consistent body force and stress, writing `BodyForce.m`, `BodyStress.m`, and
  `ImposedDisplacements.m` automatically
- `d_1`, `d_2`, `loop_end`, `sim_end` - adaptivity strategy and number of refinement loops

## Repository layout

| Directory | Contents |
|-----------|----------|
| `stiffness_matrix/` | DG stiffness assembly: volume integrals, interior/boundary surface (penalty) integrals, and hierarchical shape functions |
| `adaptivity/` | `h`, `p`, and `hp` marking and mesh-refinement routines |
| `problem_set_up/` | Geometry, boundary conditions, material model, and (symbolic) manufactured solutions |
| `mesh_generation/` | The `triangle_mesh` mesher, mesh seeding, and face-topology construction |
| `rhs/` | Force-vector integration (body, Neumann, Dirichlet, mixed) |
| `solvers/` | Linear solve |
| `error_calculation/`, `norms/` | Error estimators and `L2`/`DG` norm computations |
| `boundary_conditions/` | Boundary-condition averaging helpers |
| `plotters/` | Mesh, displacement, stress and polynomial-order plotting |
| `tests/`, `example_problems/` | Test harness, expected values, and example meshes |

Top-level entry points: `path_add.m`, `analytical_problem.m`, `example_problem_1/2/3.m`,
`main_tests.m`, `tests.m`.
