# AGENTS.md — PoissonGrids.jl

Working notes for coding agents. Read this before touching `src/`.

## What the package is

`PoissonGrids.jl` builds **one-dimensional adaptive grids** by equidistributing a
user-supplied scalar *monitor function* `M(x) > 0`. It is a small, dependency-light
mesh-generation package: one solver, three monitor constructors, ~350 lines of source.

Public API (all of it):

| Symbol | Kind | Role |
| --- | --- | --- |
| `solve_grid(xmin, xmax, M, nc; verbose, tol, maxiter)` | function | returns `nc+1` grid vertices spanning `[xmin, xmax]` |
| `gaussian_monitor(α, xc, σ)` | function | `M(x) = 1 + α exp(-(x-xc)²/σ²)` — localized bump |
| `tanh_monitor(α, κ, c; direction)` | function | monotone step from `1` to `1+α` across `c` |
| `window_monitor(α, κ, b, c)` | function | plateau `≈1+α` on `[c-b, c+b]`, `≈1` outside |
| `ConvergenceError` | exception | thrown when `maxiter` is exhausted before `tol` |

Monitors are plain closures. Anything callable, positive and finite on the grid works —
the three constructors are conveniences, not a required type.

## The math

Solve, on the computational coordinate ξ ∈ [0, 1],

```
d/dξ ( M(x(ξ)) · dx/dξ ) = 0 ,   x(0) = xmin,  x(1) = xmax
```

The first integral is `M · dx/dξ = const`, i.e. the **equidistribution principle**:
cell width is inversely proportional to `M`. Large `M` ⇒ dense cells. This is the 1-D
Winslow/Poisson grid generator, hence the package name.

There is no matrix assembly and no Newton solve. The steady state is reached by a
**damped pseudo-transient (second-order Richardson / "dynamic relaxation") iteration**
in fictitious time τ, with a Gershgorin-derived diagonal preconditioner. Both the
residual and its diagonal Jacobian entries come from the *same* local kernels, the
Jacobian path being ForwardDiff applied to the kernel — see `Val{AD}` below.

## Design picture

```
solve_grid(xmin, xmax, M, nc)
 │
 ├── validate_solver_inputs(...)                      fail fast on bad args
 ├── u ← collect(LinRange(xmin, xmax, nc+1))          physical vertices (the unknown)
 ├── Δξ ← 1/nc                                        uniform computational spacing
 ├── validate_monitor_on_grid(M, u)                   M finite & > 0 at every cell mid
 │
 ├── PRECONDITIONER PASS  (Val(true) == `derivative`)
 │     run_flux!(Gq, u, M, Δξ, derivative)            Gq.x[i] = Σ|∂q_i/∂u|
 │     Poisson!(G, q, Gq, b, Δξ, derivative)          G[i]    = Σ|∂r/∂q · ∂q/∂u|
 │     D .= G ./ 2                                    diagonal PC; α, β, Δτ from D
 │
 └── for iter in 1:maxiter                            (Val(false) == `primitive`)
       r0 ← r
       run_flux!(q, u, M, Δξ, primitive)              q_i   = -M(x_mid) (u_{i+1}-u_i)/Δξ
       Poisson!(r, q, b, Δξ, primitive)               r_{i+1} = -(q_E - q_W)/Δξ - b
       update_rate!(∂u∂τ, r, D, β)                    ∂u∂τ ← r/D + β·∂u∂τ   (momentum)
       update_variable!(u, ∂u∂τ, α)                   u[2:end-1] += α·∂u∂τ  (interior only)
       validate_monitor_on_grid(M, u)
       every 5 iters: residual = ‖r‖/√n; converged?   ← also refreshes D, Δτ, α, β
     end
```

Two structural ideas carry the whole file:

1. **`Val{AD}` kernel duality.** `run_flux!`/`qi`/`Poisson!` each take a
   `Val{true}`/`Val{false}` switch, aliased module-locally as `derivative` and
   `primitive`. `Val(false)` evaluates the residual; `Val(true)` runs
   `ForwardDiff.gradient` over the *same* local kernel (`qi_local`, `Poisson_local`)
   and sums `|∂|` to get a Gershgorin row-sum bound. One kernel, two meanings —
   the Jacobian can never drift out of sync with the residual.
2. **Kernels are `SVector`-local.** `Poisson_local` and `qi_local` see only a
   2-element stencil, so ForwardDiff is cheap, stack-allocated and non-allocating.

Boundary handling is implicit: `Poisson!` writes `u[I+1]` for `I ∈ 1:nv-2`, i.e. only
interior entries, and `update_variable!` steps only `u[2:end-1]`. Endpoints therefore
stay pinned at `xmin`/`xmax` by construction, with no explicit BC code.

## File map

| File | Contents |
| --- | --- |
| `src/PoissonGrids.jl` | module, exports, the `derivative`/`primitive` `Val` aliases |
| `src/monitors.jl` | the three monitor constructors — pure closures, no state |
| `src/solver.jl` | `solve_grid`, validation, `ConvergenceError`, and the physics kernels |
| `test/runtests.jl` | monitor algebra, solver invariants, input/monitor rejection, non-convergence |
| `docs/` | Documenter site; `index.md` narrative + `api.md` `@docs` blocks |

## Working in this repo

- Use the Julia MCP session with Revise; edits under `src/` reload automatically.
- Manifests are untracked and gitignored; `Pkg` regenerates them on demand.
- Full run: `Pkg.test()` — only before opening a PR. The suite is a few seconds.
- CI: `ci.yml` — a matrix of Julia 1.10 / 1.11 / 1.12 / 1 on Linux, macOS and Windows with
  coverage upload, plus a `downgrade` job pinning every dependency to the lowest bound
  `[compat]` declares. Also `docs.yml` and `TagBot.yml`.
- `julia` is the current release; `julia +lts` is 1.10, the lowest supported version.
  Pass `--startup-file=no` when invoking `+lts` from the shell.
- Docs: `julia --project=docs docs/make.jl` after `Pkg.develop(path=".")` + `Pkg.instantiate()`
  in that environment.
- **Every example in `README.md`, `docs/src/index.md` and the docstrings is a `jldoctest`**
  and is checked on every docs build (`makedocs` sets `checkdocs = :exports`, and a
  mismatch terminates the build). Changing the solver numerics — tolerances, the time-step
  rule, the preconditioner — will change printed cell widths and ratios and break those
  doctests. That is the intended alarm, not an unrelated failure: re-verify the new values
  in a session and update the expected output.

## Behavior worth knowing before you change anything

- **`nc` is cells, not vertices.** Output length is `nc + 1`. Off-by-one bugs live here.
- **The convergence check runs only every 5th iteration** (`iter == 1 || mod(iter,5) == 0`),
  and so does the preconditioner refresh. `maxiter` is therefore quantized in practice.
- Typical cost: the Gaussian example converges in a few hundred iterations at `nc = 32`,
  under ~800 at `nc = 128`. Iteration counts grow with monitor sharpness, not just `nc`.
- `validate_monitor_on_grid` runs **every** iteration — `nc` extra monitor evaluations
  per step. Deliberate fail-fast; do not silently drop it to buy speed.
- A constant monitor reproduces the uniform grid exactly (tested).
- Arithmetic is Float64 internally (`zeros(nv)`); a `Float32` domain still returns
  `Float32` because `update_variable!` writes back into `u`, but the iteration itself
  is not eltype-generic. Nothing here is GPU- or OffsetArray-ready.

## Known rough edges

These are real, verified in-session, and unowned — flag before "fixing" one, since each
is a behavior change rather than a bug report from a user.

1. **`λmax` is structurally constant.** `D .= G ./ 2` makes `G ./ D ≡ 2`, so
   `λmax = maximum(G ./ D)` is always exactly `2` and `Δτ = 2/√2 · CFL` never varies.
   The Gershgorin bound `G` still matters as the preconditioner `D`, but the eigenvalue
   estimate it feeds is a no-op. Any rework of the time-step logic starts here.
2. **`b` (the source term) is allocated and never written.** It is a zero vector for the
   whole solve; the `∂r∂b·∂b∂u` term named in the comment above `Poisson!` is not
   computed. It is a placeholder for a forcing term the solver does not yet support.
3. **`q` and `Gq` are single-field NamedTuples** `(; x = zeros(nc))`. That shape is
   scaffolding for a 2-D `(x, y)` extension; in 1-D it buys nothing.
4. **Residual normalization is odd:** `n = √(prod(nc))`, then `residual = ‖r‖/√n`, i.e.
   `‖r‖ / nc^(1/4)`. Whatever was intended, `tol` means something resolution-dependent.
5. `ConvergenceError.residual` is hardcoded `::Float64`.
