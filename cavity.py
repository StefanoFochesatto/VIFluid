# Solve the time-dependent Navier-Stokes equations, using backward
# Euler time stepping, in a lid-driven cavity with a stress-free base.
# Velocity and pressure are approximated using Taylor-Hood (P_2 x P_1)
# elements. The implicit step equations are solved by Newton's method
# using a monolithic MUMPS direct solution for each Newton step.

# settings for classroom demonstration; this run takes a couple of
# minutes; change to m=32 and N=50 for quicker run, for example
import matplotlib.animation as animation
from animate import adapt, RiemannianMetric
from firedrake.pyplot import tripcolor
from firedrake.output import VTKFile
import matplotlib.pyplot as plt
from viamr import VIAMR


m = 15                # resolution; m x m mesh
#N = 200    # number of time steps
N = 200
dt = 0.1                  # time step
Re = 2000.0               # Reynolds number; Re -> 0 is very viscous
outname = 'result.pvd'    # writes here; open this with Paraview


from firedrake import *
from animatePVD import *
from navierstokes import *

mesh = UnitSquareMesh(m, m)
for j in range(N):
    if j == 0:
        # mesh, function spaces, functions
        Z, V, W, up = NSFunctions(mesh)
    else:
        Z, V, W, upNew = NSFunctions(mesh)
        upNew.interpolate(up)
        up = Function(Z).interpolate(upNew)
        

        

    # Saving this for VI formulation
    # Define ufl for some inequality constraint on the pressure
    ubMixed = Function(Z)
    ubVelocity, ubPressure= ubMixed.subfunctions
    ubPressure.assign(Function(W).interpolate(PETSc.INFINITY))
    ubVelocity.assign(Function(V).interpolate(
        as_vector([PETSc.INFINITY, PETSc.INFINITY])))

    lbMixed = Function(Z)
    lbVelocity, lbPressure = lbMixed.subfunctions
    lbPressure.assign(Function(W).interpolate((PETSc.NINFINITY)))
    lbVelocity.assign(Function(V).interpolate(
        as_vector([(PETSc.NINFINITY), (PETSc.NINFINITY)])))



    # weak form for an implicit time step

    # Cross mesh interpolation
    if j == 0:
        uold = Function(V)
        uold.interpolate(as_vector([0.0, 0.0]))
        #outfile = VTKFile(outname)
        u, p = up.subfunctions
        u.rename("u (velocity)")
        p.rename("p (pressure)")
        t = 0.0
    else:
        u, p = up.subfunctions
        u.rename("u (velocity)")
        p.rename("p (pressure)")
        uold = Function(V).interpolate(uold)    
    
    F = NSTimeStepWeakForm(Z, up, uold, dt=dt, Re=Re)
    sparams = NSSolverParameters()

    # Dirichlet conditions: moving top, no-slip sides
    # Neumann conditions:   no stress bottom
    # note there is no null space; the Jacobian is invertible
    bcs = [
        DirichletBC(Z.sub(0), Constant((1.0, 0.0)), (4,)),  # top
        DirichletBC(Z.sub(0), Constant((0.0, 0.0)), (1, 2)),  # sides
    ]

    # main time-stepping solve
    print(f'running {N} time steps of length {dt} to tf={N*dt} on {m}x{m} mesh,')
    print(f'  saving velocity and pressure at each step to {outname} ...')
    print(f't = {t:.3f}:')
    
    # Solve
    problem = NonlinearVariationalProblem(F, up, bcs)
    solver = NonlinearVariationalSolver(problem, nullspace=None,  solver_parameters=sparams, options_prefix="")
    solver.solve(bounds=(lbMixed, ubMixed))
    VTKFile(f"result/step{j}.pvd").write(u, p, time=t)
    
    
    # Update for next time step
    t += dt
    uold.interpolate(u)
    
    # Setup next mesh
    s = Function(W).interpolate(sqrt(uold.sub(0)**2 + uold.sub(1)**2))
    ags = Function(W).interpolate(sqrt(dot(grad(s), grad(s))))
    
    freeboundaryMetric = RiemannianMetric(TensorFunctionSpace(mesh, "CG", 1))
    mp = {
        "target_complexity": 300.0,  # target number of nodes
        "p": 2.0,  # normalisation order
        "h_min": 1.0e-7,  # minimum allowed edge length
        "h_max": 1.0,  # maximum allowed edge length
    }
    freeboundaryMetric.set_parameters({"dm_plex_metric": mp})
    freeboundaryMetric.interpolate(ags * ufl.Identity(2)) # Isotropic metric
    freeboundaryMetric.normalise()
    mesh = adapt(mesh, freeboundaryMetric)

print(f't = {t:.3f}:')
VTKFile(f"result/step{j}.pvd").write(u, p, time=200)
print(f'  ... done writing to {outname}')

animate_pvds("result.pvd","result", num_timesteps=N)


