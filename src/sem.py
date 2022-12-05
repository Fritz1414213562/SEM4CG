import numpy as np
from fenics import *
#from mshr import *
from scipy.interpolate import RegularGridInterpolator
from sem_params import SEMParams
from gridspace import GridSpace
from const import SEMConst
import sys

class StericExclusionModel():

	TOL = 1e-14
	RELTOL = 1e-8
	MAXITR = 20000
	LOGLVL = 30

	def __init__(self, params : SEMParams):

		if SEMParams.log_flag:
			print("[1] Initialize box grids ...")
		self.space = GridSpace(params.box_params, params.grid_params)
		self.space.set_pore_params(*params.pore_params)
		self.bulk  = SEMConst.BULK
		self.cutoff = params.cutoff
		if SEMParams.log_flag:
			print("Done")

		if SEMParams.log_flag:
			print("[2] Initialize Krylov solver parameters ...")
		self.solver = KrylovSolver("gmres", "amg")
		#self.solver = KrylovSolver("gmres", "ilu")
		parameters["krylov_solver"]["nonzero_initial_guess"] = True
		parameters["krylov_solver"]["monitor_convergence"]   = True
		self.solver.parameters["relative_tolerance"] = self.RELTOL
		self.solver.parameters["maximum_iterations"] = self.MAXITR
		set_log_level(self.LOGLVL)
		if SEMParams.log_flag:
			print("Done")

		if SEMParams.log_flag:
			print("[3] Define Functional spaces and boundaries ...")
		self.mesh = self.__box_mesh(params)
		self.ground, self.terminal = self.__sub_domains(params)
		self.boundary_parts = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1)
		self.ground.mark(self.boundary_parts, 1)
		self.terminal.mark(self.boundary_parts, 2)
		self.ds = Measure("ds", domain = self.mesh, subdomain_data = self.boundary_parts)

		self.fspace_v = FunctionSpace(self.mesh,  'P', 1)
		self.fspace_f = FunctionSpace(self.mesh, 'CG', 1)
		self.bcs = [\
			DirichletBC(self.fspace_v, Constant(params.volt_bottom), self.ground),\
			DirichletBC(self.fspace_v, Constant(params.volt_top), self.terminal)\
		]
		self.u = TrialFunction(self.fspace_v)
		self.v = TestFunction(self.fspace_v)
		self.sig = Function(self.fspace_f)
		self.f = Constant(0)
		self.a = self.sig * dot(grad(self.u), grad(self.v)) * dx
		self.L = self.f * self.v * dx
		#self.DE = self.sig * dot(grad(self.u), grad(self.v)) * dx - self.f * self.v * dx
		#self.mat, self.L = lhs(self.DE), rhs(self.DE)
		self.sol = Function(self.fspace_v)
		self.flux_top = dot(Constant((0, 0, -1)), self.sig * nabla_grad(self.sol)) * self.ds(2)
		self.flux_btm = dot(Constant((0, 0, -1)), self.sig * nabla_grad(self.sol)) * self.ds(1)
		if SEMParams.log_flag:
			print("Done")

	#	self.flux_top = dot(Constant((1, 1, 1)), self.sig * nabla_grad(self.sol)) * ds(2)
	#	self.flux_btm = dot(Constant((1, 1, 1)), self.sig * nabla_grad(self.sol)) * ds(1)


	def run(self, positions, resnames):

		inbox_indices = \
			(positions[:, 2] < 0.5 * SEMParams.box_params[2]) & (positions[:, 2] > -0.5 * SEMParams.box_params[2]) &\
			(positions[:, 1] < 0.5 * SEMParams.box_params[1]) & (positions[:, 1] > -0.5 * SEMParams.box_params[1]) &\
			(positions[:, 0] < 0.5 * SEMParams.box_params[0]) & (positions[:, 0] > -0.5 * SEMParams.box_params[0])

		conds = self.space.compute_conductivity(self.cutoff, positions[inbox_indices], resnames[inbox_indices], self.bulk)
		self.intsig = RegularGridInterpolator((self.space.xs, self.space.ys, self.space.zs), conds, bounds_error = False, fill_value = self.bulk)

		self.__load()

		A, bb = assemble_system(self.a, self.L, self.bcs)
		self.solver.set_operator(A)
		self.solver.solve(self.sol.vector(), bb)
		vals = self.sol.vector().get_local()

		#A = assemble(self.a)
		#b = assemble(self.L)
		#[bc.apply(A, b) for bc in self.bcs]
		#u = Function(self.fspace_v)
		#U = self.u.vector()
		#self.solver.solve(A, U, b)

		#flux_top = dot(Constant((0, 0, 1)), self.sig * nabla_grad(self.sol)) * ds(2)
		#flux_btm = dot(Constant((0, 0, 1)), self.sig * nabla_grad(self.sol)) * ds(1)
		#ft = assemble(self.flux_top)
		#fb = assemble(self.flux_btm)

		ft = assemble(self.flux_top)
		fb = assemble(self.flux_btm)

		return ft, fb


	def __box_mesh(self, params: SEMParams):

		initials = - 0.5 * params.box_params
		finals   =   0.5 * params.box_params
		nums     = ((finals - initials) / params.grid_params).astype(int)

		return BoxMesh(Point(initials), Point(finals), nums[0], nums[1], nums[2])


	def __sub_domains(self, params: SEMParams):

		z_ini = - 0.5 * params.box_params[2]
		z_fin =   0.5 * params.box_params[2]

		ground   = CompiledSubDomain("on_boundary && near(x[2], zwall, tol)", tol = self.TOL, zwall = z_ini)
		terminal = CompiledSubDomain("on_boundary && near(x[2], zwall, tol)", tol = self.TOL, zwall = z_fin)

		return ground, terminal


	def __load(self):

		vec = self.sig.vector()
		values = vec.get_local()
		dofmap = self.fspace_f.dofmap()
		my_first, my_last = dofmap.ownership_range()

		fdim = self.fspace_f.dim()
		mdim = self.mesh.geometry().dim()
		F_dof_coords = self.fspace_f.tabulate_dof_coordinates()
		F_dof_coords.resize((fdim, mdim))

		unowned = dofmap.local_to_global_unowned()
		dofs = filter(lambda dof: dofmap.local_to_global_index(dof) not in unowned, range(my_last - my_first))
		coords = F_dof_coords[list(dofs)]
		values[:] = self.intsig(coords)
		vec.set_local(values)
		vec.apply("insert")
