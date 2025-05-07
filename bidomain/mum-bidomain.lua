
ug_load_script("ug_util.lua")
ug_load_script("util/refinement_util.lua")

--  MODEL PARAMETERS
Cm = 1e-11 --  1[μF cm^{-2}] = 1e-11 [mA ms mV^{-1] μm^{-2]]
-- Cg = 5e-12 --  0.5[μF cm^{-2}] = 5e-12 [ms Omega^{-1] µm^{-2}]
Rm = 5e11 --  5 [kOhm cm^2] = 5e11[Ohm μm^2] specific membrane resistance
Rg = 1.5e8 --  0.0015 [kOhm cm^2] = 1.5e^8 [Ohm μm^2] specific gap junction resistance

v_0 = -80 --  -80[mV] resting membrane potential
conductivity_sigma_i = 4e-7 --  4  [mS cm^{-1} = 4e-7 [Ohm^{-1} µm^{-1}]
conductivity_sigma_e = 2e-6 --  20 [mS cm^{-1} = 2e-6 [Ohm^{-1} µm^{-1}]
v_infty = v_0
ue_infty = 0.0 -- 0[mV]

number_cells_x = 25 --  Nx number of cells in x direction
number_cells_y = 25 --  Ny number of cells in y direction
number_cells_z = 1 --  Nz number of cells in z direction

cell_length_x = 120 --  lx cell length[μm] in x direction
cell_length_y = 12 --  ly cell length[μm] in y direction
cell_length_z = 12 --  lz cell length[μm] in z direction


--  volume intra- and extracellular space
total_number_cells = number_cells_x * number_cells_y * number_cells_z --  omega domain with 25x25 cells in 2D
domain_cell_length_x = number_cells_x * cell_length_x --  Lx length[μm] of omega domain in x direction
domain_cell_length_y = number_cells_y * cell_length_y --  Ly length[μm] of omega domain in y direction
domain_cell_length_z = number_cells_z * cell_length_z --  Lz length[μm] of omega domain in z direction
Vcell = cell_length_x * cell_length_y * cell_length_z --  volume[μm³] of one cell
print("Volume of one cell is: "..Vcell.."μm³")
Vomegaintra = total_number_cells * Vcell --  volume[μm³] of intracellular space
Vomega = Vomegaintra/0.75 --  volume of omega domain
Vomegaextra = Vomega - Vomegaintra --  volume of extracellular space
print("Volume of domain is: "..Vomega.."μm³")
print("Volume of intracellular space is".. Vomegaintra.."μm³")
print("Volume of extracellular space is".. Vomegaextra.."μm³")


--  compute surface to volume ratio
Acell = 2* cell_length_x * cell_length_y + 2* cell_length_y * cell_length_z + 2* cell_length_x * cell_length_z --  total area[μm²] of one cell
Agamaomega = Acell * total_number_cells --  total membrane area of all cells in the domain
chi = Agamaomega/Vomega
print("surface to volume ratio is: "..chi.." 1/μm")


--  total area of one intercalated disc
A_j_x = cell_length_y * cell_length_z * 98/100 --  area single intercalated disc perpendicular to the x direction
A_j_y = cell_length_x * cell_length_z * 80/100 --  area single intercalated disc perpendicular to the y direction
A_j_z = cell_length_x * cell_length_y * 80/100 --  area single intercalated disc perpendicular to the z direction
print("Area single intercalated disc perpendicular to the x direction: ".. A_j_x .."μm²")
print("Area single intercalated disc perpendicular to the y direction: ".. A_j_y .."μm²")
print("Area single intercalated disc perpendicular to the z direction: ".. A_j_z .."μm²")


--  total area of intercalated discs in omega domain
disc_x = number_cells_x - 1 --  num of intercalated discs in x direction
disc_y = number_cells_y - 1 --  num of intercalated discs in y direction
disc_z = number_cells_z - 1 --  num of intercalated discs in z direction

Ajomegax = number_cells_y * number_cells_z * disc_x * A_j_x --  total area of all intercalated discs perpendicular to the x direction
Ajomegay = number_cells_x * number_cells_z * disc_y * A_j_y --  total area of all intercalated discs perpendicular to the y direction
Ajomegaz = number_cells_x * number_cells_y * disc_z * A_j_z --  total area of all intercalated discs perpendicular to the y direction
print("Area intercalated discs perpendicular to the x direction: "..Ajomegax.."μm²")
print("Area intercalated discs perpendicular to the y direction: "..Ajomegay.."μm²")
print("Area intercalated discs perpendicular to the z direction: "..Ajomegaz.."μm²")


--  intracellular conductivity tensor
gxnumerator = conductivity_sigma_i * (number_cells_x -1) * (number_cells_x -1) * Rg * Vomegaintra
gxdenominator = domain_cell_length_x * domain_cell_length_x * Ajomegax
gx = gxnumerator / gxdenominator
gynumerator = conductivity_sigma_i * (number_cells_y -1) * (number_cells_y -1) * Rg * Vomegaintra
gydenominator = domain_cell_length_y * domain_cell_length_y * Ajomegay
gy = gynumerator / gydenominator
-- gznumerator = conductivity_sigma_i * (Nz-1) * (Nz-1) * Rg * Vomegaintra
-- gzdenominator = Lz * Lz * Ajomegaz
-- gz = gznumerator / gzdenominator
print("Value gx: "..gx)
print("Value gy: "..gy)
-- print("Value gz: {gz}")
conductivity_sigma_i_x = (Vomegaintra * conductivity_sigma_i) / (Vomega * (1 + gx))
conductivity_sigma_i_y = (Vomegaintra * conductivity_sigma_i) / (Vomega * (1 + gy))
-- conductivity_sigma_i_z = (Vomegaintra * conductivity_sigma_i) / (Vomega * (1 + gz))
print("Value conductivity_sigma_i_x: "..conductivity_sigma_i_x.." S/μm")
print("Value conductivity_sigma_i_y: "..conductivity_sigma_i_y.." S/μm")
-- print("Value conductivity_sigma_i_z: {conductivity_sigma_i_z} S/μm")


--  extracellular conductivity tensor
conductivity_sigma_e_x = Vomegaextra * conductivity_sigma_e / Vomega
conductivity_sigma_e_y = conductivity_sigma_e_x
conductivity_sigma_e_z = conductivity_sigma_e_x
print("Value conductivity_sigma_e_x: "..conductivity_sigma_e_x.." S/μm")
print("Value conductivity_sigma_e_y: "..conductivity_sigma_e_y.." S/μm")
print("Value conductivity_sigma_e_z: "..conductivity_sigma_e_z.." S/μm")


--  init UG with dim = 2 and CPU block algebra
dim = 2
InitUG(dim, AlgebraType("CPU", 1))
grid_name = "grids/mum_cells25x25_hv.ugx"
spatial_num_ref = 2


--  create domain (\Omega)
required_subsets = {"Inner", "V", "C", "H"}
domain = util.CreateDomain(grid_name, spatial_num_ref, required_subsets)


--  create approximation space (Ansatzraum)
approx_space = ApproximationSpace(domain)


--  add unknown functions to approximation space and define
--  interpolation method (Lagrange Interpolation first order)
approx_space:add_fct("V", "Lagrange", 1) -- add function membrane potential, V = Ui - Ue
approx_space:add_fct("Ue", "Lagrange", 1) -- add function extracellular potential
approx_space:init_levels()
approx_space:init_top_surface()
approx_space:print_statistic()




--  create element discretization (FV of convection diffusion equation)
--  we need one discretization per unknown, we can handle them
--  by putting them in a dictionary


-- <Test-Section>
local conductivity_sigma_i_x = 2e-7
local conductivity_sigma_i_y = 4e-8
local conductivity_sigma_e_x = 0.2*2.0e-6
local conductivity_sigma_e_y = 0.2*2.0e-6
local chi = 0.25 -- 0.25
local impuls = 4e-9 -- ø 4e-10
-- </Test Section>
print("conductivity ix: ".. conductivity_sigma_i_x)
print("conductivity iy: ".. conductivity_sigma_i_y)
print("conductivity ex: ".. conductivity_sigma_e_x)
print("conductivity ey: ".. conductivity_sigma_e_y)


chiCm = chi * Cm
print("chiCm =", chiCm)

inner_diffusion_tensor = ConstUserMatrix() -- create Intracellular Matrix Mi
inner_diffusion_tensor:set_entry(0, 0, conductivity_sigma_i_x)
inner_diffusion_tensor:set_entry(1, 1, conductivity_sigma_i_y)

extracellular_diffusion_tensor = ConstUserMatrix() -- create Extracellular Matrix Me
extracellular_diffusion_tensor:set_entry(0, 0, conductivity_sigma_e_x)
extracellular_diffusion_tensor:set_entry(1, 1, conductivity_sigma_e_y)

io_diffusion_tensor = ScaleAddLinkerMatrix() -- sum of extra and intracellular Matrix Mi+Me
io_diffusion_tensor:add(1, inner_diffusion_tensor)
io_diffusion_tensor:add(1, extracellular_diffusion_tensor)

minus_inner_diffusion_tensor = ScaleAddLinkerMatrix()
minus_inner_diffusion_tensor:add(-1.0, inner_diffusion_tensor)

minus_extracellular_diffusion_tensor = ScaleAddLinkerMatrix()
minus_extracellular_diffusion_tensor:add(-1.0, extracellular_diffusion_tensor)

minus_io_diffusion_tensor = ScaleAddLinkerMatrix()
minus_io_diffusion_tensor:add(-1, inner_diffusion_tensor)
minus_io_diffusion_tensor:add(-1, extracellular_diffusion_tensor)



--elem_disc_v = ConvectionDiffusionFV1("V", "Inner")
elem_disc_v = ConvectionDiffusionFV1("V", "Inner")
elem_disc_v:set_mass_scale(chiCm) -- m1 = chi*Cm
elem_disc_v:set_diffusion(inner_diffusion_tensor)
grad_v = elem_disc_v:gradient()


--elem_disc_ue = ConvectionDiffusionFV1("Ue", "Inner")
elem_disc_ue = ConvectionDiffusionFV1("Ue", "Inner")
elem_disc_ue:set_mass_scale(0.0)
elem_disc_ue:set_diffusion(io_diffusion_tensor)
elem_disc_ue:set_stationary()

grad_ue = elem_disc_ue:gradient()



flux_v = ScaleAddLinkerVectorMatrix()
flux_v:add(minus_inner_diffusion_tensor, grad_ue)
elem_disc_v:set_flux(flux_v)

flux_ue = ScaleAddLinkerVectorMatrix()
flux_ue:add(minus_inner_diffusion_tensor, grad_v)
elem_disc_ue:set_flux(flux_ue)







--  Ion Current I(V) = (V-V0)/Rm
--  define functions and derivatives describing the reaction term:

--  for V
function reaction_v(v)
    return chi*(v-v_0)/Rm
end

function reaction_v_dv(v)
    return chi/Rm
end

--  create Function Object compatible with UG4:
nonlinear_growth_v = LuaUserFunctionNumber(reaction_v, 1)
nonlinear_growth_v:set_input_and_deriv(0, elem_disc_v:value(),  reaction_v_dv)


--  add reaction term to discretization:
elem_disc_v:set_reaction(nonlinear_growth_v)


--  add source
function source_function(x, y, t, si)
    if t<=1.0 and x<=600.0 and y<=60.0 then
        return impuls
    else
        return 0.0
    end
end

elem_disc_v:set_source(source_function)

--  define initial values via functions (will be interpolated over domain!)
function initial_value_v(x, y, t, si)
    return v_infty
end

function initial_value_ue(x, y, t, si)
    return ue_infty
end

--  define boundary conditions:
dirichlet_bnd = DirichletBoundary()
-- dirichlet_bnd:add(0, "Ue", "H,V,C")
-- dirichlet_bnd:add(uinfty, "V", "Boundary")
--dirichlet_bnd:add(-80, "V", "V,C")

-- neumann_boundary = NeumannBoundaryFV1("V");
-- neumann_boundary:add(0, "V", "Boundary")
-- neumann_boundary:add(0, "V", "H,C")


--  add element discretization and boundary conditions to
--  the domain discretization
domain_disc = DomainDiscretization(approx_space)
domain_disc:add(elem_disc_v)
domain_disc:add(elem_disc_ue)
domain_disc:add(dirichlet_bnd)
-- domain_disc:add(neumann_boundary)

--  create gridFunction that stores an unknown variabel for each
--  function in our approximation space
u = GridFunction(approx_space)


--  create LU solver (direct solver -> be careful for simulations
--  with higher degrees of freedom, will get slow (O(n**3)))
--linear_base_solver = LU()
linear_base_solver = LU() --SuperLU()

--[[
local linear_conv_check = ConvCheck()
linear_conv_check:set_maximum_steps(50)
linear_conv_check:set_reduction(1e-8)
linear_conv_check:set_minimum_defect(1e-14)
linear_conv_check:set_supress_unsuccessful(true)

local smoother = Jacobi(0.5)

local	gmg = GeometricMultiGrid(approx_space)
gmg:set_discretization(domain_disc)
gmg:set_base_level(0)
gmg:set_base_solver(linear_base_solver)
gmg:set_presmoother(smoother)
gmg:set_postsmoother(smoother)
gmg:set_cycle_type("F")
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)
gmg:set_rap(true)

local	limex_linear_solver = BiCGStab()
limex_linear_solver:set_preconditioner(gmg)
limex_linear_solver:set_convergence_check(linear_conv_check)
--]]

--  interpolate initial conditions over domain
Interpolate("initial_value_v", u, "V")
Interpolate("initial_value_ue", u, "Ue")


--  LIMEX SOLVER CONFIGURATION
limex_estimator = CompositeGridFunctionEstimator()
l2_component_space_v = L2ComponentSpace("V", 2)
l2_component_space_ue = L2ComponentSpace("Ue", 2)
limex_estimator:add(l2_component_space_v)
limex_estimator:add(l2_component_space_ue)

limex_cost_strategy = LimexNonlinearCost()

limex_newton_conv_check = ConvCheck(1, 5e-10, 1e-10, true)
limex_newton_conv_check:set_supress_unsuccessful(true)

limex_newton_solver = NewtonSolver()
limex_newton_solver:set_linear_solver(linear_base_solver)
limex_newton_solver:set_convergence_check(limex_newton_conv_check)

limex_time_integrator = LimexTimeIntegrator(2)
limex_time_integrator:add_stage(1, limex_newton_solver, domain_disc)
limex_time_integrator:add_stage(2, limex_newton_solver, domain_disc)
limex_time_integrator:add_stage(3, limex_newton_solver, domain_disc)
limex_time_integrator:set_dt_min(1e-9)
limex_time_integrator:set_dt_max(0.5)
limex_time_integrator:add_error_estimator(limex_estimator)
limex_time_integrator:set_increase_factor(2.0)
limex_time_integrator:select_cost_strategy(limex_cost_strategy)
limex_time_integrator:enable_matrix_cache()
limex_time_integrator:set_time_step(0.1)
limex_time_integrator:set_tolerance(1e-3)



start_time = 0.0
end_time = 5.0




vtk = VTKOutput()
-- vtk:print("solution", u, 0, 0.0)
observer = VTKOutputObserver("output/solution",vtk)

limex_time_integrator:attach_observer(observer)
limex_time_integrator:apply(u, end_time, u, start_time)

--vtk:write_time_pvd()
