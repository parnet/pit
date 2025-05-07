-- Copyright (c) 2010-2016:  G-CSC, Goethe University Frankfurt
-- Authors: Arne Naegel
-- 
-- This file is part of UG4.
-- 
-- UG4 is free software: you can redistribute it and/or modify it under the
-- terms of the GNU Lesser General Public License version 3 (as published by the
-- Free Software Foundation) with the following additional attribution
-- requirements (according to LGPL/GPL v3 §7):
-- 
-- (1) The following notice must be displayed in the Appropriate Legal Notices
-- of covered and combined works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (2) The following notice must be displayed at a prominent place in the
-- terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
-- 
-- (3) The following bibliography is recommended for citation and must be
-- preserved in all covered files:
-- "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
--   parallel geometric multigrid solver on hierarchically distributed grids.
--   Computing and visualization in science 16, 4 (2013), 151-164"
-- "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
--   flexible software system for simulating pde based models on high performance
--   computers. Computing and visualization in science 16, 4 (2013), 165-179"
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
-- GNU Lesser General Public License for more details.


PluginRequired("ConvectionDiffusion")
PluginRequired("Limex")
PluginRequired("XBraidForUG4")
ug_load_script("ug_util.lua")
print()
print()
print()


-- ===============================================================================================
-- get parameter from terminal
-- ===============================================================================================
pit_desc = {}
pit_inst = {}

space_time_communicator = util.xbraid.split_communicator(pit_desc,pit_inst)
print("")
ug_load_script("util/refinement_util.lua")

script_method = util.GetParam("--method", "XBraid", "target method of this script")
print("<script method='"..script_method.."'/>")

-- Parse parameters and print help
local grid_name = util.GetParam("--grid", "grids/monodomain_laplace_sample_grid_2d.ugx","filename of underlying grid")
local num_refs = util.GetParamNumber("--numRefs", 4, "number of refinements")

util.CheckAndPrintHelp("Beispiel: Mono-Domain-Gleichungen (2D)");

-- ===============================================================================================
-- Create geometric grid and refine
-- ===============================================================================================

local blockSize = 2 -- select 2 for 2x2 point-block system
InitUG(2, AlgebraType("CPU", blockSize));  


-- Load a domain without initial refinements.
local mandatory_subsets = { "Inner", "Boundary"}
dom = util.CreateDomain(grid_name, 0, mandatory_subsets)

-- Refine the domain (redistribution is handled internally for parallel runs)
print("refining...")
util.refinement.CreateRegularHierarchy(dom, num_refs, true)

-- ===============================================================================================
-- Define problem parameter
-- ===============================================================================================


-----------------------------------------
-- A) Modellparameter
-- 
-- \chi (Cm \frac{\partial V}{\partial t} + I)  + \nabla \cdot [ -sigma \nabla V] = 0
-- Transformation: u = (V-Vr)/(Vp-Vr)   <->  V = (Vp-Vr)*u + Vr
-- 
-- Einheiten:
-- 
-- L = 1 cm
-- T = 1 ms
-- U = 1 mV
-- 
-----------------------------------------

local Vr = -85.0   -- [mV]
local Vs = -75.0   -- [mV]
local Vp =  15.0   -- [mV]

local scaleV = (Vp-Vr)        -- scaling [mV]
local a = (Vs-Vr)/(Vp-Vr)     -- threshold [1]


local VoltPerMilliVolt = 1e-3
local Cm  = 1.0 * VoltPerMilliVolt --  [µF*cm^{-2}] = [ms * mA * V^{-1}]

local c1 = 0.175 -- [ms^{-1}]
local c2 = 0.03  -- [ms^{-1}]
local c3 = 0.011 -- [ms^{-1}]
local  b = 0.55  -- []


local sigma = 2.0  * VoltPerMilliVolt  -- mS/cm, 1 [S = A/V => 1 mS = 1 mA / V]
local chi = 2000 -- cm^{-1}

local DDiff = sigma/(Cm*chi)     --    mA/(V*cm) * V/(mA ms) / cm^{-1} = cm^2 / ms 
local lchar = 1.0 -- characteristic length scale for diffusion
local tDiff = lchar*lchar/DDiff


print ("Activation a:"..a)
print ("ScaleV :"..scaleV)


print ("DDiff = "..DDiff.." cm^2/ms")
print ("tDiff = "..tDiff.." ms")
print ("tReact = "..1/c1.." ms")
print ("tReact = "..1/c2.." ms")
print ("tReact = "..1/b.." ms")
print ("tReactDamp = "..1/(a*c1+c3*b).." ms")

local uinfty = 0.0
local winfty = uinfty / b


-- ===============================================================================================
-- Discretization of model and approximation space and
-- ===============================================================================================

local approx_space = ApproximationSpace(dom)

approx_space:add_fct("u", "Lagrange", 1)   -- linear approximation space
approx_space:add_fct("w", "Lagrange", 1)

print("approximation space:")
approx_space:init_levels()
approx_space:init_top_surface()
approx_space:print_statistic()


-----------------------------------------
-- C) Elementdiskretisierung (FV)
-----------------------------------------
local elem_disc = {}
elem_disc_u = ConvectionDiffusion("u", "Inner", "fv1")  -- potential
elem_disc_w = ConvectionDiffusion("w", "Inner", "fv1")  -- gating variable


-- elem_disc_u:set_mass_scale(1.0)  -- default
-- elem_disc_w:set_mass_scale(1.0)  -- default

-----------------------------------------
-- C).1 diffusion-tensor
-- Df = sigma/(Cm*chi)  
-----------------------------------------
local diffTensorU = ConstUserMatrix(DDiff/4.0)
diffTensorU:set_entry(0,0, DDiff)


elem_disc_u:set_diffusion(diffTensorU)
elem_disc_w:set_diffusion(0.0)

-----------------------------------------
-- C).2 reaction terms:
-- \dot u = (c1*u*(u-a)*(1.0-u) - c2*w)
-- \dot w = b*(v - d*w)
-----------------------------------------


function reaction_u(u,w)
	return (c1*u*(u-a)*(1.0-u) - c2*w) * (-1.0)
end

function reaction_w(u,w)
	return c3*(u - b*w)*(-1.0)
end

function reaction_u_du(u,w)
	return (a*(2.0*u-1.0) + (2.0-3.0*u)*u ) * (-c1)
end

function reaction_u_dw(u,w)
	return -c2 * (-1.0) --  todo simplify to c2
end

function reaction_w_du(v,w)
	return c3*(-1.0) --  todo simplify to -1.0*c3
end

function reaction_w_dw(u,w)
	return -c3*b*(-1.0)--  todo simplify c3*b
end

nonlinear_growth_u = LuaUserFunctionNumber("reaction_u", 2)
nonlinear_growth_u:set_input(0, elem_disc_u:value())
nonlinear_growth_u:set_input(1, elem_disc_w:value())
nonlinear_growth_u:set_deriv(0, "reaction_u_du")
nonlinear_growth_u:set_deriv(1, "reaction_u_dw")


nonlinear_growth_w = LuaUserFunctionNumber("reaction_w", 2)
nonlinear_growth_w:set_input(0, elem_disc_u:value())
nonlinear_growth_w:set_input(1, elem_disc_w:value())
nonlinear_growth_w:set_deriv(0, "reaction_w_du")
nonlinear_growth_w:set_deriv(1, "reaction_w_dw")


elem_disc_u:set_reaction(nonlinear_growth_u)
elem_disc_w:set_reaction(nonlinear_growth_w)

-----------------------------------------
--  C).3 Quellen = Externe Ströme (nicht benutzt!)
-- 1.0 / (Cm*scaleV)
-----------------------------------------
-- function ISourceU(x,y,t,si)  return 0.0/(Cm*scaleV) end
-- elem_disc_u:set_source("ISourceU")


-----------------------------------------
-- D) Anfangswerte
-----------------------------------------
-- Initial values ("Anfangswerte")
function initial_value_u(x, y)
	if (x*x+y*y<4) then return uinfty+0.6
	else return uinfty end
end

function initial_value_w(x, y)
  return winfty
end


-----------------------------------------
-- E) Randwerte
-----------------------------------------
local dirichlet_bnd = DirichletBoundary()
dirichlet_bnd:add(uinfty, "u", "Boundary")
dirichlet_bnd:add(winfty, "w", "Boundary")


-----------------------------------------
-- F) Diskretisierung auf ganzem Gebiet
-----------------------------------------
local domain_disc = DomainDiscretization(approx_space)
domain_disc:add(elem_disc_u)
domain_disc:add(elem_disc_w)
domain_disc:add(dirichlet_bnd)

-----------------------------------------
-- F) Aufsetzen des Lösers
--    (using 'util/solver_util.lua')
-----------------------------------------
---
---
local newton_conv_check = ConvCheck()
newton_conv_check:set_maximum_steps(2)
newton_conv_check:set_minimum_defect(1e-12)
newton_conv_check:set_reduction(1e-6)
newton_conv_check:set_verbose(true)
newton_conv_check:set_supress_unsuccessful(true)

local linear_conv_check = ConvCheck()
linear_conv_check:set_maximum_steps(10)
linear_conv_check:set_minimum_defect(1e-12)
linear_conv_check:set_reduction(1e-6)
linear_conv_check:set_verbose(true)
linear_conv_check:set_supress_unsuccessful(true)


local smoother = SymmetricGaussSeidel ()
smoother:enable_consistent_interfaces(false)
smoother:enable_overlap(false)

local base_solver = LU()

local gmg = GeometricMultiGrid(approx_space)
gmg:set_base_solver(base_solver)
gmg:set_smoother(smoother)
gmg:set_base_level(0)
gmg:set_cycle_type("V")
gmg:set_num_presmooth(3)
gmg:set_num_postsmooth(3)
gmg:set_rap(false)
gmg:set_smooth_on_surface_rim(false)
gmg:set_emulate_full_refined_grid(false)
gmg:set_gathered_base_solver_if_ambiguous(false)


local linear_solver = BiCGStab()
linear_solver:set_convergence_check(linear_conv_check)
linear_solver:set_preconditioner(gmg)


local nonlinear_solver = NewtonSolver()
nonlinear_solver:set_linear_solver(linear_solver)
nonlinear_solver:set_convergence_check(newton_conv_check)


--[[ local solver_desc = {

  -- BiCGStab 
	linSolver = {

	 precond = {
		  type		= "gmg",
		  approxSpace	= approx_space,
		  smoother	= "sgs",
		  baseSolver	= "lu"
	 },
	
  },
	
}--]]


-----------------------------------------
-- G) solve
-----------------------------------------

-- local nonlinear_solver = util.solver.CreateSolver(solver_desc)
print (nonlinear_solver)
-- print(nonlinear_solver:config_string())

print("<script_solving>")
--local A = AssembledLinearOperator(domainDisc)
local u_start = GridFunction(approx_space)
--local b = GridFunction(approx_space)
--u_start:set(0.0)
-- domainDisc:adjust_solution(u)
-- domainDisc:assemble_linear(A, b)

Interpolate("initial_value_u", u_start, "u")
Interpolate("initial_value_w", u_start, "w")

vtk = VTKOutput()

local start_time = 0
local end_time = tDiff
local time_steps = 200

if (script_method == "Sequential") then
	print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print("    |                               Nonlinear Solve                                    |")
	print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print("<script_sequential>")
	time = BraidTimer()
	time:start()



	local dt = (end_time-start_time)/time_steps
	util.SolveNonlinearTimeProblem(u_start, domain_disc, nonlinear_solver, vtk,
			"CardioMonoDomain2D", "ImplEuler", 1, start_time, end_time, dt);


	integration_time = time:get()
	print("</script_sequential>")
elseif (script_method == "XBraid") then
	print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print("    |                                     XBraid                                       |")
	print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("<script_xbraid>")
print("<xbraid_preparation>")
pit_desc.norm ="l2"
pit_desc.use_residual = false
pit_desc.cfactor = "4_50"--"5_5_4_2"
pit_desc.level_num_ref = "4_4" --  "5_5_5_5"
pit_desc.hierarchy ={
	max_levels = 2,
}
pit_desc.spatial_coarsen_and_refine = false
pit_desc.store_values = 0
pit_desc.access_level = 0
pit_desc.print_level = 3
pit_desc.print_file = "xbraid_output"
pit_desc.skip_downcycle_work = true
pit_desc.cycle = "V"
pit_desc.cf_weight = 1
pit_inst.spatial_norm = util.xbraid.create_norm(pit_desc, pit_inst)
print("")
util.xbraid.set_domain_disc(pit_desc, pit_inst,domain_disc)
print("")

--ø  pit_desc.logger = "jobb"
pit_desc.time_interval = {
	start_time = start_time,
	end_time = end_time,
	time_steps = time_steps
}
pit_desc.conv_check = {
	max_iterations = 30,
	absolut_tol = 1e-5,
	relative_tol = 0
}
pit_desc.temporal_norm = 3

pit_desc.driver = {
	name = "Nonlinear",
}

pit_desc.integrator = {
	name = "ThetaIntegratorNL",
	nonlinear_solver = nonlinear_solver,
	theta = 1,
}


util.xbraid.create_logger(pit_desc,pit_inst)
util.xbraid.create_integrator(pit_desc, pit_inst)
print("")
util.xbraid.create_driver(pit_desc, pit_inst)
print("")
print(pit_inst)
xbraid = util.xbraid.create_instance(pit_desc, pit_inst)
print("")
print("</xbraid_preparation>")
print("<xbraid_apply>")
	--[[
    	print("<test>")
        print("<test_integrator_nl>")
        integrator = pit_inst.driver
        u_stop = u_start:clone()
        print("t_end="..pit_desc.time_interval.end_time)
        print("t_start"..pit_desc.time_interval.start_time)
        local dt = (end_time-start_time)/time_steps
        integrator:apply(u_start:clone(), pit_desc.time_interval.start_time + dt, u_start, pit_desc.time_interval.start_time)
        print("</test_integrator_nl>")
        print("<test_integrator_nl>")
        -- [ [
        integrator = pit_inst.integrator.integrator
        u_stop = u_start:clone()
        print("t_end="..pit_desc.time_interval.end_time)
        print("t_start"..pit_desc.time_interval.start_time)
            local dt = (end_time-start_time)/time_steps
        integrator:apply(u_start:clone(), pit_desc.time_interval.start_time + dt, u_start, pit_desc.time_interval.start_time)

        ] ]
print("</test_integrator_nl>")
print("</test>")
	--]]
util.xbraid.apply(pit_desc, pit_inst, u_start)
print("</xbraid_apply>")
print("</script_xbraid>")

end




walltime:stop() -- stop timer for entire script execution
if space_time_communicator:get_temporal_rank() == 0 then
	print("<integration time=\""..seconds.."\" unit=\"s\" />")
	total_time = walltime:get()
	t_rank = space_time_communicator:get_temporal_rank()
	print("<summary walltime=\"".. total_time .."\" rank=\""..t_rank.."\" />")
end
print("</script_solving>")


