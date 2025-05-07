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
space_time_communicator = util.xbraid.split_communicator()

ug_load_script("util/refinement_util.lua")

script_method = util.GetParam("--method", "XBraid", "target method of this script")

XBraidParameters = {
	time = { t_0 = util.GetParamNumber("--start-time", 0.0, "start time"),
			 t_end = util.GetParamNumber("--end-time", 1e+5, "end time"), -- input will be ignored if problem has endtime
			 numtime = util.GetParamNumber("--num-time", 256, "number of time points")},

	driver = util.GetParam("--driver", "Integrator", ""),
	integrator = util.GetParam("--integrator", "ThetaTimeStep", ""),
	redirect = util.GetParamNumber("--redirect", 0, "") == 1,
	theta = util.GetParamNumber("--theta", 1, ""),
	order = util.GetParamNumber("--order", 2, ""),
	grid_step = util.GetParamNumber("--gridstep", 1, ""),
	coarse = util.GetParamNumber("--coarse", 3, ""),

	cfactor_parameter = util.GetParam("--cfactor", "16_4_4_4", ""),
	cfactor_default = util.GetParamNumber("--cfactor-default", 2, ""),
	max_level = util.GetParamNumber("--maxlevel", 20, " "),
	mgrit_cycle_type = util.GetParam("--cycle", "V", ""),
	mgrit_relax_type = util.GetParam("--relax", "FCF", ""),
	store_values = util.GetParamNumber("--store-level", 0, ""),
	print_level = util.GetParamNumber("--printlevel", 1, ""),
	access_level = util.GetParamNumber("--accesslevel", 1, ""),
	sequential = util.GetParam("--sequential", 0, "") == 1,
	initializer = util.GetParam("--initializer", "StartValue", "") == 1,
	temporal_norm = util.GetParamNumber("--tol-norm-braid", 3, ""),
	conv_check = {
		max_iter = util.GetParamNumber("--maxiter", 50, " "),
		reduction = util.GetParamNumber("--tol-rel-braid", 0, ""), -- todo if 0 then not set
		absolute = util.GetParamNumber("--tol-abs-braid", 0, ""), -- todo if 0 then not set
	},
	skip_downcycle_work = util.GetParamNumber("--boolskipdown", 1, "") == 1,
	max_refinement = 10,
	spatial_coarsen_and_refine = util.GetParamNumber("--spatial-refine", 0, "") == 1,
	min_coarsening = 2,

	outputfile = "integrator_out",
	use_residual = util.GetParamNumber("--use-residual", 0, " 0 xbraid residual, 1 use residual") == 1,

	time_refinement = util.GetParamNumber("--time-refine", 0, "") == 1,
	max_time_refinement = util.GetParamNumber("--max-time-refine", 0, ""), -- todo
	rich_refine = util.GetParamNumber("--rich-refine", 2, ""),
	rich_bound = util.GetParamNumber("--rich-bound", 1.1, ""),
	richardson_estimation = util.GetParamNumber("--rich-est", 0, "") == 1,
	richardson_extrapolation = util.GetParamNumber("--rich-ext", 0, "") == 1,
	richardson_local_order = util.GetParamNumber("--rich-order", 2, ""),
	verbose = util.GetParamNumber("--verbose", 1, "") == 1,
	norm = util.GetParam("--norm", "displacement", ""),
	level_config = util.GetParam("--level-config", "simple", ""),
	level_numref = util.GetParam("--level-num-ref", "simple", "")
}
-- complete some settings basing on the given parameter
XBraidParameters.printfile = "000 " .. XBraidParameters.driver .. "_" .. XBraidParameters.time.numtime .. "_" .. XBraidParameters.max_level .. "_" .. XBraidParameters.mgrit_cycle_type .. "_" .. XBraidParameters.mgrit_relax_type .. ".mgrit"
XBraidParameters.sync = XBraidParameters.time_refine and XBraidParameters.rich_est
XBraidParameters.cfactor = util.xbraid.cfactor_from_string(XBraidParameters.cfactor_parameter)
XBraidParameters.level_num_ref = util.xbraid.numref_from_string(XBraidParameters.level_numref)
print(XBraidParameters.cfactor)


if XBraidParameters.redirect then
	repl = ReplaceStandardStream()
	repl:set_space_time_comm(space_time_communicator)
	repl:apply()
end



-- Parse parameters and print help
local gridName	= util.GetParam("--grid", "grids/monodomain_laplace_sample_grid_2d.ugx","filename of underlying grid")
local numRefs		= util.GetParamNumber("--numRefs", 4, "number of refinements")

util.CheckAndPrintHelp("Beispiel: Mono-Domain-Gleichungen (2D)");

-- ===============================================================================================
-- Create geometric grid and refine
-- ===============================================================================================

local blockSize = 2 -- select 2 for 2x2 point-block system
InitUG(2, AlgebraType("CPU", blockSize));  


-- Load a domain without initial refinements.
local mandatorySubsets = {"Inner", "Boundary"}
dom = util.CreateDomain(gridName, 0, mandatorySubsets)

-- Refine the domain (redistribution is handled internally for parallel runs)
print("refining...")
util.refinement.CreateRegularHierarchy(dom, numRefs, true)

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

local Vr = -85   -- mV
local Vs = -75   -- mV
local Vp =  15   -- mV

local scaleV = (Vp-Vr)        -- Skalierung [mV]
local a = (Vs-Vr)/(Vp-Vr)     -- Anregeschwellwert [1]


local VoltPerMilliVolt = 1e-3
local Cm  = 1.0 * VoltPerMilliVolt --  uF/cm^2,  where: 1 F=C/V = As / V => 1 mF = ms A/V  => 1 uF = ms * mA /V

local c1 = 0.175 -- ms^{-1}
local c2 = 0.03  -- ms^{-1}
local c3 = 0.011 -- ms^{-1}
local  b = 0.55  -- [1]


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

local approxSpace = ApproximationSpace(dom)

approxSpace:add_fct("u", "Lagrange", 1)   -- lineare Ansatzfunktionen
approxSpace:add_fct("w", "Lagrange", 1)

print("approximation space:")
approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()


-----------------------------------------
-- C) Elementdiskretisierung (FV)
-----------------------------------------
local elemDisc = {}
elemDisc["u"] = ConvectionDiffusion("u", "Inner", "fv1")  -- potential
elemDisc["w"] = ConvectionDiffusion("w", "Inner", "fv1")  -- gating variable


-- elemDisc["u"]:set_mass_scale(1.0)  -- default
-- elemDisc["w"]:set_mass_scale(1.0)  -- default

-----------------------------------------
-- C).1 diffusiontensor
-- Df = sigma/(Cm*chi)  
-----------------------------------------
local diffTensorU = ConstUserMatrix(DDiff/4.0)
diffTensorU:set_entry(0,0, DDiff)


elemDisc["u"]:set_diffusion(diffTensorU)
elemDisc["w"]:set_diffusion(0.0)

-----------------------------------------
-- C).2 Reaktionsterme:
-- \dot u = (c1*u*(u-a)*(1.0-u) - c2*w)
-- \dot w = b*(v - d*w)
-----------------------------------------


function ReactionU(u,w)
	return (c1*u*(u-a)*(1.0-u) - c2*w) * (-1.0)
end
function ReactionW(u,w)
	return c3*(u - b*w)*(-1.0)
end
function ReactionU_dU(u,w)
	return (a*(2.0*u-1.0) + (2.0-3.0*u)*u ) * (-c1)
end

function ReactionU_dW(u,w)
	return -c2 * (-1.0) --  todo simplify to c2
end

function ReactionW_dU(v,w)
	return c3*(-1.0) --  todo simplify to -1.0*c3
end
function ReactionW_dW(u,w)
	return -c3*b*(-1.0)--  todo simplify c3*b
end



local nonlinearGrowth = {}
nonlinearGrowth["u"] = LuaUserFunctionNumber("ReactionU", 2)
nonlinearGrowth["u"]:set_input(0, elemDisc["u"]:value())
nonlinearGrowth["u"]:set_input(1, elemDisc["w"]:value())
nonlinearGrowth["u"]:set_deriv(0, "ReactionU_dU")
nonlinearGrowth["u"]:set_deriv(1, "ReactionU_dW")


nonlinearGrowth["w"] = LuaUserFunctionNumber("ReactionW", 2)
nonlinearGrowth["w"]:set_input(0, elemDisc["u"]:value())
nonlinearGrowth["w"]:set_input(1, elemDisc["w"]:value())
nonlinearGrowth["w"]:set_deriv(0, "ReactionW_dU")
nonlinearGrowth["w"]:set_deriv(1, "ReactionW_dW")


elemDisc["u"]:set_reaction(nonlinearGrowth["u"])          
elemDisc["w"]:set_reaction(nonlinearGrowth["w"])  

-----------------------------------------
--  C).3 Quellen = Externe Ströme (nicht benutzt!)
-- 1.0 / (Cm*scaleV)
-----------------------------------------
-- function ISourceU(x,y,t,si)  return 0.0/(Cm*scaleV) end
-- elemDisc["u"]:set_source("ISourceU") 


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
local dirichletBND = DirichletBoundary()
dirichletBND:add(uinfty, "u", "Boundary")
dirichletBND:add(winfty, "w", "Boundary")


-----------------------------------------
-- F) Diskretisierung auf ganzem Gebiet
-----------------------------------------
local domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(elemDisc["u"])
domainDisc:add(elemDisc["w"])
domainDisc:add(dirichletBND)

-----------------------------------------
-- F) Aufsetzen des Lösers
--    (using 'util/solver_util.lua')
-----------------------------------------
---
---
local newton_conv_check = ConvCheck()
newton_conv_check:set_maximum_steps(100)
newton_conv_check:set_minimum_defect(1e-12)
newton_conv_check:set_reduction(1e-6)
newton_conv_check:set_verbose(true)



local solverDesc = {
	
	-- Newton's method for non-linear problem
	type = "newton",
	convCheck = newton_conv_check,
  -- BiCGStab 
	linSolver = {
	 type = "bicgstab",
	
	 precond = {
		  type		= "gmg",
		  approxSpace	= approxSpace,
		  smoother	= "sgs",
		  baseSolver	= "lu"
	 },
	
  },
	
}

-----------------------------------------
-- G) Lösen
-----------------------------------------

local nonlinear_solver = util.solver.CreateSolver(solverDesc)
print (nonlinear_solver)
print(nonlinear_solver:config_string())

print("\nsolving...")
local A = AssembledLinearOperator(domainDisc)
local u_start = GridFunction(approxSpace)
local b = GridFunction(approxSpace)
u_start:set(0.0)
-- domainDisc:adjust_solution(u)
-- domainDisc:assemble_linear(A, b)

Interpolate("MyInitialValueU", u_start, "u")
Interpolate("MyInitialValueW", u_start, "w")


log_job = Paralog()
log_job:set_comm(space_time_communicator)
log_job:set_filename("job")
log_job:init()


--local vtk=VTKOutput()
--vtk:select_nodal("t", "TNodal")
--out = BraidHeatCheck()
--out:set_log(log_job)
--out:set_solution_name(vtk,"solution")
--out:set_diff_name(vtk,"error")
--out:set_vtk_write_mode(false, false)

if (script_method == "SequentialNonLinearSolver") then
	time = BraidTimer()
	time:start()

	local startTime = 0
	local endTime = tDiff
	local dt = (endTime-startTime)/200.0
	util.SolveNonlinearTimeProblem(u_start, domainDisc, nonlinear_solver, VTKOutput(),
			"CardioMonoDomain2D", "ImplEuler", 1, startTime, endTime, dt);

	print("done")

	integration_time = time:get()
elseif (script_method == "XBraid") then
	print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
	print("    |                                     XBraid                                       |")
	print("    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")

	print("[info] - Created logger for job-output")
	if XBraidParameters.level_config == "simple" then -- all level use same methods
		XBraidParameters.default = {}
		XBraidParameters.default.integrator = "ThetaIntegratorNL";
		XBraidParameters.default.solver = nonlinear_solver;
		XBraidParameters.default.conv_check = newton_conv_check;

	elseif XBraidParameters.level_config == "finecoarse" then -- all level use same methods

		print("all level use same method except fine level")
		XBraidParameters.fine = {}
		XBraidParameters.fine.integrator = "ThetaTimeStep";
		XBraidParameters.fine.solver = nonlinear_solver;
		XBraidParameters.fine.conv_check = conv_check;
		XBraidParameters.fine.loose_tolerance = 1e-10;
		XBraidParameters.fine.tight_tolerance = 1e-20;
		XBraidParameters.fine.target_tolerance = XBraidParameters.conv_check.absolute;

		XBraidParameters.coarse = {}
		XBraidParameters.coarse.integrator = "ThetaTimeStep";
		XBraidParameters.coarse.solver = nonlinear_solver;
		XBraidParameters.coarse.conv_check = conv_check_coarse;

	end
	-- override norm function to use specific problem metrics
	function util.xbraid.set_norm(desc, inst)
		print("heat - set_norm ! not xbraid")
		if desc.norm == "l2" then
			print("l2" )
			norm = BraidEuclidianNorm()
		else
		end
		inst.braid:set_norm_provider(norm)
		print()
	end



	inst = inst or {}



	print("heat - xbraid.set_domain_disc")
	util.xbraid.set_domain_disc(XBraidParameters,inst, domainDisc)
	util.xbraid.set_initializer(XBraidParameters,inst, initializer)

	print("heat - xbraid.set_communicator")
	util.xbraid.set_communicator(XBraidParameters,inst, space_time_communicator)

	print("heat - xbraid.create_driver")
	util.xbraid.create_driver(XBraidParameters,inst)

	print("heat - xbraid.create_instance")
	xbraid = util.xbraid.create_instance(XBraidParameters,inst)

	-- braid:print_settings()
	time = BraidTimer()

	time:start()

	-- sgt = SpatialGridTransfer()
	-- sgt:set_approx_space(approxSpace)
	-- sgt:set_domain_disc(domainDisc)
	-- sgt:init()

	-- inst.driver:set_spatial_grid_transfer(sgt)
	-- u_start_nontop = sgt:make_nontop(u_start)
	--print("nontop")
	--vtk_output_observer = MockObserver()
	-- inst.driver:attach_observer(vtk_output_observer);
	-- xbraid:apply(u_start_nontop, parameters.end_time, u_start_nontop, parameters.start_time)
	xbraid:apply(u_start, tDiff, u_start, 0)
	time:stop()
	integration_time = time:get()
end





if XBraidParameters.redirect then
	repl:undo() -- give back the std::cout to terminal
end
walltime:stop() -- stop timer for entire script execution
if space_time_communicator:get_temporal_rank() == 0 then
	print("\n\n")
	print(integration_time .. " seconds for time integration")
	print(walltime:get() .. " seconds wall time  (rank=" .. space_time_communicator:get_temporal_rank() .. ")")
end
