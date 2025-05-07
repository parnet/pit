------------------------------------------------------------------------------
--
--   Lua - Script for poroelacticity
--
--   Author: Arne Naegel
--          (derived from on solid_mechanics app by Raphael Prohl)
--
------------------------------------------------------------------------------
PluginRequired("ConvectionDiffusion")
PluginRequired("SmallStrainMechanics")
PluginRequired("Limex")
PluginRequired("SuperLU6")
PluginRequired("XBraidForUG4") -- xbraidforug4
ug_load_script("ug_util.lua")
local walltime = BraidTimer()
print()
print()
print()

-- Expand path.
local myPath = ug_get_current_path()
package.path = package.path..";".. myPath.."/../config/?.lua;".. myPath.."/?.lua"

pit_desc = {} -- xbraidforug4
pit_inst = {} -- xbraidforug4

space_time_communicator = util.xbraid.split_communicator(pit_desc,pit_inst)
script_method = util.GetParam("--method", "XBraid", "target method of this script")
print("<script method='"..script_method.."'/>")
ug_load_script("util/load_balancing_util_2.lua")
ug_load_script("util/profiler_util.lua")
-- ug_load_script("plugins/Limex/limex_util.lua")
-- ug_load_script("generic.lua")
-- ug_load_script("cryer.lua")
-- ug_load_script("footing.lua")
-- ug_load_script("../config/barry_mercer.lua")
-- util.biot.CheckAssertions()
-- ug_load_script("mandel.lua")


-- TIMES AND TIME-STEPPING
local startTime  = util.GetParamNumber("--start", 0.0, "end time")
local endTime    = util.GetParamNumber("--end", 1e+5, "end time")
local dtFrac     = util.GetParamNumber("--dtFrac", 1.0/256.0, "time step size")
local dtMinFrac  = util.GetParamNumber("--dtminFrac", 1e-2, "minimal admissible time step size")
local dtMaxFrac  = util.GetParamNumber("--dtmaxFrac", 0.1, "minimal admissible time step size (as fraction of tend)")
local dtRed      = util.GetParamNumber("--dtred", 0.5, "time step size reduction factor on divergence")


-- REFINEMENT
local numPreRefs   = util.GetParamNumber("--numPreRefs", 0, "number of pre-Refinements (before distributing grid)")
local numRefs      = util.GetParamNumber("--num-refs", 3, "total number of refinements (incl. pre-Refinements)") --4 --
local time_steps      = util.GetParamNumber("--time-steps", 256, "number of timesteps") --4 --


-- MORE ARGUMENTS
local ARGS = {
    discConfig 			= util.GetParam("--disc-id", "P1_stab"), -- cryer3d‚
    problemID 			= util.GetParam("--problem-id", "bm2D_P1"), -- cryer3d‚
    solverID 				= util.GetParam("--solver-id", "GMGKrylov"),  --  "FixedStressEX", "UzawaMG", "UzawaSmoother","UzawaMGKrylov"

    -- Global options.
    useVTK 					=  util.HasParamOption("--with-vtk", "Plot VTK"),
    useDebugIter 		=  util.HasParamOption("--with-debug-iter", "Activate debug solver."),
    -- doCheck 			=  util.HasParamOption("--with-check", ""),

    bSteadyStateMechanics = not util.HasParamOption("--with-transient-mechanics"), -- OPTIONAL: transient mechanics

    -- Grid.
    adaptiveRefs   	= util.GetParamNumber("--num-additional-adaptive-refs", 0, "total number of refinements (incl. pre-Refinements)"),

    -- Multigrid
    MGCycleType 		= util.GetParam("--mg-cycle-type", "F", "V,F,W"),
    MGBaseLevel 		= util.GetParamNumber("--mg-base-level", 0, "some non-negative integer"),
    MGNumSmooth 		= util.GetParamNumber("--mg-num-smooth", 2, "some positive integer"),
    MGSmootherType 	= util.GetParam("--mg-smoother-type", "uzawa3", "uzawa3,uzawa,cgs"),
    MGDebugLevel 		= util.GetParam("--mg-debug-level", 0, "some non-negative integer"),


    -- Time stepping.

    -- LIMEX
    LimexTOL     		= util.GetParamNumber("--limex-tol", 1e-3, "TOL"),
    LimexNStages 		= util.GetParamNumber("--limex-num-stages", 4, "number of LIMEX stages q"),

}


print ("dtFrac        ="..dtFrac)

print ("MGSmootherType="..ARGS.MGSmootherType)
print ("MGNumSmooth   ="..ARGS.MGNumSmooth)
print ("MGCycleType   ="..ARGS.MGCycleType)
print ("MGBaseLevel   ="..ARGS.MGBaseLevel)
print ("MGDebugLevel  ="..ARGS.MGDebugLevel)

GetLogAssistant():set_debug_level("LIB_DISC_MULTIGRID", ARGS.MGDebugLevel);
--SetDebugLevel("LIB_DISC_MULTIGRID", 0)

-- Set parameters
local kperm   = 1e-0 -- m/s 1e-3
local poro    = 0.2
local nu      = 0.25

local EYoung  = 2.0 * 1e+2                  -- kPa 2.0 * 1e+4
local Kmedium = EYoung/(3.0*(1.0-2.0*nu))
local Kfluid  = 2.2 * 1e+6                  -- kPa -- 2.2 * 1e+6 --

print ("Kmedium = "..Kmedium)
print ("Kfluid  = "..Kfluid)


-- Setup for discretization.
local config = {
    ["P1_stab"] = BiotDiscConfig("ux,uy", 1, "p", 1, 1.0/12.0), -- P1-P1 + stabilization
    ["P2P1"]    = BiotDiscConfig("ux,uy", 2, "p", 1, 0.0),      -- P2-P1
}
--deleeuw2d--deleeuw2d -- cryer3d --cryer2d -- mandel3d --, mandel--, cryer3d
local problemList = {
    -- legacy style
    ["deleeuw2d"] = deleeuw2d,
    ["deleeuw3d"] = deleeuw3d,
    ["deleeuw3dTet"] = deleeuw3dTet,
    ["cryer3d"] = cryer3d,
    ["cryer3dTet"] = cryer3dTet,
    ["footing2D"] = footing2D,
    ["footing2D_tri"] = footing2D_tri,
    ["footing3D"] = footing3D,

    ["bm2D_tri"] = barrymercer2D_tri,

    -- C++ style
    ["bm2D_P1"] = BarryMercerProblem2dCPU1(config["P1_stab"]),
    ["bm2D_P2"] = BarryMercerProblem2dCPU1(config["P2P1"]),
}

local problem = problemList[ARGS.problemID]
if (not problem) then
    print ("ERROR: Problem '".. ARGS.problemID.. "' not found")
    quit()
end

-- Parse command line args
if (problem.parse_cmd_args) then problem:parse_cmd_args() end

--problem:init(kperm, poro, nu, 1.0/Kmedium, 1.0/Kfluid, 0.0)
if(problem.init) then
    problem:init(kperm, poro, nu, 1.0/Kmedium, 0.0, 0.0)
end

-- Characteristic time (implemented by C++ object)
local charTime = problem:get_char_time()
print("charTime="..charTime)

-- Start & end time
startTime = 0.0
if (problem.start_time) then startTime = problem:start_time() end

endTime   = 2.0*charTime
if (problem.end_time) then endTime = problem:end_time() end

-- Time stepping
local dt  = dtFrac*endTime
local dtMin = dtMinFrac
local dtMax = endTime
print ("dt="..dt)


if (problem == mandel) then
    local time = 1e-4
    while (time<=10.0) do
        problem:create_test_data(time*charTime)
        time = time *10.0;
    end
end


local doSteadyState = false
local doTransient = true

----------------------------------
----------------------------------
--  Settings
----------------------------------
----------------------------------

local dim = problem.dim or 2
local cpu = problem.cpu or 1    -- default: block dim+1

-- Order for Ansatz functions.
local porder = problem.porder or problem:get_porder() or 1
local uorder = problem.uorder or problem:get_uorder() or (porder+1)

InitUG(dim, AlgebraType("CPU", cpu));



-- OUTPUT-ASSISTANT FOR SEVERAL PROCESSES
GetLogAssistant():enable_file_output(true, "output_p_"..ProcRank()..
        "_Lev"..numRefs..".txt")
GetLogAssistant():set_debug_level("SchurDebug", 7);
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Domain / ApproximationSpace setup
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------


local balancerDesc = {
    hierarchy = {
        --type = "standard",
        -- maxRedistProcs = ARGS.redistProcs,

        -- minElemsPerProcPerLevel = ARGS.minElemsPerProcPerLevel,
        -- qualityRedistLevelOffset = ARGS.qualityRedistLevelOffset,
        -- intermediateRedistributions = ARGS.intermediateRedistributions,

        type            = "standard",
        minElemsPerProcPerLevel   = 32,
        maxRedistProcs        = 120,
        qualityRedistLevelOffset  = 2,
        intermediateRedistributions = true,

        {
            upperLvl = 0,
            maxRedistProcs = 40
        },

        {
            upperLvl = 2,
            maxRedistProcs = 120
        },
    }
} -- balancerDesc


-- Create, Load, Refine and Distribute Domain
local gridName = problem.gridName or problem:get_gridname()
-- local dom = problem:create_domain(numRefs, numPreRefs)

local mandatorySubsets = problem.mandatorySubsets or nil
local dom = util.CreateDomain(gridName, 0, mandatorySubsets)

-----------------------------------------------------------------
--  Special refinement
-----------------------------------------------------------------
--
local doHangingNode = false

if (ARGS.adaptiveRefs>0) then
    doHangingNode = true
    local refiner = HangingNodeDomainRefiner(dom);


    local radius = 0.1

    for i=1,numRefs do
        MarkForRefinement_All(refiner)
        refiner:refine()
    end

    local center = Vec()
    center:set_coord(0, 0.25)
    center:set_coord(1, 0.25)

    for i=1,ARGS.adaptiveRefs do
        print ("Ref")
        MarkForRefinement_FacesInSphere(dom, refiner, center, radius)
        -- MarkForAdaption_ElementsTouchingSubset(dom, refiner, dom:subset_handler(), 4, "refine")
        refiner:refine()
        radius = radius*0.7
    end
else
    -- Classic refinement
    util.refinement.CreateRegularHierarchy(dom, numRefs, true, balancerDesc)
end
-----------------------------------------------------------------
--  Approximation Space
-----------------------------------------------------------------
local approxSpace = util.biot.CreateApproxSpace(dom, dim, uorder, porder)

--------------------------------------------------------------------------------
-- Problem Setup
--------------------------------------------------------------------------------
print("FE discretization...")
local bSteadyStateMechanics = ARGS.bSteadyStateMechanics -- true

-- For computing consistent initial values.
local domainDisc0 = DomainDiscretization(approxSpace)
problem:add_elem_discs(domainDisc0, bSteadyStateMechanics)  -- implemented by C++ object
problem:add_boundary_conditions(domainDisc0, bSteadyStateMechanics)  -- implemented by C++ object

-- For time-dependent problem.
local domainDiscT = DomainDiscretization(approxSpace)
problem:add_elem_discs(domainDiscT, bSteadyStateMechanics)  -- implemented by C++ object
problem:add_boundary_conditions(domainDiscT, bSteadyStateMechanics)  -- implemented by C++ object


-- For Uzawa fixed-stress smoother.
local uzawaSchurUpdateDisc = DomainDiscretization(approxSpace)
problem:add_uzawa_discs(uzawaSchurUpdateDisc, bSteadyStateMechanics)  -- implemented by C++ object
print("done!")


if doHangingNode then
    local hangingNodeConstraints = SymP1Constraints()
    domainDisc0:add(hangingNodeConstraints)
    domainDiscT:add(hangingNodeConstraints)
    -- uzawaSchurUpdateDisc:add(hangingNodeConstraints)
    print("Added hanging nodes!")
end



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
local u = GridFunction(approxSpace)
local dbgVector=u:clone()
--------------------------------------------------------------------------------


----------------------------------------
-- create algebraic Preconditioner
----------------------------------------

local jac = Jacobi()
jac:set_damp(0.66)
local gs = GaussSeidel()
local sgs = SymmetricGaussSeidel()
local bgs = BackwardGaussSeidel()
bgs:enable_overlap(true)
gs:enable_overlap(true)
sgs:enable_overlap(true)
local ilu = ILU()
--ilu:set_beta(-0.5);
local ilut = ILUT()
ilut:set_threshold(1e-3)

--local egs_weights = u:clone();
--egs_weights:set(1.0);
--Interpolate(0.1, egs_weights, "p")

local egs = ElementGaussSeidel() -- patches per node
egs:select_schur_cmp({"p"}, 4.0)
egs:set_relax(0.125)

local cgs = ComponentGaussSeidel(1.0, {"p"}) -- patches per node
cgs:set_alpha(1.0)
cgs:set_beta(1.0) --- 0 > 0.25  (beta=0.0: no pressure change) -- 1.0: works
cgs:set_weights(true)


local ssc_vanka_space
if (dim == 2) then
    ssc_vanka_space = VertexCenteredVankaSubspace2dCPU1({"p"}, {"ux", "uy"})
else
    ssc_vanka_space = VertexCenteredVankaSubspace3dCPU1({"p"}, {"ux", "uy", "uz"})
end


local ssc = SequentialSubspaceCorrection(1.0)
ssc:set_vertex_subspace(ssc_vanka_space)

local dbgWriter = GridFunctionDebugWriter(approxSpace)
local uzawaSchurUpdateOp = AssembledLinearOperator()
uzawaSchurUpdateOp:set_discretization(uzawaSchurUpdateDisc)



--- Factory for Uzawa iteration.
-- @function createUzawaIteration
-- @param #string sSchurCmp  Schur complement will be built for this unknown.
-- @param aiForward Approximate Inverse (forward problem)
-- @param aiSchur Approximate Inverse (Schur complement)
-- @param aiBackward Approximate Inverse (backward problem)
function createUzawaIteration(sSchurCmp, aiForward, aiSchur, aiBackward, uzawaSchurUpdateOp, uzawaSchurWeight)

    local uzawa = UzawaBase(sSchurCmp)
    local weight = uzawaSchurWeight or 1.0
    if (aiForward) then uzawa:set_forward_iter(aiForward)  end
    if (aiSchur) then uzawa:set_schur_iter(aiSchur) end
    if (aiBackward) then uzawa:set_backward_iter(aiBackward)  end

    uzawa:set_schur_operator_update(uzawaSchurUpdateOp, weight)
    -- uzawa:set_debug(dbgWriter)
    return uzawa
end

local uzawaForward = {}
local uzawaBackward = {}

local uzawaWeight = 1.0
local uzawaForward1 = createUzawaIteration("p", gs, Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
local uzawaBackward1 = createUzawaIteration("p", nil, Jacobi(0.66), bgs, uzawaSchurUpdateOp, uzawaWeight)

local uzawaForward2 = createUzawaIteration("p", SymmetricGaussSeidel(), Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
local uzawaBackward2 = createUzawaIteration("p", nil, Jacobi(0.66), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)

uzawaForward[3] = createUzawaIteration("p", SymmetricGaussSeidel(), SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
uzawaBackward[3] = createUzawaIteration("p", nil, SymmetricGaussSeidel(), SymmetricGaussSeidel(), uzawaSchurUpdateOp, uzawaWeight)

uzawaForward[0] = createUzawaIteration("p", Jacobi(0.5), Jacobi(0.5), nil, uzawaSchurUpdateOp, uzawaWeight)
uzawaBackward[0] = createUzawaIteration("p", nil, Jacobi(0.5), Jacobi(0.5), uzawaSchurUpdateOp, uzawaWeight)

uzawaForward[4] = createUzawaIteration("p", CG(Jacobi()), CG(Jacobi()), nil, uzawaSchurUpdateOp, uzawaWeight)
uzawaBackward[4] = createUzawaIteration("p", nil, CG(Jacobi()), CG(Jacobi()), uzawaSchurUpdateOp, uzawaWeight)
--[[
local pi=LinearIteratorProduct()
pi:add_iterator(SymmetricGaussSeidel())
pi:add_iterator(SymmetricGaussSeidel())

local uzawaForward2 = createUzawaIteration("p", pi, SymmetricGaussSeidel(), nil, uzawaSchurUpdateOp, uzawaWeight)
local uzawaBackward2 = createUzawaIteration("p", nil, SymmetricGaussSeidel(), pi, uzawaSchurUpdateOp, uzawaWeight)
--]]

--local uzawaForward = createUzawaIteration("p", Jacobi(0.66), Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
--local uzawaBackward = createUzawaIteration("p", nil, Jacobi(0.66), Jacobi(0.66), uzawaSchurUpdateOp, uzawaWeight)
local uzawaSym = createUzawaIteration("p", gs, sgs, bgs, uzawaSchurUpdateOp, uzawaWeight)
--local uzawaBackward = createUzawaIteration("p", nil, Jacobi(0.5), Jacobi(0.66), uzawaSchurOp, uzawaWeight)
local uzawa = uzawaForward



local preSmoother
local postSmoother

if (ARGS.MGSmootherType == "uzawa") then
    preSmoother  = uzawaForward1
    postSmoother = uzawaBackward1
elseif (ARGS.MGSmootherType == "uzawa0") then
    preSmoother  = uzawaForward[0]
    postSmoother = uzawaBackward[0]
elseif (ARGS.MGSmootherType == "uzawa2") then
    preSmoother  = uzawaForward2
    postSmoother = uzawaBackward2
elseif (ARGS.MGSmootherType == "uzawa3") then
    preSmoother  = uzawaForward[3]
    postSmoother = uzawaBackward[3]
elseif (ARGS.MGSmootherType == "uzawa4") then
    preSmoother  = uzawaForward[4]
    postSmoother = uzawaBackward[4]
elseif (ARGS.MGSmootherType == "cgs") then
    preSmoother  = cgs
    postSmoother = cgs
elseif (ARGS.MGSmootherType == "vanka-ssc") then
    preSmoother  = ssc
    postSmoother = ssc
elseif (ARGS.MGSmootherType == "sgs") then
    preSmoother  = sgs
    postSmoother = sgs
else
    quit()
end

-------------------------
-- create GMG
-- TODO: Move to
-------------------------

-- Base Solver
local	baseConvCheck = ConvCheck()
baseConvCheck:set_maximum_steps(5000)
baseConvCheck:set_reduction(1e-12)
baseConvCheck:set_verbose(false)

local	base = BiCGStab()
base:set_preconditioner(jac)
base:set_convergence_check(baseConvCheck)

local	baseCG = CG()
baseCG:set_preconditioner(jac)
baseCG:set_convergence_check(baseConvCheck)

-- exact base solver
local	baseLU = LU()
local	superLU = SuperLU()


-- Geometric Multi Grid
local	gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDiscT)
gmg:set_base_level(ARGS.MGBaseLevel)
gmg:set_base_solver(superLU)
gmg:set_presmoother(preSmoother)
gmg:set_postsmoother(postSmoother)
gmg:set_cycle_type(ARGS.MGCycleType) -- 1:V, 2:W -- "F"
gmg:set_num_presmooth(ARGS.MGNumSmooth)
gmg:set_num_postsmooth(ARGS.MGNumSmooth)
gmg:set_rap(problem.bRAP or true)  -- mandatory, if set_stationary
gmg:set_smooth_on_surface_rim(true)

if (problem.bRAP) then print ("gmg:bRAP=true")
else print ("gmg:bRAP=false") end

-- gmg:set_debug(dbgWriter)


local transfer = StdTransfer()
transfer:enable_p1_lagrange_optimization(false)
--transfer:set_debug(dbgWriter)
gmg:set_transfer(transfer)

local gmgP = GeometricMultiGrid(approxSpace)
-- gmgP:set_discretization(domainDiscP)
gmgP:set_base_level(numPreRefs)  -- was 1 in Cincyj
gmgP:set_base_solver(baseLU)  -- was baseLU in Cincy
gmgP:set_presmoother(sgs)
gmgP:set_postsmoother(sgs)
gmgP:set_cycle_type("V") -- 1:V, 2:W -- "F"
gmgP:set_num_presmooth(3)
gmgP:set_num_postsmooth(3)
gmgP:set_rap(true)  -- mandatory, if set_stationary


local gmgU = GeometricMultiGrid(approxSpace)
-- gmgU:set_discretization(domainDiscU)
gmgU:set_base_level(numPreRefs)  -- was 1 in Cincyj
gmgU:set_base_solver(baseLU)  -- was baseLU in Cincy
gmgU:set_presmoother(sgs)
gmgU:set_postsmoother(sgs)
gmgU:set_cycle_type("V") -- 1:V, 2:W -- "F"
gmgU:set_num_presmooth(3)
gmgU:set_num_postsmooth(3)
gmgU:set_rap(true)  -- mandatory, if add_ionary
gmgU:set_debug(dbgWriter)

local uzawaTotal       = createUzawaIteration("p", ILUT(1e-8), ILUT(1e-8), nil, uzawaSchurUpdateOp, 1.0)      -- ???
local fixedStressLU    = createUzawaIteration("p", nil, ILUT(1e-12), ILUT(1e-12), uzawaSchurUpdateOp, 1.0)
local fixedStressSuperLU    = createUzawaIteration("p", nil, SuperLU(), SuperLU(), uzawaSchurUpdateOp, 1.0)
-- local fixedStressLU    = createUzawaIteration("p", nil, SuperLU(), SuperLU(), uzawaSchurUpdateOp, 1.0)
local fixedStressMG    = createUzawaIteration("p", nil, gmgU, gmgP, uzawaSchurUpdateOp, 1.0)


--------------------------------
-- debug solver /iter
--------------------------------

local p0 = 1.0
if (problem.modelParameter and problem.modelParameter.p0) then
    p0 = problem.modelParameter.p0
end


local cmpConvCheck = CompositeConvCheck(approxSpace)
cmpConvCheck:set_component_check("ux", p0*1e-14, 1e-6)
cmpConvCheck:set_component_check("uy", p0*1e-14, 1e-6)
if (dim==3) then
    cmpConvCheck:set_component_check("uz", p0*1e-14, 1e-6)
end
cmpConvCheck:set_component_check("p", p0*1e-14, 1e-6)
cmpConvCheck:set_maximum_steps(10)
cmpConvCheck:set_supress_unsuccessful(true)

local cmpConvCheck2 = CompositeConvCheck(approxSpace)
cmpConvCheck2:set_component_check("ux", p0*1e-12, 1e-6)
cmpConvCheck2:set_component_check("uy", p0*1e-12, 1e-6)
if (dim==3) then
    cmpConvCheck2:set_component_check("uz", p0*1e-12, 1e-6)
end
cmpConvCheck2:set_component_check("p", p0*1e-12, 1e-6)
cmpConvCheck2:set_maximum_steps(50)

cmpConvCheck2 = ConvCheck(200, 1e-25, 1e-20)

local dbgSolver = LinearSolver()
dbgSolver:set_preconditioner(gmg) -- cgs, gmg, uzawa
dbgSolver:set_convergence_check(cmpConvCheck)

local dbgIter= DebugIterator()
dbgIter:set_preconditioner(gmg)  -- gmg is the 'real' preconditioner
dbgIter:set_solver(dbgSolver)
dbgIter:set_solution(dbgVector)
dbgIter:set_random_bounds(-5e-6, 5e-6)
dbgIter:set_debug(dbgWriter)  -- print t_0 anf t_N

--------------------------------
-- create and choose a Solver
--------------------------------

local solver = {}

local convCheck = ConvCheck()
convCheck:set_maximum_steps(50)
convCheck:set_reduction(1e-10)
convCheck:set_minimum_defect(1e-14)
-- convCheck = cmpConvCheck  -- for DEBUGGING purposes

local iluSolver = LinearSolver()
iluSolver:set_preconditioner(ilut)
iluSolver:set_convergence_check(convCheck)

local jacSolver = LinearSolver()
jacSolver:set_preconditioner(jac)
jacSolver:set_convergence_check(convCheck)

local myIter = gmg
if (ARGS.useDebugIter) then myIter =  dbgIter end


solver["UzawaSmoother"] = LinearSolver()
solver["UzawaSmoother"]:set_preconditioner(uzawaForward2)
solver["UzawaSmoother"]:set_convergence_check(convCheck)

solver["GMG"] = LinearSolver()
solver["GMG"]:set_preconditioner(myIter) -- gmg, dbgIter
solver["GMG"]:set_convergence_check(convCheck) -- cmpConvCheck

solver["GMGKrylov"] = BiCGStab()
solver["GMGKrylov"]:set_preconditioner(myIter) -- gmg, dbgIter
solver["GMGKrylov"]:set_convergence_check(convCheck) -- cmpConvCheck

solver["GMRES"] = GMRES(3)
solver["GMRES"]:set_preconditioner(myIter) -- gmg, dbgIter
solver["GMRES"]:set_convergence_check(convCheck) -- cmpConvCheck

solver["FixedStressEX"] = LinearSolver()
solver["FixedStressEX"]:set_preconditioner(fixedStressSuperLU)
solver["FixedStressEX"]:set_convergence_check(convCheck)

solver["FixedStressEXKrylov"] = BiCGStab()
solver["FixedStressEXKrylov"]:set_preconditioner(fixedStressSuperLU)
solver["FixedStressEXKrylov"]:set_convergence_check(convCheck)

solver["FixedStressMG"] = LinearSolver() -- BiCGStab()
solver["FixedStressMG"]:set_preconditioner(fixedStressMG)
solver["FixedStressMG"]:set_convergence_check(convCheck)

solver["SuperLU"] = SuperLU() -- SuperLU

solver["LU"] = LinearSolver()
solver["LU"]:set_preconditioner(LU())
solver["LU"]:set_convergence_check(convCheck)

local cgSolver = CG()
cgSolver:set_preconditioner(myIter) --(gmg)
cgSolver:set_convergence_check(convCheck)


local sluSolver = SuperLU()

-- local luSolver = LinearSolver()
-- luSolver:set_preconditioner(LU())
-- luSolver:set_convergence_check(convCheck)

-- Select solver.
local lsolver = solver[ARGS.solverID]
--solver = jacSolver
--lsolver = iluSolver
--lsolver = gmgSolver
--lsolver = cgSolver
-- lsolver = bicgstabSolver
--lsolver = gmresSolver
--lsolver:set_compute_fresh_defect_when_finished(true)
-- lsolver = sluSolver

local vtk=VTKOutput()
vtk:select_nodal("p", "PNodal")
if (dim == 2) then vtk:select({"ux", "uy"}, "uNodal") end
if (dim == 3) then vtk:select({"ux", "uy", "uz"}, "uNodal") end

--vtk:select_element( displacementEqDisc:displacement(), "DispElem")
--vtk:select_element( displacementEqDisc:divergence(), "DivElem")
--vtk:select_element( flowEqDisc:gradient()dbgSolver, "GradP")
--vtk:select(massLinker, "Mass")

-- Init error estimator.
local biotErrorEst
if (problem.error_estimator) then
    biotErrorEst = problem:error_estimator()
else
    biotErrorEst = util.biot.CreateDefaultErrorEst(dim)
end



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Solve transient (linear) problem
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------


local newtonCheck = ConvCheck()
newtonCheck:set_maximum_steps(10)
newtonCheck:set_minimum_defect(1e-10)
newtonCheck:set_reduction(5e-6)
newtonCheck:set_verbose(true)

local newtonCheck2 = CompositeConvCheck(approxSpace)
newtonCheck2:set_component_check("ux", p0*1e-7, 5e-6)
newtonCheck2:set_component_check("uy", p0*1e-7, 5e-6)
if (dim==3) then
    newtonCheck2:set_component_check("uz", p0*1e-7, 5e-6)
end
newtonCheck2:set_component_check("p", p0*1e-9, 5e-6)
newtonCheck2:set_maximum_steps(2)

local newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(lsolver)
newtonSolver:set_convergence_check(newtonCheck)

--[[
local lineSearch = StandardLineSearch();
lineSearch:set_maximum_steps(6)
lineSearch:set_accept_best(true)
newtonSolver:set_line_search(lineSearch)
--]]

--newtonSolver:set_debug(dbgWriter)

local nlsolver = newtonSolver

print(lsolver:config_string())

if (problem.check) then problem:check(u) end

print("Interpolation start values")
problem:interpolate_start_values(u, startTime)

-- Create callback.
function myStepCallback0(u, step, time)
    vtk:print("PoroElasticityInitial.vtu", u, step, time)
    problem:post_processing(u, step, time)
end

print ("Integrating from 0.0 to "..endTime)


--dt =dt*1e-4*problem:get_char_time() -- smaller => more complicated

local charTime = problem:get_char_time()
-- dt = dtFrac*problem:get_char_time()
dt = dtFrac*endTime
dtMin = 1e-2*dt

print("<script_solving>")
if (script_method == "Transient") then


    --
    local myclock = CuckooClock()
    local stepClock = CuckooClock()

    if (( ARGS.LimexNStages > 0)) then
        local dt0 = charTime*1e-50
        print("Computing consistent initial value w/ dt0="..dt0)
        util.SolveNonlinearTimeProblem(u, domainDisc0, nlsolver, myStepCallback0, "PoroElasticityInitial",
                "ImplEuler", 1, startTime, dt0, dt0, dt0, dtRed);
    end



    if ( ARGS.LimexNStages==0) then

        -- Execute linear solver test suite.
        convCheck:set_reduction(1e-10)
        convCheck:set_maximum_steps(500)

        -- local dtTestSet = {1.0, 0.1, 0.01, 1e-3, 1e-4, 1e-6, 1e-8, 0.0}
        local dtTestSet = {1e-8, 0.0}
        for index,dtvalue in ipairs(dtTestSet) do
            dt = dtvalue*charTime
            endTime = dt
            print("%DTFACTOR=\t"..dtvalue.."\tindex=\t"..index)
            problem:interpolate_start_values(u, startTime)
            myclock:tic()
            util.SolveNonlinearTimeProblem(u, domainDiscT, nlsolver, myStepCallback0, "SolverTest"..index,
                    "ImplEuler", 1, startTime, endTime, dt, dt, dtRed);
            print("MYCLOCK="..myclock:cuckoo().."; "..myclock:toc())
        end


    elseif ( ARGS.LimexNStages==1) then
        print("Solving with Euler: dt="..dt.."/"..dtFrac)
        -- STANDARD (implicit Euler) time-stepping.
        local bCheckpointing = false
        --[[
        util.SolveLinearTimeProblem(u, domainDiscT, lsolver, nil, "PoroElasticityTransient",
                         "ImplEuler", 1, startTime, endTime, dt, dtMin, 0.5,
                          bCheckpointing, myStepCallback0)
         --]]
        -- [[
        util.SolveNonlinearTimeProblem(u, domainDiscT, nlsolver, myStepCallback0, "PoroElasticityTransient",
                "ImplEuler", 1, startTime, endTime, dt, dtMin, 0.5, bCheckpointing);
        --]]

    else --  ARGS.LimexNStages > 1
        -- LIMEX time-stepping.

        -- Adjust NEWTON for LIMEX.
        newtonCheck:set_maximum_steps(1)
        newtonCheck:set_supress_unsuccessful(true)

        -- Create & configure LIMEX descriptor.
        local limexDesc = {
            nstages = ARGS.LimexNStages,
            steps = {1,2,3,4,5,6,7,8,9,10},

            domainDisc=domainDiscT,
            nonlinSolver = nlsolver,

            tol = ARGS.LimexTOL,
            dt = dt,
            dtmin = dtMin,
            dtmax = dtMax,

            -- gammaDiscOPT= gammaTensorDisc,  -- no gamma for linear problem
        }


        -- Call factory.
        local limex = util.limex.CreateIntegrator(limexDesc)
        limex:add_error_estimator(biotErrorEst)
        limex:set_stepsize_safety_factor(0.25)
        limex:set_increase_factor(10)
        limex:set_stepsize_greedy_order_factor(0.0)

        limex:disable_matrix_cache()        -- This problem is linear
        --limex:set_time_derivative(udot)   --


        -- Create (& attach) observers.
        if (ARGS.useVTK) then
            local vtkFull = VTKOutput()
            local vtkobserver = VTKOutputObserver("PoroElasticityLimex.vtk", vtkFull)
            limex:attach_observer(vtkobserver)
        end

        local luaobserver = LuaCallbackObserver()

        function myLuaLimexPostProcess(step, time, currdt)
            print ("Time per step :"..stepClock:toc()) -- get time for last step
            local usol=luaobserver:get_current_solution()
            problem:post_processing(usol, step, time)
            stepClock:tic() -- reset timing
            return 0;
        end

        luaobserver:set_callback("myLuaLimexPostProcess")
        limex:attach_observer(luaobserver)



        -- Solve problem using LIMEX.
        myclock:tic()
        limex:apply(u, endTime, u, startTime)
        print("CDELTA="..myclock:toc())
    end

elseif (script_method == "SteadyState") then
        --------------------------------------------------------------------------------
        --------------------------------------------------------------------------------
        --  Solve linear, steady state problem (easy!)
        --------------------------------------------------------------------------------
        --------------------------------------------------------------------------------
    print("<steady_state>")

    local A = MatrixOperator()
    u = GridFunction(approxSpace)
    local b = GridFunction(approxSpace)
    u:set(0.0)
    b:set(0.0)

    -- 1. assemble matrix and rhs
    domainDiscT:assemble_linear(A, b)

    -- 2. set dirichlet values in start iterate
    u:set(0.0)
    domainDiscT:adjust_solution(u)

    -- 3. init solver for linear Operator
    lsolver:init(A)

    SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")

    -- 4. apply solver
    u:set_random(0.0, 1.0)
    lsolver:apply_return_defect(u,b)
    vtk:print("PoroElasticitySteadyState", u, 0, 0.0)

    print("</steady_state>")

elseif (script_method == "Sequential") then
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
    log_norms = Paralog()
    log_norms:set_comm(space_time_communicator)
    log_norms:set_filename("norms")
    log_norms:init()
    biot_norm = BiotBraidDisplacementNorm()
    biot_norm:set_log(log_norms)


    pit_desc.norm = biot_norm --"l2"
    pit_desc.use_residual = util.GetParamNumber("--residual", 0, "string that describes temporal coarsening-factors for every level") == 1
    pit_desc.cfactor = util.GetParam("--c-factor", "4_64", "string that describes temporal coarsening-factors for every level")
    pit_desc.level_num_ref = util.GetParam("--r-factor", "5_5", "string that describes spatial refinement-factors for every level")--  "5_5_5_5"
    pit_desc.hierarchy ={
        max_levels = util.GetParam("--max-level", "10", ""),
    }
    pit_desc.spatial_coarsen_and_refine = false
    pit_desc.store_values = util.GetParamNumber("--store-values", 0, "")
    pit_desc.access_level = util.GetParamNumber("--access-level", 0, "")
    pit_desc.print_level = 3
    pit_desc.print_file = "xbraid_output"
    pit_desc.skip_downcycle_work = util.GetParam("--skip-downcycle", 1, "") == 1
    pit_desc.cycle = util.GetParam("--cycle", "V", "")
    pit_desc.relaxation = util.GetParam("--relaxation", "F", "")
    pit_desc.cf_weight = 1
    util.xbraid.set_domain_disc(pit_desc, pit_inst,domainDiscT)
    pit_inst.spatial_norm = util.xbraid.create_norm(pit_desc, pit_inst)

    print("")

    --ø  pit_desc.logger = "jobb"
    pit_desc.time_interval = {
        start_time = startTime,
        end_time = endTime,
        time_steps = time_steps
    }
    pit_desc.conv_check = {
        max_iterations = util.GetParamNumber("--max-iterations", 5, ""),
        absolut_tol =  1e-50,--2.90e-6 for numref6; 1.16e-5 for numref5; 7.25e-8 for numref7
        relative_tol = 0
    }
    pit_desc.temporal_norm = 3

    pit_desc.driver = {
        name = "Nonlinear",
    }

    pit_desc.integrator = {
        name = "ThetaIntegratorNL",
        nonlinear_solver = nlsolver,
        theta = 1,
    }


    util.xbraid.create_logger(pit_desc,pit_inst)
    util.xbraid.create_integrator(pit_desc, pit_inst)
    util.xbraid.create_driver(pit_desc, pit_inst)
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
    util.xbraid.apply(pit_desc, pit_inst, u)
    print("</xbraid_apply>")
    print("</script_xbraid>")
end  --  doSteadyState


walltime:stop() -- stop timer for entire script execution
if space_time_communicator:get_temporal_rank() == 0 then
    seconds = walltime:get()
    print("<integration time=\""..seconds.."\" unit=\"s\" />")
    total_time = walltime:get()
    t_rank = space_time_communicator:get_temporal_rank()
    print("<summary walltime=\"".. total_time .."\" unit=\"s\" rank=\""..t_rank.."\" />")
end
print("</script_solving>")