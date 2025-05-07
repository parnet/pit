
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

XBraidParameters.level_num_ref = util.xbraid.numref_from_string(XBraidParameters.level_numref)





if XBraidParameters.redirect then
    repl = ReplaceStandardStream()
    repl:set_space_time_comm(space_time_communicator)
    repl:apply()
end

log_job = Paralog()
log_job:set_comm(space_time_communicator)
log_job:set_filename("job")
log_job:init()

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



if XBraidParameters.redirect then
    repl:undo() -- give back the std::cout to terminal
end
