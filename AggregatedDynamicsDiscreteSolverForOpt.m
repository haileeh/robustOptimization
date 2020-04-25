function OutputMetric = AggregatedDynamicsDiscreteSolverForOpt(Target,Chaser,tspan,deltaT,term_type)

Target.deltaT = deltaT;

[out,~] = AggregatedDynamicsDiscreteSolver(Target,Chaser,tspan,term_type);

OutputMetric = out.int.TotalLM + out.TotalDeltaV_mps;

end