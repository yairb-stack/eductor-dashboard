const ExposureCompute=ExposureRuntime.create({model:Exposure,buildHash:EDUCTOR_BUILD.codeSHA256,wallData:key=>EDUCTOR_WALLS[key]});
// Classifies failed calculations and checks alternatives with the SAME solver.
// Suggestions never change inputs unless the user selects one in the UI.
