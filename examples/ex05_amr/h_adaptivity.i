[Mesh]
  file = square.e

  # Start out with an initial refinement so we have the option to coarsen back
  # if desired.  The mesh with not coarsen more than the initial mesh (after this step)
  uniform_refine = 1
[]

[Variables]
  active = 'convected diffused'

  [./convected]
    order = FIRST
    family = LAGRANGE
  [../]

  [./diffused]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[Kernels]
  active = 'example_diff conv diff'

  [./example_diff]
    type = Diffusion
    variable = convected
  [../]

  [./conv]
    type = Convection
    variable = convected
    some_variable = diffused
  [../]

  [./diff]
    type = Diffusion
    variable = diffused
  [../]
[]

[BCs]
  active = 'left_convected right_convected left_diffused right_diffused'

  [./left_convected]
    type = DirichletBC
    variable = convected
    boundary = '1'
    value = 0
  [../]

  [./right_convected]
    type = DirichletBC
    variable = convected
    boundary = '2'
    value = 1

    some_var = diffused
  [../]

  [./left_diffused]
    type = DirichletBC
    variable = diffused
    boundary = '1'
    value = 0
  [../]

  [./right_diffused]
    type = DirichletBC
    variable = diffused
    boundary = '2'
    value = 10
  [../]

[]

[Executioner]
  type = Steady
  petsc_options = '-snes_mf_operator'

  # The adapativity block
  [./Adaptivity]
    steps = 7
    refine_fraction = 0.3
    coarsen_fraction = 0
    max_h_level = 10
    error_estimator = KellyErrorEstimator
    print_changed_info = true
  [../]
[]

[Output]
  file_base = out
  interval = 1
  exodus = true
  perf_log = true
[]
   
    
