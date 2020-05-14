AddTest(
        NAME 2D_RichardsFlow_h_us_quad_ogs5
        PATH Parabolic/Richards
        EXECUTABLE ogs
        EXECUTABLE_ARGS RichardsFlow_2d_compare_ogs5.prj
        TESTER vtkdiff
        DIFF_DATA
        h_us_quad_1000.vtu richards_pcs_0_ts_100_t_100.000000.vtu PRESSURE1 pressure 1e-1 1e-1
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME 2D_RichardsFlow_h_us_quad_small
        PATH Parabolic/Richards
        EXECUTABLE ogs
        EXECUTABLE_ARGS RichardsFlow_2d_small.prj
        TESTER vtkdiff
        DIFF_DATA
        ref_t_1600.000000.vtu richards_pcs_0_ts_1100_t_1600.000000.vtu pressure pressure 1e-8 1e-8
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
        NAME LARGE_2D_RichardsFlow_h_us_quad
        PATH Parabolic/Richards
        RUNTIME 90
        EXECUTABLE ogs
        EXECUTABLE_ARGS RichardsFlow_2d_large.prj
        TESTER vtkdiff
        DIFF_DATA
        ref_t_20000.000000.vtu richards_pcs_0_ts_18200_t_20000.000000.vtu pressure pressure 1e-8 1e-8
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME 2D_RichardsFlow_h_us_quad_small_PID_adaptive_dt
    PATH Parabolic/Richards
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_small_PID_adaptive_dt.prj
    TESTER vtkdiff
    DIFF_DATA
    ref_t_1600.000000.vtu richards_pcs_0_ts_803_t_1600.000000.vtu pressure pressure 1e-8 1e-3
# The following three comparisons are used just to check whether the output is
# made at the fixed times of 50, 100 and 500, which are given in the project
# file of RichardsFlow_2d_small_adaptive_dt.prj
    richards_pcs_0_ts_28_spec_t_50.000000.vtu richards_pcs_0_ts_28_t_50.000000.vtu pressure pressure 1e-10 1e-10
    richards_pcs_0_ts_53_spec_t_100.000000.vtu richards_pcs_0_ts_53_t_100.000000.vtu pressure pressure 1e-10 1e-10
    richards_pcs_0_ts_253_spec_t_500.000000.vtu richards_pcs_0_ts_253_t_500.000000.vtu pressure pressure 1e-10 1e-10
    REQUIREMENTS NOT OGS_USE_MPI
)

AddTest(
    NAME 2D_RichardsFlow_h_us_quad_small_iteration_adaptive_dt
    PATH Parabolic/Richards
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_small_iteration_adaptive_dt.prj
    REQUIREMENTS NOT OGS_USE_MPI
    RUNTIME 51
    # No vtkdiff comparison here, because of the different file names for
    # different machines, which again is due to the adaptive time stepping
    # scheme. When the output file format can be specified in the project files,
    # e.g. in the form Richards_%t where %t is the current time, the output will
    # no longer be ambiguous.
)

# Comparison test for richards mechanics w/o deformations.
AddTest(
    NAME RichardsMechanics_RichardsFlow_2d_richardsflow
    PATH RichardsMechanics
    EXECUTABLE ogs
    EXECUTABLE_ARGS RichardsFlow_2d_richardsflow.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB RichardsFlow_2d_richardsflow_pcs_0_ts_*.vtu pressure pressure 5e-11 1e-15
    GLOB RichardsFlow_2d_richardsflow_pcs_0_ts_*.vtu saturation saturation 1e-14 1e-15
)

AddTest(
    NAME RichardsFlow_square_1e2_flow_fully_saturated
    PATH Parabolic/Richards
    EXECUTABLE ogs
    EXECUTABLE_ARGS flow_fully_saturated.prj
    WRAPPER time
    TESTER vtkdiff
    REQUIREMENTS NOT OGS_USE_MPI
    DIFF_DATA
    GLOB flow_fully_saturated_pcs_0_ts_*.vtu pressure pressure 1e-9 0
    GLOB flow_fully_saturated_pcs_0_ts_*.vtu darcy_velocity darcy_velocity 1e-11 1e-15
)
#PETSc/MPI
#AddTest(
#    NAME 2D_RichardsFlow_h_us_quad_small_PID_adaptive_dt
#    PATH Parabolic/Richards
#    EXECUTABLE_ARGS RichardsFlow_2d_small_PID_adaptive_dt.prj
#    WRAPPER mpirun
#    WRAPPER_ARGS -np 1
#    TESTER vtkdiff
#    REQUIREMENTS OGS_USE_MPI
#    RUNTIME 220
#    DIFF_DATA
#    ref_t_1600.000000.vtu richards_pcs_0_ts_803_t_1600_000000_0.vtu pressure pressure 1e-8 1e-3
#)

# Comparison test for richards mechanics w/o deformations.
AddTest(
    NAME Parallel_RichardsMechanics_RichardsFlow_2d_richardsflow
    PATH RichardsMechanics
    EXECUTABLE_ARGS RichardsFlow_2d_richardsflow.prj
    WRAPPER mpirun
    WRAPPER_ARGS -np 1
    TESTER vtkdiff
    REQUIREMENTS OGS_USE_MPI
    DIFF_DATA
    RichardsFlow_2d_richardsflow_pcs_0_ts_99_t_1900.000000.vtu RichardsFlow_2d_richardsflow_pcs_0_ts_99_t_1900_000000_0.vtu pressure pressure 5e-8 1e-10
    RichardsFlow_2d_richardsflow_pcs_0_ts_99_t_1900.000000.vtu RichardsFlow_2d_richardsflow_pcs_0_ts_99_t_1900_000000_0.vtu saturation saturation 1e-10 1e-11
)

# 1d infiltration with RichardsFlow
AddTest(
    NAME 1d-infiltration
    PATH Parabolic/Richards/infiltration_1D
    EXECUTABLE ogs
    EXECUTABLE_ARGS 1d-infiltration.prj
    TESTER vtkdiff
    DIFF_DATA
    out_ex1_pcs_0_ts_0_t_0.000000_expected.vtu out_ex1_ts_0.vtu pressure pressure 1e-6 1e-6
    out_ex1_pcs_0_ts_50_t_2303.496783_expected.vtu out_ex1_ts_50.vtu pressure pressure 1e-6 1e-6
    out_ex1_pcs_0_ts_100_t_7303.496783_expected.vtu out_ex1_ts_100.vtu pressure pressure 1e-6 1e-6
    out_ex1_pcs_0_ts_150_t_12303.496783_expected.vtu out_ex1_ts_150.vtu pressure pressure 1e-6 1e-6
    out_ex1_pcs_0_ts_200_t_17303.496783_expected.vtu out_ex1_ts_200.vtu pressure pressure 1e-6 1e-6
    out_ex1_pcs_0_ts_250_t_22303.496783_expected.vtu out_ex1_ts_250.vtu pressure pressure 1e-6 1e-6
    out_ex1_pcs_0_ts_300_t_27303.496783_expected.vtu out_ex1_ts_300.vtu pressure pressure 1e-6 1e-6
    out_ex1_pcs_0_ts_350_t_32303.496783_expected.vtu out_ex1_ts_350.vtu pressure pressure 1e-4 1e-5
    out_ex1_pcs_0_ts_400_t_37303.496783_expected.vtu out_ex1_ts_400.vtu pressure pressure 1e-4 1e-5
    out_ex1_pcs_0_ts_450_t_42303.496783_expected.vtu out_ex1_ts_450.vtu pressure pressure 1e-4 1e-5
    out_ex1_pcs_0_ts_500_t_47303.496783_expected.vtu out_ex1_ts_500.vtu pressure pressure 1e-4 1e-5
    out_ex1_pcs_0_ts_550_t_52303.496783_expected.vtu out_ex1_ts_550.vtu pressure pressure 1e-4 1e-5
    out_ex1_pcs_0_ts_600_t_57303.496783_expected.vtu out_ex1_ts_600.vtu pressure pressure 1e-4 1e-5
    out_ex1_pcs_0_ts_650_t_62293.496783_expected.vtu out_ex1_ts_650.vtu pressure pressure 1e-4 1e-5
    out_ex1_pcs_0_ts_700_t_67293.496783_expected.vtu out_ex1_ts_700.vtu pressure pressure 1e-4 1e-5
    out_ex1_pcs_0_ts_750_t_72293.496783_expected.vtu out_ex1_ts_750.vtu pressure pressure 1e-4 1e-5
    out_ex1_pcs_0_ts_800_t_77181.778033_expected.vtu out_ex1_ts_800.vtu pressure pressure 1e-4 1e-5
    out_ex1_pcs_0_ts_850_t_82181.778033_expected.vtu out_ex1_ts_850.vtu pressure pressure 1e-4 1e-5
    out_ex1_pcs_0_ts_894_t_86400.000000_expected.vtu out_ex1_ts_894.vtu pressure pressure 1e-4 1e-5
    REQUIREMENTS NOT OGS_USE_MPI
)
