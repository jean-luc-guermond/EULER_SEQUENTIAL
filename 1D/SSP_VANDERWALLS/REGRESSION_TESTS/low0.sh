run="srun-serial ../EXECUTABLE/ssp_euler_1D_vanderwalls.exe regression"

cp data_low_order_case0 data
cp mesh_case0.FEM mesh.FEM
cp reference_regression_test_low_order_case0 reference_regression_test
$run
