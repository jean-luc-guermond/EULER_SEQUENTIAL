
echo ${CMAKE_CURRENT_BINARY_DIR}
cp data_galerkin data
cp reference_regression_test_galerkin reference_regression_test
srun-serial ../../EXECUTABLE/idp_euler_1D_smooth.exe regression
