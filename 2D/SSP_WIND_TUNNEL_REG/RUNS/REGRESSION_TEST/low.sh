
echo ${CMAKE_CURRENT_BINARY_DIR}
cp data_low_order data
cp reference_regression_test_low_order reference_regression_test
srun-serial ../../EXECUTABLE/ssp_euler.exe regression
