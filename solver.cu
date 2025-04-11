
#include <cusolverSp.h>
#include <cusparse.h>
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cublas_v2.h>
#include <device_launch_parameters.h>

#include "solver.h"


__global__ void cuSpPrecondition(double* diagonal_inv, double* data, long long int n) {
    long long int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < n) {
        double d = data[idx];
        if (diagonal_inv[idx] != 0)
            data[idx] = d / diagonal_inv[idx];
        else
            data[idx] = d * 0.1;
    }
}

__global__ void cuProject(double* data, long long int n) {
    long long int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx < n) {
        double d = data[idx];
        if (d < 0)
            data[idx] = 0;
    }
}

void _amg_gmres(double* csrval, int* csrcol, int* csrrow, int n, int nnz, double* yb, double* dx) {
    AMGX_SAFE_CALL(AMGX_initialize());
    AMGX_SAFE_CALL(AMGX_initialize_plugins());
    AMGX_config_handle cfg;
    AMGX_SAFE_CALL(AMGX_config_create_from_file(&cfg, "E:/code/NanoporePDE/configs/FGMRES_CLASSICAL_AGGRESSIVE_PMIS.json"));
    AMGX_resources_handle rsrc;
    AMGX_SAFE_CALL(AMGX_resources_create_simple(&rsrc, cfg));
    
    AMGX_matrix_handle A;
    AMGX_vector_handle x;
    AMGX_vector_handle b;
    AMGX_solver_handle solver;

    AMGX_SAFE_CALL(AMGX_matrix_create(&A, rsrc, AMGX_mode_dDDI));
    AMGX_SAFE_CALL(AMGX_vector_create(&x, rsrc, AMGX_mode_dDDI));
    AMGX_SAFE_CALL(AMGX_vector_create(&b, rsrc, AMGX_mode_dDDI));
    AMGX_SAFE_CALL(AMGX_solver_create(&solver, rsrc, AMGX_mode_dDDI, cfg));

    AMGX_SAFE_CALL(AMGX_matrix_upload_all(A, n, nnz, 1, 1, csrrow, csrcol, csrval, nullptr));
    AMGX_SAFE_CALL(AMGX_vector_upload(b, n, 1, yb));
    AMGX_SAFE_CALL(AMGX_vector_upload(x, n, 1, dx));
    AMGX_SAFE_CALL(AMGX_solver_setup(solver, A));
    AMGX_SAFE_CALL(AMGX_solver_solve(solver, b, x));
    AMGX_SAFE_CALL(AMGX_vector_download(x, dx));
}

void _tfqmr(double* csrval, long long int* csrcol, long long int* csrrow, void** d_buffer, double* diag, long long int n, long long int nnz, bool* buffer,
    double* y, double* x, double* yk1, double* yk2, double* r0, double* uk1, double* uk2, double* wk, double* vk, double* rk, double* dk) {
    //sove J*dv = -F;
    dim3 griddim((n + 255) / 256);
    dim3 blockdim(256);
    double epsi = 0.0;
    double theta = 0.0;
    double rho;
    double tau;
    double normal_b;
    double rho_new;
    double const_none = -1.0;
    double const_one = 1.0;
    double const_zero = 0.0;
    cublasHandle_t handle;
    cublasCreate(&handle);

    cusparseHandle_t sphandle;
    cusparseSpMatDescr_t matJ;
    cusparseDnVecDescr_t vecX, vecY;
    cusparseCreate(&sphandle);
    cusparseCreateDnVec(&vecX, n, x, CUDA_R_64F); // x
    cusparseCreateDnVec(&vecY, n, r0, CUDA_R_64F); // y 
    cusparseCreateCsr(&matJ, n, n, nnz, csrrow, csrcol, csrval, CUSPARSE_INDEX_64I, CUSPARSE_INDEX_64I, CUSPARSE_INDEX_BASE_ZERO, CUDA_R_64F); // J
    size_t bufferSize = 0;
    if (*buffer) {
        cusparseSpMV_bufferSize(sphandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            &const_none, matJ, vecX, &const_one, vecY, CUDA_R_64F, CUSPARSE_MV_ALG_DEFAULT, &bufferSize);
        cudaMalloc(d_buffer, bufferSize);
        *buffer = false;
    }
    //cudaMemset(x, 0.5, n * sizeof(double));
    cublasDcopy(handle, n, y, 1, r0, 1);
    cudaDeviceSynchronize();
    cuSpPrecondition << <griddim, blockdim >> > (diag, y, n);
    cudaDeviceSynchronize();
    cublasDnrm2(handle, n, y, 1, &normal_b);
    /* ==================start solve the linear equations ================*/
    // calculte r0 = yb - Ax 
    cusparseSpMV(sphandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
        &const_none, matJ, vecX, &const_one, vecY, CUDA_R_64F, CUSPARSE_MV_ALG_DEFAULT, *d_buffer);
    cudaDeviceSynchronize();
    cuSpPrecondition << <griddim, blockdim >> > (diag, r0, n);
    cudaDeviceSynchronize();
    cublasDcopy(handle, n, r0, 1, yk1, 1);
    cublasDcopy(handle, n, r0, 1, wk, 1);
    cublasDcopy(handle, n, r0, 1, rk, 1);
    cudaMemset(dk, 0, n * sizeof(double));
    cusparseDnVecSetValues(vecX, yk1);
    cusparseDnVecSetValues(vecY, uk1);
    cusparseSpMV(sphandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
        &const_one, matJ, vecX, &const_zero, vecY, CUDA_R_64F, CUSPARSE_MV_ALG_DEFAULT, *d_buffer);
    cudaDeviceSynchronize();
    cuSpPrecondition << <griddim, blockdim >> > (diag, uk1, n);
    cudaDeviceSynchronize();
    cublasDcopy(handle, n, uk1, 1, vk, 1);
    cublasDdot(handle, n, r0, 1, r0, 1, &rho);
    cublasDnrm2(handle, n, rk, 1, &tau);
    cudaDeviceSynchronize();
    // start;
    int k = 0;
    bool terminate = false;
    double residual;

    while (k < 20000 && !terminate) {
        k++;
        // 1. sigma = <r0,v>; alpha = rho / simga; y2 = y1 - alpha*v; u2 = Ay2;
        double sigma;
        cublasDdot(handle, n, r0, 1, vk, 1, &sigma);
        cudaDeviceSynchronize();
        double alpha = rho / sigma;
        double temporary = -alpha;
        cublasDcopy(handle, n, yk1, 1, yk2, 1);
        cublasDaxpy(handle, n, &temporary, vk, 1, yk2, 1);
        cudaDeviceSynchronize();
        cusparseDnVecSetValues(vecX, yk2);
        cusparseDnVecSetValues(vecY, uk2);
        cusparseSpMV(sphandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            &const_one, matJ, vecX, &const_zero, vecY, CUDA_R_64F, CUSPARSE_MV_ALG_DEFAULT, *d_buffer);
        cudaDeviceSynchronize();
        cuSpPrecondition << <griddim, blockdim >> > (diag, uk2, n);
        cudaDeviceSynchronize();
        //2. odd or even
        for (int j = 1; j <= 2; j++) {
            int m = 2 * k - 2 + j;
            double coeff;
            if (j == 1) {
                cublasDaxpy(handle, n, &temporary, uk1, 1, wk, 1);
                coeff = theta * theta * epsi / alpha;
                cublasDscal(handle, n, &coeff, dk, 1);
                cublasDaxpy(handle, n, &const_one, yk1, 1, dk, 1);
            }
            else {
                cublasDaxpy(handle, n, &temporary, uk2, 1, wk, 1);
                coeff = theta * theta * epsi / alpha;
                cublasDscal(handle, n, &coeff, dk, 1);
                cublasDaxpy(handle, n, &const_one, yk2, 1, dk, 1);
            }
            cublasDnrm2(handle, n, wk, 1, &theta);
            cudaDeviceSynchronize();
            theta = theta / tau;
            double c = 1 / sqrt(1 + theta * theta);
            tau = tau * theta * c;
            epsi = c * c * alpha;
            cublasDaxpy(handle, n, &epsi, dk, 1, x, 1);
            cudaDeviceSynchronize();
            residual = tau * sqrt(double(m + 1)) / normal_b;
            if (k % 1 == 0 && j == 1)
                std::cout << "step: " << k << " tfres: " << residual << " rb: " << normal_b << std::endl;
            if (residual <= 1e-4 || rho == 0) {
                terminate = true;
                if (rho == 0)
                    std::cout << "rho 0 stop" << std::endl;
                else {
                    std::cout << "step: " << k << " tfres: " << residual << " rb: " << normal_b << std::endl;
                }
                break;
            }
        }
        //3. rho = <r0, w>; beta = rho / rho_old
        double rho_new;
        cublasDdot(handle, n, r0, 1, wk, 1, &rho_new);
        cudaDeviceSynchronize();
        double beta = rho_new / rho;
        rho = rho_new;
        //4. y1 = w + beta*y2, u1 = Ay1;
        cublasDcopy(handle, n, yk2, 1, yk1, 1);
        cublasDscal(handle, n, &beta, yk1, 1);
        cublasDaxpy(handle, n, &const_one, wk, 1, yk1, 1);
        cudaDeviceSynchronize();
        cusparseDnVecSetValues(vecX, yk1);
        cusparseDnVecSetValues(vecY, uk1);
        cusparseSpMV(sphandle, CUSPARSE_OPERATION_NON_TRANSPOSE,
            &const_one, matJ, vecX, &const_zero, vecY, CUDA_R_64F, CUSPARSE_MV_ALG_DEFAULT, *d_buffer);
        cudaDeviceSynchronize();
        cuSpPrecondition << <griddim, blockdim >> > (diag, uk1, n);
        cudaDeviceSynchronize();
        //5. v = u1 + beta(u2 + beta*v);
        cublasDscal(handle, n, &beta, vk, 1);
        cublasDaxpy(handle, n, &const_one, uk2, 1, vk, 1);
        cublasDscal(handle, n, &beta, vk, 1);
        cublasDaxpy(handle, n, &const_one, uk1, 1, vk, 1);
        cudaDeviceSynchronize();
    }
    //std::cout << k << " " << residual << std::endl;
    cusparseDestroyDnVec(vecX);
    cusparseDestroyDnVec(vecY);
    cusparseDestroySpMat(matJ);
    cublasDestroy(handle);
    cusparseDestroy(sphandle);
}


void PNPNSSolver::solve(PNPNS::CSRMatrix& A, double* x, double* yc, long long int nn, long long int nnz) {
    /*
    allocate the memory on gpu for tfqmr
    */
    double* csrval_gpu;
    long long int* csrcol_gpu;
    long long int* csrrow_gpu;
    double* x_gpu;
    double* yc_gpu;
    double* yk1_gpu;
    double* yk2_gpu;
    double* r0_gpu;
    double* uk1_gpu;
    double* uk2_gpu;
    double* wk_gpu;
    double* vk_gpu;
    double* rk_gpu;
    double* dk_gpu;
    double* diaginv_gpu;

    cudaMalloc((void**)&csrval_gpu, nnz * sizeof(double));
    cudaMalloc((void**)&csrcol_gpu, nnz * sizeof(long long  int));
    cudaMalloc((void**)&csrrow_gpu, (nn + 1) * sizeof(long long int));
    cudaMalloc((void**)&x_gpu, nn * sizeof(double));
    cudaMalloc((void**)&yc_gpu, nn * sizeof(double));
    cudaMalloc((void**)&yk1_gpu, nn * sizeof(double));
    cudaMalloc((void**)&yk2_gpu, nn * sizeof(double));
    cudaMalloc((void**)&r0_gpu, nn * sizeof(double));
    cudaMalloc((void**)&uk1_gpu, nn * sizeof(double));
    cudaMalloc((void**)&uk2_gpu, nn * sizeof(double));
    cudaMalloc((void**)&wk_gpu, nn * sizeof(double));
    cudaMalloc((void**)&vk_gpu, nn * sizeof(double));
    cudaMalloc((void**)&rk_gpu, nn * sizeof(double));
    cudaMalloc((void**)&dk_gpu, nn * sizeof(double));
    cudaMalloc((void**)&diaginv_gpu, nn * sizeof(double));

    cudaMemcpy((void*)csrrow_gpu, A.outerIndexPtr(), (nn + 1) * sizeof(long long int), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)csrcol_gpu, A.innerIndexPtr(), nnz * sizeof(long long int), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)csrval_gpu, A.valuePtr(), nnz * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)yc_gpu, (void*)yc, nn * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy((void*)x_gpu, (void*)x, nn * sizeof(double), cudaMemcpyHostToDevice);

    Eigen::VectorXd diagonal;
    diagonal = A.diagonal();
    cudaMemcpy((void*)diaginv_gpu, (void*)diagonal.data(), nn * sizeof(double), cudaMemcpyHostToDevice);
    void* d_buffer = NULL;
    bool isallocatebuffer = true;
    /*----------------------------------------------------------------*/
    //_tfqmr(csrval_gpu, csrcol_gpu, csrrow_gpu, &d_buffer, diaginv_gpu, nn, nnz, &isallocatebuffer,
    //    yc_gpu, x_gpu, yk1_gpu, yk2_gpu, r0_gpu, uk1_gpu, uk2_gpu, wk_gpu, vk_gpu, rk_gpu, dk_gpu);
    int* col_data = new int[nnz];
    int* row_data = new int[nn + 1];
    for (int i = 0; i <= nn; i++) {
        row_data[i] = A.outerIndexPtr()[i];
    }
    for (int i = 0; i < nnz; i++) {
        col_data[i] = A.innerIndexPtr()[i];
    }
    
    _amg_gmres(A.valuePtr(), col_data, row_data, nn, nnz, yc, x);
    //cudaMemcpy((void*)x, x_gpu, nn * sizeof(double), cudaMemcpyDeviceToHost);
    /*free memory*/
    cudaFree(csrval_gpu);
    cudaFree(csrcol_gpu);
    cudaFree(csrrow_gpu);
    cudaFree(diaginv_gpu);
    cudaFree(x_gpu);
    cudaFree(yc_gpu);
    cudaFree(yk1_gpu);
    cudaFree(yk2_gpu);
    cudaFree(r0_gpu);
    cudaFree(uk1_gpu);
    cudaFree(uk2_gpu);
    cudaFree(wk_gpu);
    cudaFree(vk_gpu);
    cudaFree(rk_gpu);
    cudaFree(dk_gpu);
    cudaFree(d_buffer);
}
