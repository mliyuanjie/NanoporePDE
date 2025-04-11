#pragma once
#include "operator.h" 
#include <amgx_c.h>

#define AMGX_SAFE_CALL_WITH_DEBUG(call)              \
    {                                                \
        AMGX_RC err = call;                          \
        if (err != AMGX_RC_OK) {                     \
            const char *msg;                         \
            AMGX_get_error_string(err, &msg);        \
            std::cerr << "AMGX error: " << msg << std::endl; \
        }                                            \
    }

class PNPNSSolver {
public:
	PNPNSSolver() {};
	void solve(PNPNS::CSRMatrix& A, double* x, double* yc, long long int nn, long long int nnz);
};
