#include "user_mat.h"

//fixed sized matrix of dimension 1
#define fmat1(type, name, N1) type name##_1D[(N1)]; ViewCMatrixMat <type> name(&name##_1D[0],N1)
//fixed sized matrix of dimension 2
#define fmat2(type, name, N1, N2) type name##_1D[(N1)*(N2)]; ViewCMatrixMat <type> name(&name##_1D[0],N1,N2)


template <size_t N> void matMul(ViewCArrayMat <real_t> a, ViewCArrayMat <real_t> b, ViewCArrayMat <real_t> c)
{
	for (size_t j = 0; j < N; j++) {
		for (size_t i = 0; i < N; i++) {
			real_t dum = 0.0;

			for (size_t k = 0; k < N; k++) {
				dum += a(j, k) * b(k, i);
			}

			c(j, i) += dum;
		}
	}
	return;
}


void rodrigues(ViewCMatrixMat <real_t>& arot, ViewCMatrixMat <real_t>& c) {
	// **********************************************************************
	//     rodrigues --> version 7-2021
	//
	//     builds incremental rotation matrix 'arot' based on rodrigues formu
	//     'c' is the incremental lattice spin.
	//     'arot' transforms from initial to final orientation.
	// **********************************************************************

	fmat2(real_t, th, 3, 3);
	fmat2(real_t, th2, 3, 3);
	fmat1(real_t, v, 3);

	real_t snorm, snorm1;

	//v(1) = c(3, 2);
	//v(2) = c(1, 3);
	//v(3) = c(2, 1);

	v(1) = c(8);
	v(2) = c(3);
	v(3) = c(4);

	snorm = sqrt(v(1) * v(1) + v(2) * v(2) + v(3) * v(3));

	snorm1 = tan(snorm * 0.5);

	if (snorm < 1.e-06) snorm = 1.0;

	for (size_t i = 1; i <= 3; i++) {
		v(i) = snorm1 * v(i) / snorm;
	}

	snorm = v(1) * v(1) + v(2) * v(2) + v(3) * v(3);

	th(1) = 0.0;
	th(1, 2) = -v(3);
	th(1, 3) = v(2);
	th(2, 1) = v(3);
	th(5) = 0.0;
	th(2, 3) = -v(1);
	th(3, 1) = -v(2);
	th(3, 2) = v(1);
	th(9) = 0.0;

	th2 = 0.0;
	matMul<3>(th, th, th2);



	//for (size_t i = 1; i <= 3; i++) {
	//    for (size_t j = 1; j <= 3; j++) {
	//        arot(i, j) = 2. * (th(i, j) + th2(i, j)) / (1.0 + snorm);
	//    }
	//    arot(i, i) += 1;
	//}

	//no longer snorm
	snorm = (2.0 / (1.0 + snorm));

	for (size_t i = 1; i <= 9; i++) {
		arot(i) = (th(i) + th2(i)) * snorm;
	}

	arot(1) += 1.0;
	arot(5) += 1.0;
	arot(9) += 1.0;
	return;
}

int user_material_number_vars() {
	// more complicated calculation of size goes here
	// if necessary
	return 7;
};

//inputs	: velocity gradient(vel_grad), delta time(dt), increment number(kinc), density(den), internal energy(ie)
//outputs	: pressure(p), sound speed(sspd), temperature(temp)
void user_mat_model(real_t* state_vars, real_t* vel_grad, real_t dt, size_t kinc, real_t den, real_t ie, real_t& p, real_t& sspd, real_t& temp) {
	material_t umat;

	ViewCMatrixMat <real_t> L(&vel_grad[0], 3, 3);

	fmat2(real_t, deps, 3, 3);
	fmat2(real_t, dW, 3, 3);
	fmat2(real_t, drot, 3, 3);


	for (size_t i = 1; i <= 3; i++) {
		for (size_t j = 1; j <= 3; j++) {
			deps(i, j) = 0.5 * (L(i, j) + L(j, i)) * dt;
			dW(i, j) = 0.5 * (L(i, j) - L(j, i)) * dt;
		}
	}

	rodrigues(drot, dW);


	umat.read(&state_vars[0]);
	// ---- Enter material model ---- //

	p = (umat.g - 1.0) * ie * den;

	sspd = umat.g * (umat.g - 1.0) * ie;
	if (sspd < umat.csmin) sspd = umat.csmin;
	sspd = sqrt(sspd);

	temp = ie / umat.cv;

	umat.t = temp;
	umat.d = den;

	// ---- Finished with material model ---- //
	umat.write(&state_vars[0]);
	return;
}


