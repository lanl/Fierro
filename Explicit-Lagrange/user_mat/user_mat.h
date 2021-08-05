#ifndef USER_MAT_H
#define USER_MAT_H

#include <iostream>
#include <math.h>
#include "utilities.h"
#include "matar.h"

using namespace utils;

int  user_material_number_vars();
void user_mat_model(real_t* state_vars, real_t* vel_grad, real_t dt, size_t kinc, real_t den, real_t ie, real_t& p, real_t& sspd, real_t& temp);


class material_t {

public:

	real_t cv;    // specific heat
	real_t t;     // temperature
	real_t d;     // density
	real_t en;    // specific internal energy
	real_t g;     // gamma
	real_t csmin; // minimum sound speed
	real_t b1;    // linear coefficient in Riemann solver

	// ViewCArray <double> mat_arr; // example syntax for arrays

	material_t& read(real_t* state_vars) {

		int cnt = 0;
		cv = state_vars[cnt];
		
		cnt++;
		t = state_vars[cnt];
		cnt++;
		d = state_vars[cnt];
		cnt++;
		en = state_vars[cnt];
		cnt++;
		g = state_vars[cnt];
		cnt++;
		csmin = state_vars[cnt];
		cnt++;
		b1 = state_vars[cnt];
		cnt++;

		// mat_arr.set(ViewCArray <double>(&state_vars[cnt], 3, 3, 3, 3));
		// cnt += mat_arr.size();

		return *this;

	};

	material_t& write(real_t* state_vars) {
		int cnt = 0;
		state_vars[cnt] = cv;
		cnt++;
		state_vars[cnt] = t;
		cnt++;
		state_vars[cnt] = d;
		cnt++;
		state_vars[cnt] = en;
		cnt++;
		state_vars[cnt] = g;
		cnt++;
		state_vars[cnt] = csmin;
		cnt++;
		state_vars[cnt] = b1;
		cnt++;

		// no need to write matrices as they are still 
		// pointing to the original state_vars

		return *this;

	};

	material_t& init(real_t* state_vars) {
		read(&state_vars[0]);

		cv = 1.0;
		g = (5.0 / 3.0);
		csmin = 1.0e-14;
		b1 = 1.3333;

		t = 0.0;
		d = 0.0;
		en = 0.0;


		write(&state_vars[0]);
		return *this;
	};

	~material_t() {

	}

};




#endif