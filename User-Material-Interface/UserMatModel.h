#pragma once
#include <stdio.h> //size_t definition
#include <memory> //smart pointers

// UserMatModel is a pure virual class that when inherited,
// the solve method must be defined.
class UserMatModel
{
public:
    UserMatModel() = default;
    virtual ~UserMatModel() = default;
    virtual void solve(double* vel_grad, double* stress, double dt, size_t cycle, size_t elem_gid) = 0;
    virtual void write_texture() = 0;

};

void init_user_mat_model(std::shared_ptr<UserMatModel>* elem_user_mat_model,
                         const size_t* elem_mat_id,
                         const size_t num_elems);


