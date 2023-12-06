#include <gtest/gtest.h>
#include "matar.h"
#include "Geometry.h"
#include <cassert>
#include <set>
#include <iostream>
#include <string>
#include <cstdlib>

// If you change these, you will have to change like most of the test functions.
const static size_t N = 11;
const static double length = 1.0, width = 1.0, height = 1.0;

size_t indices_to_index(size_t i, size_t j, size_t k) {
    return i + j * N + k * N * N;
}
void index_to_indices(size_t index, size_t& i, size_t& j, size_t& k) {
    k = index / (N * N);
    j = (index % (N * N)) / N;
    i = (index % N);
}

mtr::DCArrayKokkos<double> make_regular_grid() {
    mtr::DCArrayKokkos<double> m(N * N * N, 3);
    FOR_ALL(i, 0, N, j, 0, N, k, 0, N, {
        size_t index = indices_to_index(i, j, k);
        m(index, 0) = i * length / (N - 1);
        m(index, 1) = j * width  / (N - 1);
        m(index, 2) = k * height / (N - 1);
    });

    return m;
}

mtr::DCArrayKokkos<bool> paint_points(const mtr::DCArrayKokkos<double> points, const volume_t* v) {
    mtr::DCArrayKokkos<bool> contained(points.dims(0));

    for (size_t i = 0; i < points.dims(0); i++) {
        const double p[3] = {points(i, 0), points(i, 1), points(i, 2)};
        contained(i) = v->contains(p);
    }
    return contained;
}

void assert_masked_true(mtr::DCArrayKokkos<bool> mask, const std::vector<std::vector<size_t>> indices) {
    std::set<size_t> true_indices;
    for (auto ijk : indices)
        true_indices.insert(indices_to_index(ijk[0], ijk[1], ijk[2]));
    
    size_t i, j, k;
    for (size_t index = 0; index < mask.size(); index++) {
        index_to_indices(index, i, j, k);
        EXPECT_EQ(mask(index), (true_indices.count(index) > 0)) << "(" << i << "," << j << "," << k << ")";
    }
}

void assert_masked_false(mtr::DCArrayKokkos<bool> mask, const std::vector<std::vector<size_t>> indices) {
    std::set<size_t> false_indices;
    for (auto ijk : indices)
        false_indices.insert(indices_to_index(ijk[0], ijk[1], ijk[2]));
    
    size_t i, j, k;
    for (size_t index = 0; index < mask.size(); index++) {
        index_to_indices(index, i, j, k);
        EXPECT_EQ(mask(index), (false_indices.count(index) == 0)) << "(" << i << "," << j << "," << k << ")";
    }
}

class RegionInclusionTest : public testing::Test {
protected:
    RegionInclusionTest() { }

    mtr::DCArrayKokkos<double> points;
    void SetUp() override {
        points = make_regular_grid();
    }
};

TEST_F(RegionInclusionTest, Global) {
    auto v = global();
    auto mask = paint_points(points, &v);

    for (size_t i = 0; i < mask.size(); i++)
        assert(mask(i));
}

TEST_F(RegionInclusionTest, SphereGlobal) {
    double o[3] {0.5, 0.5, 0.5};
    auto v = spherical_shell(0, 2, o);
    auto mask = paint_points(points, &v);

    for (size_t i = 0; i < mask.size(); i++)
        assert(mask(i));
}

TEST_F(RegionInclusionTest, SphereCorner) {
    double o[3] {1.0, 1.0, 1.0};
    auto v = spherical_shell(0, 0.01, o);
    auto mask = paint_points(points, &v);

    assert_masked_true(mask, {{10, 10, 10}});
}

TEST_F(RegionInclusionTest, SphereCornerShell) {

    double o[3] {1.0, 1.0, 1.0};
    auto v = spherical_shell(0.01, 0.19, o);
    auto mask = paint_points(points, &v);

    assert_masked_true(
        mask, 
        {   // These are the top right corner indices right around the corner
            {9, 10, 10},
            {10, 9, 10},
            {10, 10, 9},
            {9, 9, 10},
            {9, 10, 9},
            {10, 9, 9},
            {9, 9, 9}
        }
    );
}

TEST_F(RegionInclusionTest, CylindricalCenter) {

    double o[3] {0.5, 0.5, 0.5};
    auto v = cylindrical_shell(0, 0.05, 0.51, o);
    auto mask = paint_points(points, &v);

    assert_masked_true(
        mask, 
        {   // These are the top right corner indices right around the corner
            {5, 5, 0},
            {5, 5, 1},
            {5, 5, 2},
            {5, 5, 3},
            {5, 5, 4},
            {5, 5, 5},
            {5, 5, 6},
            {5, 5, 7},
            {5, 5, 8},
            {5, 5, 9},
            {5, 5, 10},
        }
    );
}

TEST_F(RegionInclusionTest, CylindricalCenterShell) {
    double o[3] {0.5, 0.5, 0.5};
    auto v = cylindrical_shell(0.05, 2, 0.51, o);
    auto mask = paint_points(points, &v);

    assert_masked_false(
        mask, 
        {   // These are the top right corner indices right around the corner
            {5, 5, 0},
            {5, 5, 1},
            {5, 5, 2},
            {5, 5, 3},
            {5, 5, 4},
            {5, 5, 5},
            {5, 5, 6},
            {5, 5, 7},
            {5, 5, 8},
            {5, 5, 9},
            {5, 5, 10},
        }
    );
}

TEST_F(RegionInclusionTest, BoxLower) {
    auto v = box(0, 0.5, 0, 0.5, 0, 0.5);
    auto mask = paint_points(points, &v);

    std::vector<std::vector<size_t>> indices;
    for (size_t i = 0; i <= 5; i++)
        for (size_t j = 0; j <= 5; j++)
            for (size_t k = 0; k <= 5; k++)
                indices.push_back({i, j, k});
    
    assert_masked_true(mask, indices);
}

TEST_F(RegionInclusionTest, VoxelGridRandom) {
    mtr::DCArrayKokkos<bool> grid(N - 1, N - 1, N - 1, "grid");
    std::srand(1);
    std::vector<std::vector<size_t>> indices;
    for (size_t i = 0; i < N-1; i++) {
        for (size_t j = 0; j < N-1; j++) {
            for (size_t k = 0; k < N-1; k++) {
                bool filled = (double)std::rand() / RAND_MAX > 0.5;
                grid(i, j, k) = filled;
                if (filled)
                    indices.push_back({i, j, k});
            }
        }
    }

    double sizes[3] = { length / (N-1), width / (N-1), height / (N-1) };
    // Slide the grid so there is only 1 point per voxel
    // And the indices of the points line up with `indices`
    double o[3] = { -sizes[0] / 2, -sizes[1] / 2, -sizes[2] / 2 }; 
    auto v = voxel_grid(sizes, o, grid);

    auto mask = paint_points(points, &v);
    assert_masked_true(mask, indices);
}


TEST_F(RegionInclusionTest, VoxelGridRandomVolume) {
    mtr::DCArrayKokkos<bool> grid(N - 1, N - 1, N - 1, "grid");
    std::srand(1);

    double volume = 0;

    for (size_t i = 0; i < N-1; i++) {
        for (size_t j = 0; j < N-1; j++) {
            for (size_t k = 0; k < N-1; k++) {
                bool filled = (double)std::rand() / RAND_MAX > 0.5;
                grid(i, j, k) = filled;
                if (filled)
                    volume += 1.0;
            }
        }
    }

    double sizes[3] = { length / (N-1), width / (N-1), height / (N-1) };
    // Slide the grid so there is only 1 point per voxel
    // And the indices of the points line up with `indices`
    double o[3] = { -sizes[0] / 2, -sizes[1] / 2, -sizes[2] / 2 }; 
    auto v = voxel_grid(sizes, o, grid);
    v.cache_volume();
    double cell_volume = sizes[0] * sizes[1] * sizes[2];
    EXPECT_DOUBLE_EQ(volume * cell_volume, v.volume);
}


int main(int argc, char **argv) {
    Kokkos::initialize();
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}