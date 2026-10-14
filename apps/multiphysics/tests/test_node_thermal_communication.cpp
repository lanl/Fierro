#include "state.hpp"

int main(int argc, char** argv)
{
    int provided = 0;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &provided);
    if (provided < MPI_THREAD_FUNNELED) {
        MPI_Finalize();
        return 1;
    }
    Kokkos::initialize(argc, argv);
    int result = 0;
    try {
        int rank = 0;
        int size = 1;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        int next = (rank + 1) % size;
        int previous = (rank + size - 1) % size;
        CommunicationPlan plan;
        plan.initialize(MPI_COMM_WORLD);
        plan.initialize_graph_communicator(1, &next, 1, &previous);
        CArrayKokkos<size_t> counts(1, "thermal_test_neighbor_counts");
        Kokkos::parallel_for("thermal_test_counts", 1, KOKKOS_LAMBDA(int) {
            counts(0) = 1;
        });
        Kokkos::fence();
        DRaggedRightArrayKokkos<int> sends(counts, "thermal_test_send_indices");
        DRaggedRightArrayKokkos<int> receives(counts, "thermal_test_recv_indices");
        sends.host(0, 0) = 0;
        receives.host(0, 0) = 1;
        sends.update_device();
        receives.update_device();
        plan.setup_send_recv(sends, receives);

        node_t nodes;
        const std::vector<node_state> fields = {node_state::temp, node_state::coords, node_state::velocity};
        nodes.initialize(2, 3, fields, plan);
        if (nodes.temp.size() != 2 || nodes.temp_n0.size() != 2 ||
            nodes.coords.size() != 6 || nodes.vel.size() != 6)
            throw std::runtime_error("Node initialization did not allocate the requested fields");
        const auto temperature = nodes.temp;
        const auto initial = nodes.temp_n0;
        Kokkos::parallel_for("thermal_test_fill", 2, KOKKOS_LAMBDA(int node) {
            temperature(node) = node == 0 ? 300.0 + rank : -1.0;
            initial(node) = node == 0 ? 200.0 + rank : -1.0;
        });
        Kokkos::fence();
        nodes.initialize(2, 3, fields, plan);
        if (nodes.temp.device_pointer() != temperature.device_pointer() ||
            nodes.temp_n0.device_pointer() != initial.device_pointer())
            throw std::runtime_error("Repeated initialization replaced allocated temperature storage");
        nodes.temp.communicate();
        nodes.temp_n0.communicate();
        nodes.temp.update_host();
        nodes.temp_n0.update_host();
        // MATAR deliberately skips exchanges in a one-process communicator.
        const double expected_temperature = size > 1 ? 300.0 + previous : -1.0;
        const double expected_initial = size > 1 ? 200.0 + previous : -1.0;
        if (nodes.temp.host(0) != 300.0 + rank || nodes.temp.host(1) != expected_temperature ||
            nodes.temp_n0.host(0) != 200.0 + rank || nodes.temp_n0.host(1) != expected_initial)
            throw std::runtime_error("Temperature fields did not retain values and exchange the expected ghost node");
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        result = 1;
    }
    int global_result = 0;
    MPI_Allreduce(&result, &global_result, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    Kokkos::finalize();
    MPI_Finalize();
    return global_result;
}
