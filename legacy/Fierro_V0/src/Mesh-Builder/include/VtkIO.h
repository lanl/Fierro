#include <vector>

namespace {
    /**
     * \brief Given (i,j,k) coordinates within the Lagrange hex, return an offset into the local connectivity (PointIds) array.
     *
     * The \a order parameter must point to an array of 3 integers specifying the order
     * along each axis of the hexahedron.
     */
    inline int ijk_to_lagrange_hex(int i, int j, int k, const int* order) {
        bool ibdy = (i == 0 || i == order[0]);
        bool jbdy = (j == 0 || j == order[1]);
        bool kbdy = (k == 0 || k == order[2]);
        // How many boundaries do we lie on at once?
        int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0) + (kbdy ? 1 : 0);

        if (nbdy == 3) // Vertex DOF
        {              // ijk is a corner node. Return the proper index (somewhere in [0,7]):
            return (i ? (j ? 2 : 1) : (j ? 3 : 0)) + (k ? 4 : 0);
        }

        int offset = 8;
        if (nbdy == 2) // Edge DOF
        {
            if (!ibdy)
            { // On i axis
            return (i - 1) + (j ? order[0] + order[1] - 2 : 0) + (k ? 2 * (order[0] + order[1] - 2) : 0) +
                offset;
            }
            if (!jbdy)
            { // On j axis
            return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) +
                (k ? 2 * (order[0] + order[1] - 2) : 0) + offset;
            }
            // !kbdy, On k axis
            offset += 4 * (order[0] - 1) + 4 * (order[1] - 1);
            return (k - 1) + (order[2] - 1) * (i ? (j ? 3 : 1) : (j ? 2 : 0)) + offset;
        }

        offset += 4 * (order[0] + order[1] + order[2] - 3);
        if (nbdy == 1) // Face DOF
        {
            if (ibdy) // On i-normal face
            {
            return (j - 1) + ((order[1] - 1) * (k - 1)) + (i ? (order[1] - 1) * (order[2] - 1) : 0) +
                offset;
            }
            offset += 2 * (order[1] - 1) * (order[2] - 1);
            if (jbdy) // On j-normal face
            {
            return (i - 1) + ((order[0] - 1) * (k - 1)) + (j ? (order[2] - 1) * (order[0] - 1) : 0) +
                offset;
            }
            offset += 2 * (order[2] - 1) * (order[0] - 1);
            // kbdy, On k-normal face
            return (i - 1) + ((order[0] - 1) * (j - 1)) + (k ? (order[0] - 1) * (order[1] - 1) : 0) +
            offset;
        }

        // nbdy == 0: Body DOF
        offset += 2 * (
            (order[1] - 1) * (order[2] - 1) + 
            (order[2] - 1) * (order[0] - 1) +
            (order[0] - 1) * (order[1] - 1)
        );
        return offset + 
            (i - 1) + (order[0] - 1) * (
                (j - 1) + (order[1] - 1) * (
                    (k - 1)));
    }

    inline int ijk_to_lagrange_quad(int i, int j, const int* order) {
        bool ibdy = (i == 0 || i == order[0]);
        bool jbdy = (j == 0 || j == order[1]);
        // How many boundaries do we lie on at once?
        int nbdy = (ibdy ? 1 : 0) + (jbdy ? 1 : 0);

        if (nbdy == 2) // Vertex DOF
        {              // ijk is a corner node. Return the proper index (somewhere in [0,7]):
            return (i ? (j ? 2 : 1) : (j ? 3 : 0));
        }

        int offset = 4;
        if (nbdy == 1) // Edge DOF
        {
            if (!ibdy)
            { // On i axis
            return (i - 1) + (j ? order[0] - 1 + order[1] - 1 : 0) + offset;
            }
            if (!jbdy)
            { // On j axis
            return (j - 1) + (i ? order[0] - 1 : 2 * (order[0] - 1) + order[1] - 1) + offset;
            }
        }

        offset += 2 * (order[0] - 1 + order[1] - 1);
        // nbdy == 0: Face DOF
        return offset + (i - 1) + (order[0] - 1) * ((j - 1));
    }

    inline std::vector<int> ijk_to_lagrange_hex(const int p_order) {
        static std::vector<int> map;

        int order[3] = {p_order, p_order, p_order};
        for (int k = 0; k <= p_order; k++)
            for (int j = 0; j <= p_order; j++)
                for (int i = 0; i <= p_order; i++)
                    map.push_back(ijk_to_lagrange_hex(i, j, k, order));
        
        return map;
    }

    inline std::vector<int> lagrange_hex_to_ijk(const int p_order) {
        return MeshIO::_Impl::construct_inverse(ijk_to_lagrange_hex(p_order));
    }

    inline std::vector<int> ijk_to_lagrange_quad(const int p_order) {
        std::vector<int> map;

        int order[3] = {p_order, p_order};
        for (int j = 0; j <= p_order; j++)
            for (int i = 0; i <= p_order; i++)
                map.push_back(ijk_to_lagrange_quad(i, j, order));
        
        return map;
    }

    inline std::vector<int> lagrange_quad_to_ijk(const int p_order) {
        return MeshIO::_Impl::construct_inverse(ijk_to_lagrange_quad(p_order));
    }
}

namespace VtkIO {
    inline bool is_lagrange_ordered(int vtkCellType) {
        return (vtkCellType >= 68) && (vtkCellType <= 74);
    }

    inline std::vector<int> lagrange_to_ijk(const int num_dim, const int p_order) {
        switch (num_dim) {
        case 2:
            return lagrange_quad_to_ijk(p_order);
        case 3:
            return lagrange_hex_to_ijk(p_order);
        default:
            throw std::runtime_error("Invalid number of dimensions: " + num_dim);
        }
    }

    inline std::vector<int> ijk_to_lagrange(const int num_dim, const int p_order) {
        return MeshIO::_Impl::construct_inverse(lagrange_to_ijk(num_dim, p_order));
    }
}
