#ifndef __AF_BOND_LIST_H__
#define __AF_BOND_LIST_H__

#include <vector>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include "../types.h"

#define MAX_BONDS 16

using namespace std;

namespace af
{
    typedef thrust::host_vector<int> BondListHost; 
    typedef thrust::device_vector<int> BondListDevice;

    class BondListBuilder
    {
        bool base_one;
        vector<vector<uint>> bonds;

    public:
        BondListBuilder(bool base_one_indices = false) 
        : base_one(base_one_indices){}


        void add_bond(uint particleA_id, uint particleB_id)
        {
            // if base one, shift to base zero
            if (base_one)
            {
                particleA_id--;
                particleB_id--;
            }

            // figure out what's the biggest idx
            uint max_id = (particleA_id > particleB_id ? particleA_id : particleB_id);

            // extend the bond array if needed
            size_t size = bonds.size();
            if (bonds.size() <= max_id)
                for (int i = 0; i <= max_id - size; i++)
                {
                    bonds.push_back(vector<uint>());
                }

            // get number of current bonds
            size_t sizeA = bonds.at(particleA_id).size();
            size_t sizeB = bonds.at(particleB_id).size();

            // add bonds only if both can take one more
            if (sizeA <= MAX_BONDS && sizeB <= MAX_BONDS)
            {
                bonds.at(particleA_id).push_back(particleB_id);
                bonds.at(particleB_id).push_back(particleA_id);
                std::cout << "Added bond " << particleA_id << "-" << particleB_id << std::endl;
            }
        }

        BondListHost build_for_host()
        {
            cout << "Building bond list for host..." << endl;
            auto list = build_list();
            return BondListHost(list.begin(), list.end()); // build host then copy to device
        }

        BondListDevice build_for_device()
        {
            cout << "Building bond list for device..." << endl;
            auto list = build_list();
            return BondListDevice(list.begin(), list.end()); // build host then copy to device
        }

    private:

        size_t bond_list_size()
        {
            return bonds.size() * MAX_BONDS;
        }

        vector<int> build_list()
        {
            size_t size = bond_list_size();
            size_t nparts = bonds.size(); // number of particles
            vector<int> list(size, -1); // create flattened array with -1

            for (int i = 0; i < nparts; i++)
            {
                size_t p_bonds = bonds.at(i).size(); // bonds ith particle has
                for (int j = 0; j < p_bonds; j++)
                {
                    list[j + i*MAX_BONDS] = bonds[i][j]; // ith particle bonding to its jth bond
                }
            }

            return list;
        }
    };
}

#endif