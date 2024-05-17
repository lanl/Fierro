#pragma once

#include "matar.h"

using namespace mtr;

template<typename T>
void matar_print(const T &array)
{
    int ord = array.order();

    if (ord == 1)
    {
        int n = array.dims(0);
        for (int i = 0; i < n; i++)
        {
            std::cout << array(i) << " ";
        }
        std::cout << std::endl;
    } else if (ord == 2)
    {
        int m = array.dims(0);
        int n = array.dims(1);
        for (int i = 0; i < m; i++)
        {
            std::cout << i << " ---> ";
            for (int j = 0; j < n; j++)
            {
                std::cout << array(i, j) << " ";
            }
            std::cout << std::endl;
        }
    } else
    {
        int l = array.dims(0);
        int m = array.dims(1);
        int n = array.dims(2);
        for (int i = 0; i < l; i++)
        {
            std::cout << "Level " << i << std::endl;
            for (int j = 0; j < m; j++)
            {
                std::cout << j << " ---> ";
                for (int k = 0; k < n; k++)
                {
                    std::cout << array(i, j, k) << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
}
