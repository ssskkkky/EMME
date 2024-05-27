#ifndef GRID_H  // Replace MATRIX_H with your unique guard macro name
#define GRID_H

#include <vector>

template <typename T>
struct Grid {
    Grid(T leni, unsigned int npointsi)
        : len(leni),
          npoints(npointsi),
          dx((2 * leni) / (npoints - 1)),
          grid(npoints) {
        for (unsigned int i = 0; i < npoints; i++) { grid[i] = -len + i * dx; }
    }

    T len;
    unsigned int npoints;
    T dx;
    std::vector<T> grid;
};

#endif
