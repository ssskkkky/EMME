#ifndef GRID_H  // Replace MATRIX_H with your unique guard macro name
#define GRID_H

#include <vector>

template <typename T>
struct Grid {
    Grid(T leni, int npointsi)
        : len(leni),
          npoints(npointsi),
          dx((2 * leni) / (npoints - 1)),
          grid(npoints) {
        for (int i; i < npoints; i++) { grid[i] = -leni + i * dx; }
    }

    T len;
    int npoints;
    T dx;
    std::vector<T> grid;
};

#endif
