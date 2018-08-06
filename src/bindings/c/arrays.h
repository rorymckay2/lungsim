
#ifndef AETHER_ARRAYS_H
#define AETHER_ARRAYS_H

#include "symbol_export.h"

SHO_PUBLIC void set_node_field_value(int row, int col, double value);
SHO_PUBLIC double check_node_xyz_2d(int row, int col);
SHO_PUBLIC double check_elem_2d_version(int row, int col);
SHO_PUBLIC double check_node_dxyz_2d(int nderiv, int row, int col);
SHO_PUBLIC int check_nodes_in_elem_2d(int row, int col);

#endif /* AETHER_ARRAYS_H */
