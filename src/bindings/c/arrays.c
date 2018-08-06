
#include "arrays.h"

void set_node_field_value_c(int *row, int *col, double *value);
void check_node_xyz_2d_c(int *row, int *col, double *value);
void check_elem_2d_version_c(int *row, int *col, double *value);
void check_node_dxyz_2d_c(int *nderiv, int *row, int *col, double *value);
void check_nodes_in_elem_2d_c(int *row, int *col, int *value);

void set_node_field_value(int row, int col, double value)
{
  set_node_field_value_c(&row, &col, &value);
}

double check_node_xyz_2d(int row, int col)
{
  double value;
  check_node_xyz_2d_c(&row, &col, &value);
  return value;
}

double check_elem_2d_version(int row, int col)
{
  double value;
  check_elem_2d_version_c(&row, &col, &value);
  return value;
}

double check_node_dxyz_2d(int nderiv, int row, int col)
{
  double value;
  check_node_dxyz_2d_c(&nderiv, &row, &col, &value);
  return value;
}

int check_nodes_in_elem_2d(int row, int col)
{
  int value;
  check_nodes_in_elem_2d_c(&row, &col, &value);
  return value;
}
