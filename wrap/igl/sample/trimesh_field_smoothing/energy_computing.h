#ifndef ENERGY_COMPUTING_H
#define ENERGY_COMPUTING_H

#include <iostream>
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <cube_style_data.h>
#include <orthogonal_procrustes.h>
#include <shrinkage.h>
#include <igl/slice.h>
#include <igl/parallel_for.h>
#include <math.h>

double energy_computing(
    Eigen::MatrixXd U,
    Eigen::MatrixXd RAll,
    cube_style_data &data,
    int verts[]);
#endif // ENERGY_COMPUTING_H
