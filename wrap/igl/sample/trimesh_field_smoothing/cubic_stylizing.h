#ifndef CUBIC_STYLIZING_H
#define CUBIC_STYLIZING_H

#include <cube_style_data.h>
#include <cube_style_precomputation.h>
#include <cube_style_single_iteration.h>
#include <normalize_unitbox.h>
#include <exporter_cubic.h>
#include <edge_flip_cubization.h>

#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <vcg/complex/allocate.h>

#include <conversionMeshes.h>

namespace vcg{
namespace tri{

template<class MeshType >
void Stylize_Cubic( MeshType& m, MeshType&o, double cubeness, Eigen::VectorXd &energyVertexes, std::string outName, int isFlip)
{
    // check requirements
    vcg::tri::VertexVectorHasPerVertexTexCoord( m.vert );
    vcg::tri::VertexVectorHasPerVertexFlags( m.vert );

    cube_style_data data;
    data.lambda = cubeness;

    Eigen::MatrixXd V;
    Eigen::MatrixXd U;
    Eigen::MatrixXi F;
    Eigen::MatrixXd V_uv;
    Eigen::VectorXd prev_energyVertexes;

    //convert mesh in matrix
    Mesh2Matrix(m, V, F);

    normalize_unitbox(V);
    Eigen::RowVector3d meanV = V.colwise().mean();
    V = V.rowwise() - meanV;
    U = V;

    data.bc.resize(1,3);
    data.bc << V.row(F(0,0));

    data.b.resize(1);
    data.b << F(0,0);

    // precomputation ARAP and initialize ADMM parameters
    cube_style_precomputation(V,F,data);

    // apply cubic stylization
    int maxIter = 1000;
    double stopReldV = 1e-3; // stopping criteria for relative displacement
    for (int iter=0; iter<maxIter; iter++)
    {
        //prev_energyVertexes = energyVertexes;
        Eigen::MatrixXd R;
        cube_style_single_iteration(V,U,data,energyVertexes,R);

        //apply edge flips
        if(isFlip)
            Edge_flip_cubization(o,F,U,data,R);

        std::cout << "iteration: " << iter << " Total energy: " << std::scientific <<  data.objVal << std::endl;
        if (data.reldV < stopReldV) break;
        if(iter%30 == 0){
            //convert matrixes in mesh
            Matrix2Mesh(o, U, F);

            std::string outfile = outName + "_qualityVertex_" + std::to_string(iter);
            exporter_cubic_colorize(o, energyVertexes, outfile);
        }
    }

    std::cout << "Total energy: " << std::scientific <<  data.objVal << std::endl;

    //convert matrixes in mesh
    Matrix2Mesh(o, U, F);
}
}}
#endif // CUBIC_STYLIZING_H
