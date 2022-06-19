#ifndef COMPUTING_ENERGY_H
#define COMPUTING_ENERGY_H

#include <cube_style_data.h>
#include <cube_style_precomputation.h>
#include <cube_style_single_iteration.h>
#include <normalize_unitbox.h>

#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <vcg/complex/allocate.h>

namespace vcg{
namespace tri{

template<class MeshType >
void Stylize_Cubic( MeshType& m, double cubeness)
{

    typedef typename MeshType::VertexPointer  VertexPointer;
    typedef typename MeshType::FaceIterator   FaceIterator;
    typedef typename tri::MeshToMatrix<MeshType>::MatrixXm MatrixXm;

    // check requirements
    vcg::tri::VertexVectorHasPerVertexTexCoord( m.vert );
    vcg::tri::VertexVectorHasPerVertexFlags( m.vert );

    cube_style_data data;
    data.lambda = cubeness;

    Eigen::MatrixXd V;
    Eigen::MatrixXd U;
    Eigen::MatrixXi F;
    Eigen::MatrixXd V_uv;
    MatrixXm V_temp;

    //convert mesh in matrix
    vcg::tri::MeshToMatrix< MeshType >::GetTriMeshData(m, F, V_temp);
    V = V_temp.template cast<double>();

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
        std::cout << "iteration: " << iter << " reldV: " << std::scientific << data.reldV << " energy: " << std::scientific <<  data.objVal << std::endl;
        cube_style_single_iteration(V,U,data);
        if (data.reldV < stopReldV) break;
    }

    std::cout << "reldV: " << std::scientific << data.reldV << " energy: " << std::scientific <<  data.objVal << std::endl;
}
}}
#endif // COMPUTING_ENERGY_H
