#ifndef CUBIC_STYLIZING_H
#define CUBIC_STYLIZING_H

#include <cube_style_data.h>
#include <cube_style_precomputation.h>
#include <cube_style_single_iteration.h>
#include <normalize_unitbox.h>
#include <exporter_cubic.h>

#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <vcg/complex/allocate.h>

namespace vcg{
namespace tri{

template<class MeshType >
void Stylize_Cubic( MeshType& m, MeshType&o, double cubeness, Eigen::VectorXd &energyVertexes, std::string outName)
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
    Eigen:: MatrixXd energyVects;

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
        cube_style_single_iteration(V,U,data,energyVertexes);
        std::cout << "iteration: " << iter << " Total energy: " << std::scientific <<  data.objVal << std::endl;
        if (data.reldV < stopReldV) break;
        if(iter%30 == 0){

            //reconvert V and F matrixes in a mesh
            o.Clear();
            Allocator<MeshType>::AddVertices(o,U.rows());
            Allocator<MeshType>::AddFaces(o,F.rows());
            VertexPointer ivp[U.rows()];

            int i;
            for (i=0; i < U.rows(); i++){
                for (int j = 0; j < 3; j++)
                    o.vert[i].P()[j] = U(i,j);
                ivp[i]=&o.vert[i];
            }

            FaceIterator fi;
            for (i=0,fi=o.face.begin();fi!=o.face.end();i++,fi++)
                for (int j = 0; j < 3; j++)
                    (*fi).V(j) = ivp[F(i, j)];

            std::string outfile = outName + "_qualityVertex_" + std::to_string(iter);
            exporter_cubic_colorize(o, energyVertexes, outfile);
        }
    }

    std::cout << "Total energy: " << std::scientific <<  data.objVal << std::endl;

    //reconvert V and F matrixes in a mesh
    o.Clear();
    Allocator<MeshType>::AddVertices(o,U.rows());
    Allocator<MeshType>::AddFaces(o,F.rows());
    VertexPointer ivp[U.rows()];

    int i;
    for (i=0; i < U.rows(); i++){
        for (int j = 0; j < 3; j++)
            o.vert[i].P()[j] = U(i,j);
        ivp[i]=&o.vert[i];
    }

    FaceIterator fi;
    for (i=0,fi=o.face.begin();fi!=o.face.end();i++,fi++)
        for (int j = 0; j < 3; j++)
            (*fi).V(j) = ivp[F(i, j)];
}
}}
#endif // CUBIC_STYLIZING_H
