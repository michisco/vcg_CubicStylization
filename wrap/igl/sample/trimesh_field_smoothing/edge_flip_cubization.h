#ifndef EDGE_FLIP_CUBIZATION_H
#define EDGE_FLIP_CUBIZATION_H

#include <iostream>
#include <ctime>
#include <vector>
#include <Eigen/Core>

#include <vcg/simplex/face/topology.h>
#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <vcg/complex/allocate.h>
#include <wrap/io_trimesh/export.h>

namespace vcg{
namespace tri{

template<class MeshType >
void Edge_flip_cubization(
        MeshType &out,
        Eigen::MatrixXi &F,
        Eigen::MatrixXd &U,
        Eigen::VectorXd energyVerts)
{
    typedef typename MeshType::FaceIterator   FaceIterator;
    //typedef typename MeshType::FaceType FaceType;
    typedef typename MeshType::CoordType CoordType;

    Matrix2Mesh(out, U, F);

    FaceIterator fi;
    tri::UpdateTopology<MeshType>::FaceFace(out);
    int i;

    //initialize to inf
    double best_energy = std::numeric_limits<double>::infinity();

    for(i=0, fi=out.face.begin();fi!=out.face.end();i++,fi++){
        for(int j = 0; j < 3; j++){
            if(vcg::face::CheckFlipEdge((*fi), (*fi).FFi(j))){
                if( math::ToDeg( Angle(((*fi).FFp(j))->cN(), (*fi).cN()) ) <= 0.1f ){
                    CoordType v0, v1, v2, v3;
                    int i = (*fi).FFi(j);

                    v0 = (*fi).P0(i);
                    v1 = (*fi).P1(i);
                    v2 = (*fi).P2(i);
                    v3 = ((*fi).FFp(i))->P2((*fi).FFi(i));

                    // Take the parallelogram formed by the adjacent faces of edge
                    // If a corner of the parallelogram on extreme of edge to flip is >= 180
                    // the flip produce two identical faces - avoid this
                    if( (Angle(v2 - v0, v1 - v0) + Angle(v3 - v0, v1 - v0) < M_PI) &&
                        (Angle(v2 - v1, v0 - v1) + Angle(v3 - v1, v0 - v1) < M_PI))
                    {
                        // if any of two faces adj to edge in non writable, the flip is unfeasible
                        if((*fi).IsW() && ((*fi).FFp(i))->IsW())
                            vcg::face::FlipEdge((*fi), (*fi).FFi(j));
                    }
                }
            }

            //calcolare energia dopo flip
            //if e1 < best_energy -> best_energy = e1 && flip_edge = (*fi, index_edge)

            //re-flip (annulla flip)
        }
    }

    //eseguire solo flip migliore

    //convert again the mesh in U and F matrixes
    Matrix2Mesh(out, U, F);
}
}}
#endif // EDGE_FLIP_CUBIZATION_H
