#ifndef EDGE_FLIP_CUBIZATION_H
#define EDGE_FLIP_CUBIZATION_H

#include "cube_style_precomputation.h"
#include <iostream>
#include <ctime>
#include <vector>
#include <Eigen/Core>

#include <vcg/simplex/face/topology.h>
#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <vcg/complex/allocate.h>
#include <wrap/io_trimesh/export.h>
#include <energy_computing.h>

#include <conversionMeshes.h>

#ifndef ANGLE_NORM
#define ANGLE_NORM (2.f/3.f)*M_PI //120°
#endif
#ifndef RATIO_THRESHOLD
#define RATIO_THRESHOLD 0.3f
#endif
namespace vcg{
namespace tri{

template<class MeshType >
void Edge_flip_cubization(MeshType &out,
        Eigen::MatrixXi &F,
        Eigen::MatrixXd &U,
        cube_style_data &data)
{
    typedef typename MeshType::FaceIterator FaceIterator;
    typedef typename MeshType::CoordType CoordType;

    Eigen::MatrixXd U_temp;
    Eigen::MatrixXi F_temp;

    Matrix2Mesh(out, U, F);

    FaceIterator fi;
    tri::UpdateTopology<MeshType>::FaceFace(out);

    //initialize to inf
    double best_energy = std::numeric_limits<double>::infinity();
    FaceIterator best_face;
    int best_edge = -1;

    for(fi=out.face.begin();fi!=out.face.end(); fi++){
        for(int j = 0; j < 3; j++){
            if(vcg::face::CheckFlipEdge((*fi), j)){
                if( math::ToDeg( Angle(((*fi).FFp(j))->cN(), (*fi).cN()) ) <= 0.1f ){
                    CoordType v0, v1, v2, v3;
                    int z = (*fi).FFi(j);

                    v0 = (*fi).P0(z);
                    v1 = (*fi).P1(z);
                    v2 = (*fi).P2(z);
                    v3 = ((*fi).FFp(z))->P2((*fi).FFi(z));

                    // Take the parallelogram formed by the adjacent faces of edge
                    // If a corner of the parallelogram on extreme of edge to flip is >= 180
                    // the flip produce two identical faces - avoid this
                    if( (Angle(v2 - v0, v1 - v0) + Angle(v3 - v0, v1 - v0) < M_PI) &&
                        (Angle(v2 - v1, v0 - v1) + Angle(v3 - v1, v0 - v1) < M_PI))
                    {

                        if(vcg::face::CheckFlipEdgeNormal((*fi), j, ANGLE_NORM)){


                            vcg::face::FlipEdge((*fi), j);

                            double meanRatio = vcg::QualityMeanRatio((*fi).P(0), (*fi).P(1), (*fi).P(2));
                            double meanRatio_adj_0 = vcg::QualityMeanRatio(((*fi).FFp(0))->P(0), ((*fi).FFp(0))->P(1), ((*fi).FFp(0))->P(2));
                            double meanRatio_adj_1 = vcg::QualityMeanRatio(((*fi).FFp(1))->P(0), ((*fi).FFp(1))->P(1), ((*fi).FFp(1))->P(2));
                            double meanRatio_adj_2 = vcg::QualityMeanRatio(((*fi).FFp(2))->P(0), ((*fi).FFp(2))->P(1), ((*fi).FFp(2))->P(2));

                            if(meanRatio >= RATIO_THRESHOLD && meanRatio_adj_0 >= RATIO_THRESHOLD && meanRatio_adj_1 >= RATIO_THRESHOLD && meanRatio_adj_2 >= RATIO_THRESHOLD){
                                U_temp.setZero(U.rows(), U.cols());
                                F_temp.setZero(F.rows(), F.cols());
                                Mesh2Matrix(out, U_temp, F_temp);

                                //get vertex indexes
                                int vertexes [4] = { tri::Index(out, (*fi).V0(z)),
                                                     tri::Index(out, (*fi).V1(z)),
                                                     tri::Index(out, (*fi).V2(z)),
                                                     tri::Index(out, ((*fi).FFp(z))->V2((*fi).FFi(z)))};

                                //compute energy after the flip
                                double e1 = energy_computing(U_temp, data, vertexes);

                                if(e1 < best_energy){
                                    best_energy = e1;
                                    best_face = fi;
                                    best_edge = j;
                                }
                            }

                            //remove flip
                            vcg::face::FlipEdge((*fi), (j+1)%3);
                        }                       
                    }
                }
            }
        }
    }

    if(best_edge >= 0){
        //eseguire solo flip migliore
        if(vcg::face::CheckFlipEdge((*best_face), best_edge)){
            if( math::ToDeg( Angle(((*best_face).FFp(best_edge))->cN(), (*best_face).cN()) ) <= 0.1f ){
                CoordType v0, v1, v2, v3;
                int z = (*best_face).FFi(best_edge);

                v0 = (*best_face).P0(z);
                v1 = (*best_face).P1(z);
                v2 = (*best_face).P2(z);
                v3 = ((*best_face).FFp(z))->P2((*best_face).FFi(z));

                if( (Angle(v2 - v0, v1 - v0) + Angle(v3 - v0, v1 - v0) < M_PI) &&
                    (Angle(v2 - v1, v0 - v1) + Angle(v3 - v1, v0 - v1) < M_PI)){

                    if(vcg::face::CheckFlipEdgeNormal((*best_face), best_edge, ANGLE_NORM)){
                        vcg::face::FlipEdge((*best_face), best_edge);

                        double meanRatio = vcg::QualityMeanRatio((*best_face).P(0), (*best_face).P(1), (*best_face).P(2));
                        double meanRatio_adj_0 = vcg::QualityMeanRatio(((*best_face).FFp(0))->P(0), ((*best_face).FFp(0))->P(1), ((*best_face).FFp(0))->P(2));
                        double meanRatio_adj_1 = vcg::QualityMeanRatio(((*best_face).FFp(1))->P(0), ((*best_face).FFp(1))->P(1), ((*best_face).FFp(1))->P(2));
                        double meanRatio_adj_2 = vcg::QualityMeanRatio(((*best_face).FFp(2))->P(0), ((*best_face).FFp(2))->P(1), ((*best_face).FFp(2))->P(2));

                        if(meanRatio >= RATIO_THRESHOLD && meanRatio_adj_0 >= RATIO_THRESHOLD && meanRatio_adj_1 >= RATIO_THRESHOLD && meanRatio_adj_2 >= RATIO_THRESHOLD)
                            std::cout<<"FLIPPED " +  std::to_string(best_energy)<<std::endl;
                        else
                            vcg::face::FlipEdge((*best_face), (best_edge+1)%3);
                    }
                }
            }
        }
        //convert again the mesh in U and F matrixes
        Mesh2Matrix(out, U, F);
    }
}
}}
#endif // EDGE_FLIP_CUBIZATION_H
