#ifndef CONVERSIONMESHES_H
#define CONVERSIONMESHES_H

#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <vcg/complex/allocate.h>

namespace vcg{
namespace tri{
    template<class MeshType >
    void Mesh2Matrix(MeshType& convertingMesh, Eigen::MatrixXd &output_vertexes, Eigen::MatrixXi &output_faces){
        typedef typename tri::MeshToMatrix<MeshType>::MatrixXm MatrixXm;
        MatrixXm V_temp;

        //convert mesh in matrix
        vcg::tri::MeshToMatrix< MeshType >::GetTriMeshData(convertingMesh, output_faces, V_temp);
        output_vertexes = V_temp.template cast<double>();

    }

    template<class MeshType >
    void Matrix2Mesh(MeshType& convertedMesh, Eigen::MatrixXd vertexes, Eigen::MatrixXi faces){
        typedef typename MeshType::VertexPointer  VertexPointer;
        typedef typename MeshType::FaceIterator   FaceIterator;

        //reconvert V and F matrixes in a mesh
        convertedMesh.Clear();
        Allocator<MeshType>::AddVertices(convertedMesh,vertexes.rows());
        Allocator<MeshType>::AddFaces(convertedMesh,faces.rows());
        VertexPointer ivp[vertexes.rows()];

        int i;
        for (i=0; i < vertexes.rows(); i++){
            for (int j = 0; j < 3; j++)
                convertedMesh.vert[i].P()[j] = vertexes(i,j);
            ivp[i]=&convertedMesh.vert[i];
        }

        FaceIterator fi;
        for (i=0,fi=convertedMesh.face.begin();fi!=convertedMesh.face.end();i++,fi++)
            for (int j = 0; j < 3; j++)
                (*fi).V(j) = ivp[faces(i, j)];
    }
}}

#endif // CONVERSIONMESHES_H
