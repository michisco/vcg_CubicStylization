#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/snap_points.h>

#include <cube_style_data.h>
#include <cube_style_precomputation.h>
#include <cube_style_single_iteration.h>
#include <normalize_unitbox.h>
#include <cubic_stylizing.h>

#include <vcg/complex/algorithms/mesh_to_matrix.h>
#include <wrap/io_trimesh/import.h>
#include <wrap/io_trimesh/export.h>

#include <cstring>

#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>

class MyVertex; class MyEdge; class MyFace;
struct MyUsedTypes : public vcg::UsedTypes<vcg::Use<MyVertex>   ::AsVertexType,
                                           vcg::Use<MyEdge>     ::AsEdgeType,
                                           vcg::Use<MyFace>     ::AsFaceType>{};

class MyVertex  : public vcg::Vertex< MyUsedTypes, vcg::vertex::Coord3f, vcg::vertex::Normal3f, vcg::vertex::VFAdj, vcg::vertex::BitFlags  >{};
class MyFace    : public vcg::Face<   MyUsedTypes, vcg::face::FFAdj, vcg::face::VFAdj, vcg::face::VertexRef, vcg::face::Color4b, vcg::face::BitFlags > {};
class MyEdge    : public vcg::Edge<   MyUsedTypes> {};
class MyMesh    : public vcg::tri::TriMesh< std::vector<MyVertex>, std::vector<MyFace> , std::vector<MyEdge>  > {};

using namespace vcg;

using namespace Eigen;
using namespace std;

#ifndef MESH_PATH
#define MESH_PATH "../../../../meshes/"
#endif

// to run the code, type "./trimesh_field_smoothing [meshName] [lambda] [outputMeshName]"
int main(int argc, char *argv[])
{
    // load mesh and lambda
    MatrixXd V, U;
    MatrixXi F;
    double lambda;
    string meshName;
    string outputName;

    if (argc == 1)
    {
        meshName = "spot.obj"; // default mesh
        lambda = 0.2; // default lambda
        outputName = "cubic_spot.obj"; //default output mesh
    }
    else if (argc == 2)
    {
        meshName = argv[1];
        lambda = 0.2; // default lambda
        outputName = "cubic_spot.obj"; //default output mesh
    }
    else if(argc == 3){
        meshName = argv[1];
        lambda = stod(argv[2]);
        outputName = "cubic_spot.obj"; //default output mesh
    }
    else
    {
        meshName = argv[1];
        lambda = stod(argv[2]);
        outputName = argv[3];
    }

    string file = MESH_PATH + meshName;

    int n = file.length();
    // declaring character array
    char char_array[n + 1];

    // copying the contents of the
    // string to char array
    strcpy(char_array, file.c_str());

    MyMesh mesh_obj;
    MyMesh output_mesh;
    tri::io::Importer<MyMesh>::Open(mesh_obj,char_array);

    tri::Stylize_Cubic(mesh_obj, output_mesh, lambda);

    // write output mesh
    string outputFile = MESH_PATH + outputName;
    const char * outFile = outputFile.c_str();
    tri::io::Exporter<MyMesh>::Save(output_mesh,outFile);

    return 0;
}
