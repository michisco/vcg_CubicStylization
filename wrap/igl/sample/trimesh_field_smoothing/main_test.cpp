#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/snap_points.h>

#include <cube_style_data.h>
#include <cube_style_precomputation.h>
#include <cube_style_single_iteration.h>
#include <normalize_unitbox.h>

#include <vcg/complex/algorithms/mesh_to_matrix.h>

#include <wrap/io_trimesh/import_off.h>
#include <wrap/io_trimesh/export_off.h>

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

#ifndef OUTPUT_PATH
#define OUTPUT_PATH "../"
#endif

// to run the code, type "./cubeStyle_bin [meshName] [lambda]"
int main(int argc, char *argv[])
{
    // load mesh and lambda
    Eigen::MatrixXd V, U;
    MatrixXi F;
    cube_style_data data;
    string meshName;
    {
        if (argc == 1)
        {
            meshName = "spot.obj"; // default mesh
            data.lambda = 0.2; // default lambda
        }
        else if (argc == 2)
        {
            meshName = argv[1];
            data.lambda = 0.2; // default lambda
        }
        else
        {
            meshName = argv[1];
            data.lambda = stod(argv[2]);
        }
        string file = MESH_PATH + meshName;

        int n = file.length();
        // declaring character array
        char char_array[n + 1];

        // copying the contents of the
        // string to char array
        strcpy(char_array, file.c_str());

        MyMesh mesh_obj;
        tri::io::ImporterOFF<MyMesh>::Open(mesh_obj,char_array);
        tri::MeshToMatrix<MyMesh>::MatrixXm V_temp;

        tri::MeshToMatrix<MyMesh>::GetTriMeshData(mesh_obj, F, V_temp);
        V = V_temp.cast<double>();

        //igl::readOBJ(file, V, F);
        normalize_unitbox(V);
        RowVector3d meanV = V.colwise().mean();
        V = V.rowwise() - meanV;
        U = V;
    }

    // set a constrained point F(0,0)
    {
        data.bc.resize(1,3);
        data.bc << V.row(F(0,0));

        data.b.resize(1);
        data.b << F(0,0);
    }

    // precomputation ARAP and initialize ADMM parameters
    cube_style_precomputation(V,F,data);

    // cubic stylization
    int maxIter = 1000;
    double stopReldV = 1e-3; // stopping criteria for relative displacement
    for (int iter=0; iter<maxIter; iter++)
    {
        cout << "iteration: " << iter << endl;
        cube_style_single_iteration(V,U,data);
        if (data.reldV < stopReldV) break;
    }

    // write output mesh
    {
        string outputName = "cubic_";
        outputName.append("bunny_test.obj");
        string outputFile = OUTPUT_PATH + outputName;
        igl::writeOBJ(outputFile,U,F);
        //tri::io::ExporterOFF<MyMesh>::Save(meshName,"Test.off",tri::io::Mask::IOM_FACECOLOR);

        /*string inputName = "input_";
        inputName.append(meshName);
        string inputFile = OUTPUT_PATH + inputName;
        igl::writeOBJ(inputFile,V,F);*/
    }

    return 0;
}
