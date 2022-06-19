#ifndef EXPORTER_CUBIC_H
#define EXPORTER_CUBIC_H

#include <iostream>
#include <ctime>
#include <vector>
#include <Eigen/Core>

#include <vcg/complex/allocate.h>
#include <vcg/complex/algorithms/update/color.h>
#include <wrap/io_trimesh/export.h>

namespace vcg{
namespace tri{

template<class MeshType >
void exporter_cubic(
        MeshType &o,
        std::string outName)
{
    using namespace Eigen;
    using namespace std;
    using namespace vcg;
    using namespace tri;

    string outputFile = outName + ".obj";
    const char * outFile = outputFile.c_str();
    tri::io::Exporter<MeshType>::Save(o,outFile);
}

template<class MeshType >
void exporter_cubic_colorize(
        MeshType &o,
        Eigen::VectorXd &energyVertexes,
        std::string outName)
{
    using namespace Eigen;
    using namespace std;
    using namespace vcg;
    using namespace tri;

    assert(energyVertexes.size() == o.vert.size());
    for(int i = 0; i < energyVertexes.size(); i++){
        o.vert[i].Q() = energyVertexes[i];
    }

    //colorize quality vertexes
    UpdateColor<MeshType>::PerVertexQualityRamp(o, 0.00000418, 0.000644);

    string outputFileAlt = outName + ".ply";
    const char * outFileAlt = outputFileAlt.c_str();
    io::ExporterPLY<MeshType>::Save(o, outFileAlt, io::Mask::IOM_VERTQUALITY+tri::io::Mask::IOM_VERTCOLOR);
}

}}
#endif // EXPORTER_CUBIC_H
