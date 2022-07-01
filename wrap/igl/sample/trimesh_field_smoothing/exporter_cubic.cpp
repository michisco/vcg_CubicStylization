#include <exporter_cubic.h>

template<class MeshOutType >
void exporter_cubic(
        MeshOutType &o,
        std::string outName)
{
    using namespace Eigen;
    using namespace std;
    using namespace vcg;
    using namespace tri;

    string outputFile = outName + ".obj";
    const char * outFile = outputFile.c_str();
    tri::io::Exporter<MeshOutType>::Save(o,outFile);
}

template<class MeshColorType >
void exporter_cubic_colorize(
        MeshColorType &o,
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

    string outputFileAlt = outName + ".ply";
    const char * outFileAlt = outputFileAlt.c_str();
    io::ExporterPLY<MeshColorType>::Save(o, outFileAlt, io::Mask::IOM_VERTQUALITY);
}

