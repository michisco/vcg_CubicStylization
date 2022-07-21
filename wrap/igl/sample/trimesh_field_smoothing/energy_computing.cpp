#include "energy_computing.h"

double energy_computing(
    Eigen::MatrixXd U,
    cube_style_data &data,
    int verts[]){

    using namespace Eigen;
    using namespace std;

    Eigen::VectorXd energyVec;
    energyVec.setZero(U.rows());

    igl::parallel_for(
        4,
        [&data, &U, &energyVec, &verts](const int ii)
        {
            // warm start parameters
            VectorXd n = data.N.row(verts[ii]).transpose();
            VectorXd z = data.zAll.col(verts[ii]);
            VectorXd u = data.uAll.col(verts[ii]);
            double rho = data.rhoAll(verts[ii]);
            Matrix3d R; //Matrix3d R = RAll.block(0,3*verts[ii],3,3);

            // get energy parameters
            // Note: dVn = [dV n], dUn = [dU z-u]
            MatrixXi hE = data.hEList[verts[ii]];
            MatrixXd dU(3,hE.rows());
            {
                MatrixXd U_hE0, U_hE1;
                igl::slice(U,hE.col(0),1,U_hE0);
                igl::slice(U,hE.col(1),1,U_hE1);
                dU = (U_hE1 - U_hE0).transpose();
            }

            MatrixXd dV = data.dVList[verts[ii]];
            VectorXd WVec = data.WVecList[verts[ii]];
            Matrix3d Spre = dV * WVec.asDiagonal() * dU.transpose();

            // R step
            Matrix3d S = Spre + (rho * n * (z-u).transpose());
            orthogonal_procrustes(S, R);

            // save objective
            double objVal =
                0.5*((R*dV-dU)*WVec.asDiagonal()*(R*dV-dU).transpose()).trace()
               + data.lambda * data.VA(verts[ii]) * (R*n).cwiseAbs().sum();
            energyVec(verts[ii]) = objVal;
        }
    ,1000);

    return energyVec.sum();
}
