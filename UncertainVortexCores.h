#ifndef __UNCERTAINVORTEXCORES_VTKVORTEXCORES__
#define __UNCERTAINVORTEXCORES_VTKVORTEXCORES__

#include "vtkImageAlgorithm.h"
#include "vtkMultiTimeStepAlgorithm.h"
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointData.h"
#include "vtkDoubleArray.h"
#include "vtkMultiBlockDataSet.h"
#include "chrono"
#include "vector"
#include <set>
#include <random>

//Override Eigen standard allocation limit
#ifndef EIGEN_STACK_ALLOCATION_LIMIT
#define EIGEN_STACK_ALLOCATION_LIMIT 1000000
#endif
#include <Eigen/Eigen>
#include <Eigen/unsupported/MatrixFunctions>
#include <omp.h>
#include "math.h"
#include "linalg.h"

//#include "VCGJacobi.h"
//#include "parallelVectors.h"
//#include "vtkFloatArray.h"
//#include "vtkDataArray.h"
//#include "vtkDataSet.h"
//#include "vtkRectilinearGrid.h"
//#include "vtkPolyData.h"
//#include "exception"
//#include "vtkObjectFactory.h"

typedef Eigen::Matrix<double,3,1> Vector3d;
typedef Eigen::Matrix<double,3,3> Matrix3d;
typedef Eigen::Matrix<double,96,1> Vector96d;
typedef Eigen::Matrix<std::complex<double>,96,1> Vector96c;
typedef Eigen::Matrix<double,96,96> Matrix96d;
typedef Eigen::Matrix<std::complex<double>,96,96> Matrix96c;
typedef std::chrono::high_resolution_clock nanoClock;

class UncertainVortexCores : public vtkMultiTimeStepAlgorithm {
public:
    static UncertainVortexCores *New();
    vtkTypeMacro(UncertainVortexCores, vtkMultiTimeStepAlgorithm);
    void PrintSelf(ostream &os, vtkIndent indent) {};

    vtkSetMacro(numSamples, int);
    vtkGetMacro(numSamples, int);

    vtkSetMacro(useRandomSeed, bool);
    vtkGetMacro(useRandomSeed, bool);

    vtkSetMacro(useCholesky, bool);
    vtkGetMacro(useCholesky, bool);

protected:
    UncertainVortexCores();
    ~UncertainVortexCores() {};

    // Make sure the pipeline knows what type we expect as input
    int FillInputPortInformation( int port, vtkInformation* info );
    int FillOutputPortInformation( int port, vtkInformation* info );
    int RequestInformation(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector); //the function that makes this class work with the vtk pipeline
    // Generate output
    int RequestData(vtkInformation *vtkNotUsed(request), vtkInformationVector **inputVector, vtkInformationVector *outputVector);
    int RequestUpdateExtent(vtkInformation*,vtkInformationVector** inputVector,vtkInformationVector* outputVector);

private:
    vtkSmartPointer<vtkImageData> data;
    vtkDoubleArray *cellValues;

    int numFields;
    bool NewtonQuad;
    int numSamples;
    double *spacing;
    int *gridResolution;
    int arrayLength;
    int offsetY;
    int offsetZ;
    nanoClock::time_point beginning;
    char *cellValuesName;
    bool useRandomSeed;
    bool useCholesky;
    std::mt19937 gen;

    Vector96d generateNormalDistributedVec();
    bool computeParVectorsInCell(std::vector<std::tuple<Vector3d, Matrix3d>> cell);
    std::vector<Vector3d> calculateJacobiTimesVelocity(std::vector<std::tuple<Vector3d, Matrix3d>> cell);
    std::vector<std::tuple<Vector3d, Matrix3d>> calculateJacobiAndCellNodes(Vector96d sampleVector);
    int computeParallelOnCellface(std::vector<Vector3d> cellface, std::vector<Vector3d> cellfaceAcc, double *s, double *t);
    bool interpolatedJacobiIsComplex(Matrix3d jac1, Matrix3d jac2, Matrix3d jac3, double s, double t);
    bool computeParallelVectors(Vector96d sampleVector);
    bool isCloseToEdge(int index);
};

#endif