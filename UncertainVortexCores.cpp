#include "UncertainVortexCores.h"

vtkStandardNewMacro(UncertainVortexCores);

UncertainVortexCores::UncertainVortexCores() {

    this->cellValuesName = (char *) "Vortex Core Probability";
    this->SetNumberOfInputPorts(1);
    this->SetNumberOfOutputPorts(1);
    this->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS,
                                 vtkDataSetAttributes::VECTORS);
}

int UncertainVortexCores::FillInputPortInformation(int port, vtkInformation *info) {
    
    if (port == 0) {
        info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataObject");
        return 1;
    }
    return 0;
}

int UncertainVortexCores::FillOutputPortInformation(int port, vtkInformation *info) { //add more outputs if needed
    
    if (port == 0) {
        info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkImageData");
        return 1;
    }
    return 0;
}

int UncertainVortexCores::RequestUpdateExtent(vtkInformation *, vtkInformationVector **inputVector, vtkInformationVector *outputVector) {
    
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);

    if (outInfo->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP())) {
        double upTime = outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP());
        double *inTimes = inInfo->Get(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        int numInTimes = inInfo->Length(vtkStreamingDemandDrivenPipeline::TIME_STEPS());
        inInfo->Set(vtkMultiTimeStepAlgorithm::UPDATE_TIME_STEPS(), inTimes, numInTimes);
    }
    
    return 1;
}

int UncertainVortexCores::RequestInformation(vtkInformation *, vtkInformationVector **inputVector,
                                                 vtkInformationVector *outputVector) {
    return 1;
}


int UncertainVortexCores::RequestData(vtkInformation *, vtkInformationVector **inputVector,
                                          vtkInformationVector *outputVector) {
    
    vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
    //moved here from RequestUpdateExtent, maybe move somewhere else 
    //inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_EXTENT(), extent, 6);
    vtkSmartPointer<vtkMultiBlockDataSet> input = vtkMultiBlockDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
    this->numFields = input->GetNumberOfBlocks();

    this->data = vtkSmartPointer<vtkImageData>::New();
    data->ShallowCopy(vtkImageData::SafeDownCast(input->GetBlock(0)));

    this->spacing = data->GetSpacing();
    this->gridResolution = data->GetDimensions();
    this->arrayLength = gridResolution[0] * gridResolution[1] * gridResolution[2];
    this->offsetY = gridResolution[0];
    this->offsetZ = gridResolution[0]*gridResolution[1];

    int coutStep = int(double(arrayLength) / 100.0);
    if(coutStep == 0) coutStep = 1;
    
    /* std::vector<std::vector<Vector96d>> accumulatedVecField(numFields, std::vector<Vector96d>(arrayLength, Vector96d::Zero()));
    std::vector<Vector96d> meanVectors(arrayLength, Vector96d::Zero());
    std::vector<Matrix96c> eigenvectorField(arrayLength, Matrix96c::Zero());
    std::vector<std::vector<Vector96d>> sampleField(arrayLength, std::vector<Vector96d>(numSamples, Vector96d::Zero()));
    std::vector<Matrix96d> choleskyField(arrayLength, Matrix96d::Zero()); */

    //clock for random seed and calculation time
    beginning = nanoClock::now();

    //std::vector<Matrix96d> covarianceField(arrayLength, Matrix96d::Zero());
    //std::vector<std::vector<Vector96d>> sampleField2(arrayLength, std::vector<Vector96d>(numSamples, Vector96d::Zero()));
    //std::vector<std::vector<Vector96d>> normalVecField(arrayLength, std::vector<Vector96d>(numSamples, Vector96d::Zero()));
    //std::vector<Eigen::EigenSolver<Matrix96d>> eigensolverField(arrayLength);
    //std::vector<std::vector<std::vector<double>>> testVectors(numFields, std::vector<std::vector<double>>(arrayLength, std::vector<double>(3, 0.0)));    

    /* for(int field = 0; field < numFields; field++){
        vtkSmartPointer<vtkDoubleArray> temptest = vtkSmartPointer<vtkDoubleArray>::New();
        temptest->DeepCopy((vtkImageData::SafeDownCast(input->GetBlock(field))->GetPointData()->GetScalars()));
        for(int i = 0; i < arrayLength; i++){
            double transferTest[3];
            temptest->GetTuple(i, transferTest);
            testVectors[field][i][0] = transferTest[0];
            testVectors[field][i][1] = transferTest[1];
            testVectors[field][i][2] = transferTest[2];
        }
    } */

    cellValues = vtkDoubleArray::New();
    cellValues->SetNumberOfComponents(1);
    cellValues->SetNumberOfTuples(arrayLength);
    cellValues->SetName(this->cellValuesName);
    int calcMean = 0;
    #pragma omp parallel for
    for(int pointIndex = 0; pointIndex < arrayLength; pointIndex++){
        std::vector<Vector96d> neighborhood = std::vector<Vector96d>(numFields, Vector96d::Zero());
        Vector96d meanVector = Vector96d::Zero();
        Matrix96d decomposition = Matrix96d::Zero();

        ++calcMean;
        if(calcMean % coutStep == 0){
            std::chrono::duration<double> calcTime = beginning - nanoClock::now();
            int prog = int((double(calcMean) / double(arrayLength))*100);
            double ETR = (double(calcTime.count())/prog) * (100 - prog);
            cout << '\r' << std::flush;
            cout << "Prog:" << prog << "%, ETR:" << ETR << "s";
        }

        if(isCloseToEdge(pointIndex)){
            cellValues->SetTuple1(pointIndex, 0.0);
            continue;
        }

        std::set<int> indices; //set only contains a single instance of any entitiy

        int nodes[8] = {pointIndex, pointIndex+1, pointIndex+offsetY, pointIndex+offsetY+1, pointIndex+offsetZ,
                            pointIndex+offsetZ+1, pointIndex+offsetY+offsetZ, pointIndex+offsetY+offsetZ+1};
        
        for (int i = 0; i < 8; i++){
            int point = nodes[i];
            int adjacentPoints[7] = {point, point-1, point+1, point-offsetY, point+offsetY, point-offsetZ, point+offsetZ};
            for(int j = 0; j < 7; j++){
                indices.insert(adjacentPoints[j]);
            }
        }
        
        for(int fieldIndex = 0; fieldIndex < numFields; fieldIndex++){
            int c = 0;
            for(auto it = indices.begin(); it != indices.end(); it++, c++){
                double transfer[3];
                int ind = *it;
                //Calculating mean vec
                (vtkImageData::SafeDownCast(input->GetBlock(fieldIndex))->GetPointData()->GetScalars())->GetTuple(ind, transfer);

                meanVector[c*3] += transfer[0];
                meanVector[(c*3)+1] += transfer[1];
                meanVector[(c*3)+2] += transfer[2];

                neighborhood[fieldIndex][c*3]= transfer[0];
                neighborhood[fieldIndex][(c*3)+1]= transfer[1];
                neighborhood[fieldIndex][(c*3)+2]= transfer[2];
            }
        }
        meanVector = meanVector / double(numFields);
        Matrix96d covarMat = Matrix96d::Zero();
        
        for(int fieldIndex = 0; fieldIndex < numFields; fieldIndex++){
            Vector96d transferVec = neighborhood[fieldIndex] - meanVector;
            covarMat += transferVec * transferVec.transpose();
        }

        covarMat = covarMat / double(numFields);

        if(useCholesky){
            Eigen::LDLT<Matrix96d> cholesky(covarMat);

            Matrix96d D = Matrix96d::Identity();
            D.diagonal() = cholesky.vectorD();

            decomposition = cholesky.matrixL() * D;
        } else {
            Eigen::EigenSolver<Matrix96d> eigenMat(covarMat, true);

            Matrix96c evDiag = Matrix96c::Identity();
            evDiag.diagonal() = eigenMat.eigenvalues();

            decomposition = (eigenMat.eigenvectors() * evDiag).real(); //very small (e.g. 10^-16) complex parts are possible, taking real part just to be sure
        }

        int prob = 0;
        for (int sampleIteration = 0; sampleIteration < numSamples; sampleIteration++){
            Vector96d normalVec = generateNormalDistributedVec();
            Vector96d sample = decomposition * normalVec + meanVector;
            if(computeParallelVectors(sample)) prob++;
        }
        double frequency = double(prob) / double(numSamples);
        cellValues->SetTuple1(pointIndex, frequency); 
    }

    vtkSmartPointer<vtkImageData> celldata = vtkSmartPointer<vtkImageData>::New();
    celldata->CopyStructure(data);
    celldata->GetPointData()->AddArray(cellValues);

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkDataObject *output = vtkDataObject::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    output->ShallowCopy(celldata);
    auto t_end = nanoClock::now();
    std::chrono::duration<double> durationTime = t_end - beginning;

    std::cout << "Uncertain Vortex Core Line calculation finished in " << durationTime.count() << " s." << std::endl;

    return 1;
}
    

    /* //Summing up the vectors
    cout << "Calculating mean vectors..." << endl;
    for(int fieldIndex = 0; fieldIndex < numFields; fieldIndex++){

        vtkSmartPointer<vtkDoubleArray> temp = vtkSmartPointer<vtkDoubleArray>::New();
        temp->ShallowCopy((vtkImageData::SafeDownCast(input->GetBlock(fieldIndex))->GetPointData()->GetScalars()));
        
        for(int currentVec = (offsetZ+ offsetY + 1); currentVec < (arrayLength - (((offsetZ * 2) + (offsetY*2)) + 2)); currentVec++){
            
            if((currentVec+2) % gridResolution[0] == 0 or (currentVec+1) % gridResolution[0] == 0 or currentVec % gridResolution[0] == 0 or 
                currentVec % offsetZ < gridResolution[0] or currentVec % offsetZ >= (offsetZ - (gridResolution[0]*2))){
                continue; // Cell at the edge of the field, do nothing
            }
            std::set<int> cellIndices; //set only contains a single instance of any entitiy
            int cellNodeIndices[8] = {currentVec, currentVec+1, currentVec+offsetY, currentVec+offsetY+1, currentVec+offsetZ, currentVec+offsetZ+1,
                                        currentVec+offsetZ+offsetY, currentVec+offsetZ+offsetY+1};
            for (int currentPoint = 0; currentPoint < 8; currentPoint++){
                int point = cellNodeIndices[currentPoint];
                int adjacentPoints[7] = {point, (point - 1), (point + 1), (point - offsetY), (point + offsetY), (point - offsetZ), (point + offsetZ)};
                //get indices of cell nodes and adjacent points
                for(int i = 0; i < 7; i++){
                    cellIndices.insert(adjacentPoints[i]);
                }
            }
            int c = 0;
            for(auto i = cellIndices.begin(); i != cellIndices.end(); i++, c++){
                double transferVector[3];
                int ind = *i;
                //Calculating mean vec
                temp->GetTuple(ind, transferVector);
                meanVectors[currentVec][c*3] += (transferVector[0] / numFields);
                meanVectors[currentVec][(c*3)+1] += (transferVector[1] / numFields);
                meanVectors[currentVec][(c*3)+2] += (transferVector[2] / numFields);

                accumulatedVecField[fieldIndex][currentVec](c*3) = transferVector[0];
                accumulatedVecField[fieldIndex][currentVec]((c*3)+1) = transferVector[1];
                accumulatedVecField[fieldIndex][currentVec]((c*3)+2) = transferVector[2];
            }    
        }
    }
    cout << "Done." << endl;

    cout << "Calculating Eigenvector Matrices..." << endl;
    //int skipped = 0;
    int calculated = 0;
    nanoClock::time_point eigenTimeStart = nanoClock::now();
    #pragma omp parallel for
    for(int cellIndex = 0; cellIndex < arrayLength; cellIndex++){
        calculated++;
        if(calculated % 1000 == 0){
            nanoClock::time_point eigenTimeEnd = nanoClock::now();
            std::chrono::duration<double> eigenTime = eigenTimeEnd - eigenTimeStart;
            cout << "Calculated " << calculated << " Eigenvectors/Cholesky decompositions. The last 1k calculations took " << eigenTime.count() << "s." << endl;
            eigenTimeStart = nanoClock::now();
        }
        Matrix96d covarMat = Matrix96d::Zero();
        
        for(int fieldIndex = 0; fieldIndex < numFields; fieldIndex++){
            Vector96d transferVec = accumulatedVecField[fieldIndex][cellIndex] - meanVectors[cellIndex];
            covarMat += (transferVec * transferVec.transpose()) / numFields;
        }

        //covarianceField[cellIndex] = covarMat;

        if(useCholesky){
            Eigen::LLT<Matrix96d> cholesky(covarMat);
            choleskyField[cellIndex] = cholesky.matrixU();
        
        } else {
            Eigen::EigenSolver<Matrix96d> eigenMat(covarMat, true);
            //eigensolverField[cellIndex] = eigenMat;
            Matrix96c scaledEigenvecs = Matrix96c::Zero();
            for(int eigenVec = 0; eigenVec < 96; eigenVec++){
                //eigenvectors scaled by their eigenvalues
                scaledEigenvecs.col(eigenVec) = eigenMat.eigenvectors().col(eigenVec) * eigenMat.eigenvalues().row(eigenVec);
            }
            eigenvectorField[cellIndex] = scaledEigenvecs;
        }
    }
    cout << "Done." << endl;
    //cout << "Skipped " << skipped << " calculations, " << (double(skipped)/arrayLength)*100 << "%" << endl;

    if(useRandomSeed){
        nanoClock::duration d = nanoClock::now() - beginning;
        unsigned seed = d.count();
        gen.seed(seed);
    } else {
        this->gen.seed(12);
    }
    //Monte Carlo Sampling
    cout << "Generating Monte Carlo Samples..." << endl;
    #pragma omp parallel for
    for(int cellIndex = 0; cellIndex < arrayLength; cellIndex++){
        
        for (int sampleIteration = 0; sampleIteration < numSamples; sampleIteration++){
            
            Vector96d normalVec = generateNormalDistributedVec();
            //normalVecField[cellIndex]d[sampleIteration] = normalVec;
            if(useCholesky){
                Vector96d sample = (choleskyField[cellIndex].transpose() * normalVec) + meanVectors[cellIndex];
                sampleField[cellIndex][sampleIteration] = sample;
            } else {
                Vector96c sample = (eigenvectorField[cellIndex] * normalVec.cast<std::complex<double>>()) + meanVectors[cellIndex].cast<std::complex<double>>();
                sampleField[cellIndex][sampleIteration] = sample.real();
            }
        }
    }
    cout << "Done." << endl;



    cellValues = vtkDoubleArray::New();
    cellValues->SetNumberOfComponents(1);
    cellValues->SetNumberOfTuples(arrayLength);
    cellValues->SetName(this->cellValuesName);
    int parallel;
    double frequency;

    //Calculating parallel vectors in cells
    cout << "Computing Parallel Vectors..." << endl;
    for(int cellIndex = 0; cellIndex < arrayLength; cellIndex++){
        parallel = 0;
        frequency = 0.0;
        for(int sample = 0; sample < numSamples; sample++){
            bool hasParallel = false;
            hasParallel = computeParallelVectors(sampleField[cellIndex][sample]);
            if(hasParallel){
                parallel++;
                //cout << "Found parallel vectors" << endl;
            }
        }
        frequency = double(parallel) / double(numSamples);
        if (frequency > 0){
            cout << cellIndex << ": " << frequency << " - " << parallel << " - " << numSamples << endl;
        }
        cellValues->SetTuple1(cellIndex, frequency);
    }
    cout << "Done." << endl;

    vtkSmartPointer<vtkImageData> celldata = vtkSmartPointer<vtkImageData>::New();
    celldata->CopyStructure(data);
    celldata->GetPointData()->AddArray(cellValues);

    vtkInformation *outInfo = outputVector->GetInformationObject(0);
    vtkDataObject *output = vtkDataObject::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

    output->ShallowCopy(celldata);
    auto t_end = nanoClock::now();
    std::chrono::duration<double> durationTime = t_end - beginning;

    std::cout << "Uncertain Vortex Core Line calculation finished in " << durationTime.count() << " s." << std::endl;

    return 1; */
//}

bool UncertainVortexCores::isCloseToEdge(int index){
    bool isClose = false;

    //check if index is close to an edge in x direction
    if(((index+2) % this->offsetY == 0) or (((index+1) % this->offsetY) == 0) or ((index % this->offsetY) == 0)){
        isClose = true;
    }
    //check if index is close to an edge in y direction
    if(((index % this->offsetZ) < this->offsetY) or ((index % this->offsetZ) >= (this->offsetZ - (this->offsetY*2)))){
        isClose = true;
    }
    //check if index is close to an edge in z direction
    if((index < offsetZ) or (index >= (arrayLength - (offsetZ*2)))){
        isClose = true;
    }

    return isClose;
}

Vector96d UncertainVortexCores::generateNormalDistributedVec(){

    std::uniform_real_distribution<double> randomVecEntry(0.00001, 1.0); 
    
    Vector96d normalVec = Vector96d::Zero();

    for(int pair = 0; pair < 48; pair++){
            
        double u1 = randomVecEntry(this->gen);
        double u2 = randomVecEntry(this->gen);
        //Box Muller Transformation
        double z1 = sqrt(-2.0*log(u1))*cos((2*M_PI)*u2);
        double z2 = sqrt(-2.0*log(u1))*sin((2*M_PI)*u2);
        normalVec(pair*2) = z1;
        normalVec((pair*2)+1) = z2;
    }

    return normalVec;
}

std::vector<std::tuple<Vector3d, Matrix3d>> UncertainVortexCores::calculateJacobiAndCellNodes(Vector96d sampleVector){
    //Takes a 96 dimensional sample vector, splits them into the cell nodes and calculates the respective jacobis
    std::vector<Vector3d> points(32, Vector3d::Zero());
    std::vector<std::tuple<Vector3d, Matrix3d>> result(8);
    
    for (int point = 0; point < 32; point++){
        points[point][0] = sampleVector[(point*3)];
        points[point][1] = sampleVector[(point*3)+1];
        points[point][2] = sampleVector[(point*3)+2];
    }
    //vectors in sampleVector are sorted ascending by their indices in the original vectorfields, nodeIndices has the indices of the cell nodes at the first place of every subarray,
    //followed by x+1, x-1, y+1, y-1, z+1, z-1 for easy jacobi calculation, the order of cell nodes followes the spatial arrangement from bottom left to top right
    int nodeIndices[8][7] = {{7,8,6,11,4,19,0},{8,9,7,12,5,20,1},{11,12,10,14,7,23,2},{12,13,11,15,8,24,3},
                            {19,20,18,23,16,28,7},{20,21,19,24,17,29,8},{23,24,22,26,19,30,11},{24,25,23,27,20,31,12}};

    for(int node = 0; node < 8; node++){
        Matrix3d jacobi;
        jacobi.col(0) = (points[nodeIndices[node][1]] - points[nodeIndices[node][2]]) / (2 * spacing[0]);
        jacobi.col(1) = (points[nodeIndices[node][3]] - points[nodeIndices[node][4]]) / (2 * spacing[1]);
        jacobi.col(2) = (points[nodeIndices[node][5]] - points[nodeIndices[node][6]]) / (2 * spacing[2]);
        result[node] = std::make_tuple(points[nodeIndices[node][0]], jacobi);
    }

    return result;
}

std::vector<Vector3d> UncertainVortexCores::calculateJacobiTimesVelocity(std::vector<std::tuple<Vector3d, Matrix3d>> cell) {
    //calculate acceleration of cell nodes
    std::vector<Vector3d> acceleration(8, Vector3d::Zero());
    
    for(int i = 0; i < 8; i++){
        Vector3d velocity = std::get<0>(cell[i]);
        Matrix3d jacobi = std::get<1>(cell[i]);
        Vector3d accel = jacobi * velocity;
        acceleration[i] = accel;
    }
    
    return acceleration;
}

bool UncertainVortexCores::computeParVectorsInCell(std::vector<std::tuple<Vector3d, Matrix3d>> cell){
    
    double s[3];
    double t[3];
    int numParal = 0;
    int faceWithPar = 0;
    std::vector<Vector3d> acceleration = calculateJacobiTimesVelocity(cell);
    
    //indices of every triangular face of a hexahedron with 0-7 being the indices of the cell nodes, 6 cellfaces with 2 triangles each, 3 nodes per triangle
    int cellfaceIndices[6][2][3] = {{{0,1,3},{0,2,3}},
                                    {{0,1,5},{0,4,5}},
                                    {{0,2,6},{0,4,6}},
                                    {{1,3,5},{3,5,7}},
                                    {{2,3,6},{3,6,7}},
                                    {{4,5,6},{5,6,7}}};
    
    std::vector<Vector3d> cellfaceVel(3, Vector3d::Zero());
    std::vector<Vector3d> cellfaceAcc(3, Vector3d::Zero());
    
    for(int cellface = 0; cellface < 6; cellface++){
        
        bool isComplex = false;
        for(int tria = 0; tria < 2; tria++){
            
            cellfaceVel[0] = std::get<0>(cell[cellfaceIndices[cellface][tria][0]]);
            cellfaceVel[1] = std::get<0>(cell[cellfaceIndices[cellface][tria][1]]);
            cellfaceVel[2] = std::get<0>(cell[cellfaceIndices[cellface][tria][2]]);
            cellfaceAcc[0] = acceleration[cellfaceIndices[cellface][tria][0]];
            cellfaceAcc[1] = acceleration[cellfaceIndices[cellface][tria][1]];
            cellfaceAcc[2] = acceleration[cellfaceIndices[cellface][tria][2]];

            numParal = computeParallelOnCellface(cellfaceVel, cellfaceAcc, s, t);
            
            for(int par = 0; par < numParal; par++){
                Matrix3d jac1 = std::get<1>(cell[cellfaceIndices[cellface][tria][0]]);
                Matrix3d jac2 = std::get<1>(cell[cellfaceIndices[cellface][tria][1]]);
                Matrix3d jac3 = std::get<1>(cell[cellfaceIndices[cellface][tria][2]]);
                isComplex = interpolatedJacobiIsComplex(jac1, jac2, jac3, s[par], t[par]);
            }
        }

        if(isComplex){
            faceWithPar++;
        }
    }

    if (faceWithPar > 1){
        return true;
    } else {
        return false;
    }   
}

int UncertainVortexCores::computeParallelOnCellface(std::vector<Vector3d> cellfaceVel, std::vector<Vector3d> cellfaceAcc, double *s, double *t){

    //now using linalg because of original code and efficiency
    vec3 v0 = {cellfaceVel[0][0], cellfaceVel[0][1], cellfaceVel[0][2]};
    vec3 v1 = {cellfaceVel[1][0], cellfaceVel[1][1], cellfaceVel[1][2]};
    vec3 v2 = {cellfaceVel[2][0], cellfaceVel[2][1], cellfaceVel[2][2]};

    vec3 w0 = {cellfaceAcc[0][0], cellfaceAcc[0][1], cellfaceAcc[0][2]};
    vec3 w1 = {cellfaceAcc[1][0], cellfaceAcc[1][1], cellfaceAcc[1][2]};
    vec3 w2 = {cellfaceAcc[2][0], cellfaceAcc[2][1], cellfaceAcc[2][2]};
    
    vec3 v01, v02;
    vec3 w01, w02;
    mat3 V, Vinv;
    mat3 W, Winv;
    mat3 M; //Matrix for eigenvector problem, either W^-1 * V or V^-1 * W

    double eigenvalues[3];
    vec3 eigenvectors[3];

    double detV, detW;
    double absdetV, absdetW, absdetmax;
    double nx, ny, nz;
    double ss, tt;
    int numParal;
    int numEigen;
    int ok ,take;

    // The vectors v0->v1 and v0->v2 span the triangle.
	// The vectors v0,v01,v02 are the columns of the V matrix.

    vec3sub(v1, v0, v01);
    vec3sub(v2, v0, v02);
    vec3sub(w1, w0, w01);
    vec3sub(w2, w0, w02);

    mat3setcols(V, v0, v01, v02);
    mat3setcols(W, w0, w01, w02);

    detW = mat3det(W);
	detV = mat3det(V);

	absdetW = fabs(detW);
	absdetV = fabs(detV);

    //Taking matrix with larger determinant
    take = 0;
    absdetmax = 0.0;

    if (absdetW > absdetmax){
			take = 1;
			absdetmax = absdetW;
	}
    if (absdetV > absdetmax) take = 2;

    switch (take) {
        case 0:
            //Matrices not invertible
            return 0;

        case 1:
            mat3invdet(W, detW, Winv);
            mat3mul(Winv, V, M);
            break;
        case 2:
            mat3invdet(V, detV, Vinv);
            mat3mul(Vinv, W, M);
            break;
    }

    numParal = 0;
    numEigen = mat3eigenvalues(M, eigenvalues);
    
    for (int i = 0; i < numEigen; i++){

        ok = mat3realEigenvector(M, eigenvalues[i], eigenvectors[i]);
        //invert eigenvalues if V got inverted
        if(take == 2) {
            if (eigenvalues[i] == 0.0){
                ok = 0;
            } else eigenvalues[i] = 1.0 / eigenvalues[i];
        }

        if(ok){
            //scale the normed eigenvector (nx,ny,nz) to length (1,s,t)
            nx = eigenvectors[i][0];
            ny = eigenvectors[i][1];
            nz = eigenvectors[i][2];

            if (nx != 0.0){
                //local coords in triangle

                ss = ny / nx;
                tt = nz / nx;

                //check if point is inside the triangle
                if ((ss >= 0) and (tt >= 0) and (ss + tt <= 1)){
                    s[numParal] = ss;
                    t[numParal] = tt;
                    numParal++;
                }

            }
        }
    }
    
    return numParal;
}

bool UncertainVortexCores::interpolatedJacobiIsComplex(Matrix3d jac1, Matrix3d jac2, Matrix3d jac3, double s, double t){


    vec3 A;
    vec3 B;
    vec3 C;
    vec3 result[3];
    mat3 interpolated;
    double eigenvalues[3];
    for(int i = 0; i < 3; i++){
        A[0] = jac1.col(i)[0]; A[1] = jac1.col(i)[1]; A[2] = jac1.col(i)[2];
        B[0] = jac2.col(i)[0]; B[1] = jac2.col(i)[1]; B[2] = jac2.col(i)[2];
        C[0] = jac3.col(i)[0]; C[1] = jac3.col(i)[1]; C[2] = jac3.col(i)[2];
        vec3lerp3(A, B, C, s, t, result[i]);
    }

    mat3setcols(interpolated, result[0], result[1], result[2]);
    int numRealEigenvalues = mat3eigenvalues(interpolated, eigenvalues);
    bool isComplex = (numRealEigenvalues == 3 ? false : true);
    return isComplex;
}

bool UncertainVortexCores::computeParallelVectors(Vector96d sampleVector){
    
    bool hasParallel = false;
    std::vector<std::tuple<Vector3d, Matrix3d>> nodesAndJacobi = calculateJacobiAndCellNodes(sampleVector);
    hasParallel = computeParVectorsInCell(nodesAndJacobi);
    
    return hasParallel;
}

