float calculateLocalD(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, YE, 2, 1, m));
    row1.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row2.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, EQUIS, 4, 1, m));
    row3.push_back(calcularTenedor(e, YE, 4, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

void calculateAlpha(int i,Matrix &A,mesh m){
    zeroes(A,3);
    element e = m.getElement(i);
    
    A.at(0).at(0) = OperarRestaTenedor(e, YE, ZETA, 3, 4, m);
    A.at(0).at(1) = OperarRestaTenedor(e, YE, ZETA, 4, 2, m);
    A.at(0).at(2) = OperarRestaTenedor(e, YE, ZETA, 2, 3, m);

    A.at(1).at(0) = OperarRestaTenedor(e, EQUIS, ZETA, 4, 3, m);
    A.at(1).at(1) = OperarRestaTenedor(e, EQUIS, ZETA, 2, 4, m);
    A.at(1).at(2) = OperarRestaTenedor(e, EQUIS, ZETA, 3, 2, m);

    A.at(2).at(0) = OperarRestaTenedor(e, EQUIS, YE, 3, 4, m);
    A.at(2).at(1) = OperarRestaTenedor(e, EQUIS, YE, 4, 2, m);
    A.at(2).at(2) = OperarRestaTenedor(e, EQUIS, YE, 2, 3, m);

}

void calculateBeta(Matrix &B){
    zeroes(B,3,12);

    B.at(0).at(0) = -1; 
    B.at(0).at(1) =  1; 
    B.at(0).at(2) =  0; 
    B.at(0).at(3) =  0; 
    B.at(0).at(4) = -1; 
    B.at(0).at(5) =  1; 
    B.at(0).at(6) =  0; 
    B.at(0).at(7) =  0; 
    B.at(0).at(8) = -1; 
    B.at(0).at(9) =  1; 
    B.at(0).at(10) =  0; 
    B.at(0).at(11) =  0; 
    
    B.at(1).at(0) = -1; 
    B.at(1).at(1) = 0; 
    B.at(1).at(2) = 1;
    B.at(1).at(3) = 0;
    B.at(1).at(4) = -1; 
    B.at(1).at(5) = 0; 
    B.at(1).at(6) = 1;
    B.at(1).at(7) = 0;
    B.at(1).at(8) = -1; 
    B.at(1).at(9) = 0; 
    B.at(1).at(10) = 1;
    B.at(1).at(11) = 0;

    B.at(2).at(0) = -1;
    B.at(2).at(1) = 0;
    B.at(2).at(2) = 0;
    B.at(2).at(3) = 1;
    B.at(2).at(4) = -1;
    B.at(2).at(5) = 0;
    B.at(2).at(6) = 0;
    B.at(2).at(7) = 1;
    B.at(2).at(8) = -1;
    B.at(2).at(9) = 0;
    B.at(2).at(10) = 0;
    B.at(2).at(11) = 1;
}

void calculateOmega(Matrix &C){
    zeroes(C,3,4);
    C.at(0).at(0) = -1; C.at(0).at(1) = 1; C.at(0).at(2) = 0; C.at(0).at(3) = 0;
    C.at(1).at(0) = -1; C.at(1).at(1) = 0; C.at(1).at(2) = 1; C.at(1).at(3) = 0;
    C.at(2).at(0) = -1; C.at(2).at(1) = 0; C.at(2).at(2) = 0; C.at(2).at(3) = 1;
}


void calculateGamma(Matrix &m){
	zeroes(m,12,3);

	m.at(0).at(0) = 1;   m.at(0).at(1) = 0;   m.at(0).at(2) = 0;
	m.at(1).at(0) = 1;   m.at(1).at(1) = 0;   m.at(1).at(2) = 0; 
    m.at(2).at(0) = 1;   m.at(2).at(1) = 0;   m.at(2).at(2) = 0;
	m.at(3).at(0) = 1;   m.at(3).at(1) = 0;   m.at(3).at(2) = 0; 
    m.at(4).at(0) = 0;   m.at(4).at(1) = 1;   m.at(4).at(2) = 0;
	m.at(5).at(0) = 0;   m.at(5).at(1) = 1;   m.at(5).at(2) = 0; 
    m.at(6).at(0) = 0;   m.at(6).at(1) = 1;   m.at(6).at(2) = 0;
	m.at(7).at(0) = 0;   m.at(7).at(1) = 1;   m.at(7).at(2) = 0; 
    m.at(8).at(0) = 0;   m.at(8).at(1) = 0;   m.at(8).at(2) = 1;
	m.at(9).at(0) = 0;   m.at(9).at(1) = 0;   m.at(9).at(2) = 1; 
    m.at(10).at(0) = 0;  m.at(10).at(1) = 0;  m.at(10).at(2) = 1;
	m.at(11).at(0) = 0;  m.at(11).at(1) = 0;  m.at(11).at(2) = 1; 
	
}
//Calculando componente A
void calculateFinal(int i,Matrix &R,mesh m){
	zeroes(R,12,3);
    element e = m.getElement(i);
    float a,b,c,d;
    a = calcularA1(e,m);
    b = calcularA2(e,m);
    c = calcularA3(e,m);
    d = calcularA4(e,m);
    R.at(0).at(0) = a; R.at(0).at(1) = 0; R.at(0).at(2) = 0;
    R.at(1).at(0) = b; R.at(1).at(1) = 0; R.at(1).at(2) = 0;
    R.at(2).at(0) = c; R.at(2).at(1) = 0; R.at(2).at(2) = 0;
    R.at(3).at(0) = d; R.at(3).at(1) = 0; R.at(3).at(2) = 0;
    R.at(4).at(0) = 0; R.at(4).at(1) = a; R.at(4).at(2) = 0;
    R.at(5).at(0) = 0; R.at(5).at(1) = b; R.at(5).at(2) = 0;
    R.at(6).at(0) = 0; R.at(6).at(1) = c; R.at(6).at(2) = 0;
    R.at(7).at(0) = 0; R.at(7).at(1) = d; R.at(7).at(2) = 0;
    R.at(8).at(0) = 0; R.at(8).at(1) = 0; R.at(8).at(2) = a;
    R.at(9).at(0) = 0; R.at(9).at(1) = 0; R.at(9).at(2) = b;
    R.at(10).at(0) = 0; R.at(10).at(1) = 0; R.at(10).at(2) = c;
    R.at(11).at(0) = 0; R.at(11).at(1) = 0; R.at(11).at(2) = d; 
	
}

///componente k, es solo una integral
float calculateK(int i,mesh m){
    element e = m.getElement(i);
	float a, resultado, x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4;
	node n1=selectNode(1,e,m); 

	x1=selectCoord(EQUIS,n1);
	x2 = calcularTenedor(e,EQUIS,2,1,m);
	x3 = calcularTenedor(e,EQUIS,3,1,m);
	x4 = calcularTenedor(e,EQUIS,4,1,m);
	y1=selectCoord(YE,n1);
	y2 = calcularTenedor(e,YE,2,1,m);
	y3 = calcularTenedor(e,YE,3,1,m);
	y4 = calcularTenedor(e,YE,4,1,m);
	z1=selectCoord(ZETA,n1);
	z2 = calcularTenedor(e,ZETA,2,1,m);
	z3 = calcularTenedor(e,ZETA,3,1,m);
	z4 = calcularTenedor(e,ZETA,4,1,m);
	
	a=3*x3-15*y1*y1-6*y1*(y3+y2+y4)-y3*y3-y3*(y2+y4)-y2*y2-y4*y2+3*x2-y4*y4+3*x4+15*x1-z2*x2-z2*(z3+z4+6*z1)-z3*z3-z3*(z4+6*z1)-z3*z3-6*z4*z1-15*z1*z1;

	resultado = a/360;
	return resultado;
	
}
//componenteG, N*t

void calculateG(int i,Matrix &R,mesh m){
	zeroes(R,12,3);
    element e = m.getElement(i);
    float a,b,c,d;
    a = calcularG1(e,m);
    b = calcularG2(e,m);
    c = calcularG3(e,m);
    d = calcularG4(e,m);
    R.at(0).at(0) = a; R.at(0).at(1) = 0; R.at(0).at(2) = 0;
    R.at(1).at(0) = b; R.at(1).at(1) = 0; R.at(1).at(2) = 0;
    R.at(2).at(0) = c; R.at(2).at(1) = 0; R.at(2).at(2) = 0;
    R.at(3).at(0) = d; R.at(3).at(1) = 0; R.at(3).at(2) = 0;
    R.at(4).at(0) = 0; R.at(4).at(1) = a; R.at(4).at(2) = 0;
    R.at(5).at(0) = 0; R.at(5).at(1) = b; R.at(5).at(2) = 0;
    R.at(6).at(0) = 0; R.at(6).at(1) = c; R.at(6).at(2) = 0;
    R.at(7).at(0) = 0; R.at(7).at(1) = d; R.at(7).at(2) = 0;
    R.at(8).at(0) = 0; R.at(8).at(1) = 0; R.at(8).at(2) = a;
    R.at(9).at(0) = 0; R.at(9).at(1) = 0; R.at(9).at(2) = b;
    R.at(10).at(0) = 0; R.at(10).at(1) = 0; R.at(10).at(2) = c;
    R.at(11).at(0) = 0; R.at(11).at(1) = 0; R.at(11).at(2) = d; 
	
}


float calculateLocalJ(int i,mesh m){
    Matrix matrix;
    Vector row1, row2, row3;

    element e = m.getElement(i);
    row1.push_back(calcularTenedor(e, EQUIS, 2, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 3, 1, m));
    row1.push_back(calcularTenedor(e, EQUIS, 4, 1, m));

    row2.push_back(calcularTenedor(e, YE, 2, 1, m));
    row2.push_back(calcularTenedor(e, YE, 3, 1, m));
    row2.push_back(calcularTenedor(e, YE, 4, 1, m));

    row3.push_back(calcularTenedor(e, ZETA, 2, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 3, 1, m));
    row3.push_back(calcularTenedor(e, ZETA, 4, 1, m));

    matrix.push_back(row1);
    matrix.push_back(row2);
    matrix.push_back(row3);

    return determinant(matrix);
}

Matrix createLocalM(int e,mesh &m){
    Matrix matrixA,matrixK,matrixG,matrixD;
    float J,Determinant;
    
    /* [ A+K  G ]
       [  D   0 ]
    */

    //Matrix A
    /*
    	1/D * J * A * alpha * Beta
    */
    Matrix matrizR, Alpha, Beta;

    Determinant = calculateLocalD(e,m);
    J = calculateLocalJ(e,m);
	
    if(Determinant == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    float real_a = (float) (J)/(Determinant);
    calculateFinal(e,matrizR,m);
    calculateAlpha(e,Alpha,m);
    calculateBeta(Beta);
    // real_a(J/D) * A * alpha * Beta 
    productRealMatrix(real_a, productMatrixMatrix(matrizR,productMatrixMatrix(Alpha,Beta,3,3,12),12,3,12),matrixA);
    
    //Matrix K
    Matrix Alpha_t,Beta_t;
	
    float kmatrix = calculateK(e,m);
	// K' / D*D    
    float real_k = (float) (kmatrix)/(Determinant*Determinant);

    transpose(Alpha,Alpha_t);
    transpose(Beta,Beta_t);
	// (K' / D*D) * Betat*alphaT * alpha *beta  
    productRealMatrix(real_k,productMatrixMatrix(Beta_t,productMatrixMatrix(Alpha_t,productMatrixMatrix(Alpha,Beta,3,3,12),3,3,12),12,3,12),matrixK);

    
    //Matrix G
    Matrix Omega,matrizg;
	
    float real_g = (float) (J/Determinant);
    calculateOmega(Omega);
    calculateG(e,matrizg,m);
    // (J/D) * G * alpha * omega
    productRealMatrix(real_g,productMatrixMatrix(matrizg,productMatrixMatrix(Alpha,Omega,3,3,4),12,3,4),matrixG);
    ///12x4
    
    //Matrix D No cambia al del codigo 
    Matrix g_matrix,g_matrix_t,Omega_t;
    // J/24 
    float real_d = (float)(J/(24*Determinant));
	calculateGamma(g_matrix);
    transpose(Omega, Omega_t);
    transpose(g_matrix,g_matrix_t);
    // (J/24*D) * omegaT * alphaT * gammaT
    productRealMatrix(real_d,productMatrixMatrix(Omega_t,productMatrixMatrix(Alpha_t,g_matrix_t,3,3,12),4,3,12),matrixD);

    //Matrix local Global
    Matrix M;
    zeroes(M,16);
    //primeras filas y luego columnas
    //16x15 A+K

    ubicarSubMatriz(M,0,11,0,11, sumMatrix(matrixA,matrixK,12,12));
    ubicarSubMatriz(M,0,11,12,15,matrixG);
    ubicarSubMatriz(M,12,15,0,11,matrixD);

    return M;
}

void calculateF(Vector &f, mesh &m){
    zeroes(f,3);

    f.at(0) += m.getParameter(EXTERNAL_FORCE_X);
    f.at(1) += m.getParameter(EXTERNAL_FORCE_Y);
    f.at(2) += m.getParameter(EXTERNAL_FORCE_Z);

}
//f vectors 
Vector createLocalb(int e,mesh &m){
    float J;
    element ele = m.getElement(e);
    
    Vector b,b_aux,f;
    Matrix g_matrix;
	
    calculateF(f, m);
    calculateGamma(g_matrix);

    J = calculateLocalJ(e,m);

    if(J == 0){
        cout << "\n!---CATASTROPHIC FAILURE---!\n";
        exit(EXIT_FAILURE);
    }
    
    
    zeroes(b_aux,16);
    productMatrixVector(g_matrix,f,b_aux);
    productRealVector(J/24,b_aux,b);
    return b;
}
