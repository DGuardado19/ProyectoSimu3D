node selectNode(int i, element e,mesh &m){
	node n;
	switch(i){
		case 1: n = m.getNode(e.getNode1()-1); break;
		case 2: n = m.getNode(e.getNode2()-1); break;
		case 3: n = m.getNode(e.getNode3()-1); break;
        case 4: n = m.getNode(e.getNode4()-1); break;
	}
	return n;
}

float selectCoord(int c, node n){
	float v;
	switch(c){
		case EQUIS: v = n.getX(); break;
		case YE:    v = n.getY(); break;
        case ZETA:  v = n.getZ(); break;
	}
	return v;
}

float calcularTenedor(element e, int coord, int i, int j,mesh &m){
	node n1=selectNode(i,e,m), n2=selectNode(j,e,m);
	return selectCoord(coord,n1) - selectCoord(coord,n2);
}

float calculateMagnitude(float v1, float v2, float v3){
    return sqrt(pow(v1,2)+pow(v2,2)+pow(v3, 2));
}

float OperarRestaTenedor(element e, int coord1, int coord2,  float value_a, float value_b, mesh &m){
    float a, b, c, d;

    a = calcularTenedor(e, coord1, value_a, 1 ,m);
    b = calcularTenedor(e, coord2, value_b, 1 ,m);
    c = calcularTenedor(e, coord1, value_b, 1 ,m);
    d = calcularTenedor(e, coord2, value_a, 1 ,m);

    return (a*b)-(c*d);
}

//Calculando partes de matriz A
float calcularA1(element e,mesh &m){
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

	a=4*(x3*x3)+4*(x2+x4+6*x1)*x3-15*y1-3*y3-3*y2+4*x2*x2+4*x2*(x4+6*x1)-3*y4+4*x4*x2+24*x1*x1+3*(20*x1*x1-z2-z3-z4-5*z1);
	resultado = a/360;
	return resultado;
}
float calcularA2(element e,mesh &m){
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

	a=4*x3*x3+4*(2*x2+x4+6*x1)*x3-15*y1-3*y3-6*y2+12*x2*x2+8*x2*(x4+6*x1)-3*y4+4*x4*x4+24*x4*x1+3*(20*x1*x1-2*z2-z3-z4-5*z1);

	resultado = a/360;
	return resultado;
}
float calcularA3(element e,mesh &m){
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

	a=12*x3*x3+8*(x2+x4+6*x1)*x3-15*y1-6*y3-3*y2+4*x2*x2+4*x2*(x4+6*x1)-3*y4+4*x4*x4+24*x4*x1+3*(20*x1*x1-z2-2*z3-z4-5*z1);

	resultado = a/360;
	return resultado;
}
float calcularA4(element e,mesh &m){
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

	a=4*x3*x3+4*(x2+2*(x4+3*x1))*x3-15*y1-3*y3-3*y2+4*x2*x2+8*x2*(x4+3*x1)-3*(2*y4-4*x4*x4-16*x4*x1-20*x1*x1+z2+z3+2*z4+5*z1);

	resultado = a/360;
	return resultado;
}

//Calculando partes de Matriz G
float calcularG1(element e,mesh &m){
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

	a=x3-5*y1-y3-y2+x2-y4+x4+5*x1+z2+z3+z4+5*z1;
	resultado = a/120;
	return resultado;
}
float calcularG2(element e,mesh &m){
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

	a=x3-5*y1-y3-2*y2+2*x2-y4+x4+5*x1+2*z2+z3+z4+5*z1;

	resultado = a/120;
	return resultado;
}
float calcularG3(element e,mesh &m){
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

	a=2*x3-5*y1-2*y3-y2+x2-y4+x4+5*x1+z2+2*z3+z4+5*z1;

	resultado = a/120;
	return resultado;
}
float calcularG4(element e,mesh &m){
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

	a=x3-5*y1-y3-y2+x2-2*y4+2*x4+5*x1+z2+z3+2*z4+5*z1;

	resultado = a/120;
	return resultado;
}

float calculateLocalVolume(int i,mesh m){

    double Ve, u, v, w, U, V, W, a, b, c, d, X, x, Y, y, Z, z;
    
    element e = m.getElement(i);
    node n1 = m.getNode(e.getNode1()-1);
    node n2 = m.getNode(e.getNode2()-1);
    node n3 = m.getNode(e.getNode3()-1);
    node n4 = m.getNode(e.getNode4()-1);

    U = calculateMagnitude(n2.getX()-n1.getX(), n2.getY()-n1.getY(), n2.getZ()-n1.getZ() );
    V = calculateMagnitude(n3.getX()-n2.getX(), n3.getY()-n2.getY(), n3.getZ()-n2.getZ() );
    W = calculateMagnitude(n3.getX()-n1.getX(), n3.getY()-n1.getY(), n3.getZ()-n1.getZ() );

    u = calculateMagnitude(n4.getX()-n3.getX(), n4.getY()-n3.getY(), n4.getZ()-n3.getZ() );
    v = calculateMagnitude(n4.getX()-n1.getX(), n4.getY()-n1.getY(), n4.getZ()-n1.getZ() );
    w = calculateMagnitude(n4.getX()-n2.getX(), n4.getY()-n2.getY(), n4.getZ()-n2.getZ() );

    X = (w-U+v)*(U+v+w);
    x = (U-v+w)*(v-w+U);

    Y = (u-V+w)*(V+w+u);
    y = (V-w+u)*(w-u+V);

    Z = (v-W+u)*(W+u+v);
    z = (W-u+v)*(u-v+W);

    a = sqrt(x*Y*Z);
    b = sqrt(y*Z*X);
    c = sqrt(z*X*Y);
    d = sqrt(x*y*z);

    Ve = sqrt( (-a+b+c+d)*(a-b+c+d)*(a+b-c+d)*(a+b+c-d) ) / (192*u*v*w);
    
    return Ve;

}

void ubicarSubMatriz(Matrix &K,int fi,int ff,int ci,int cf,Matrix M){
    int n = 0, m= 0;
    for(int i=fi;i<=ff;i++){
        for(int j=ci;j<=cf;j++){
            K.at(i).at(j) = M.at(n).at(m);
            m++;
        }
        n++; m = 0;
    }
}
