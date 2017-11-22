/*
 *  BesselJ.cpp
 *  Genesis
 *
 *  Created by Sven Reiche on 12/5/11.
 *  Copyright 2011 Paul Scherrer Institut. All rights reserved.
 *
 */

#include "BesselJ.h"

BesselJ::BesselJ(){}
BesselJ::~BesselJ(){}


double BesselJ::BesselJ0(double x)
{
    const double     p1=1.e0;
    const double     p2=-.1098628627e-2;
    const double     p3=.2734510407e-4;
    const double     p4= -.2073370639e-5;
    const double     p5=.2093887211e-6;
    const double     q1=-.1562499995e-1;
    const double     q2=.1430488765e-3;
    const double     q3=-.6911147651e-5;
    const double     q4=.7621095161e-6;
    const double     q5=-.934945152e-7;
    const double     r1=57568490574.e0;
    const double     r2=-13362590354.e0;
    const double     r3=651619640.7e0;
    const double     r4=-11214424.18e0;
    const double     r5=77392.33017e0;
    const double     r6=-184.9052456e0;
    const double     s1=57568490411.e0;
    const double     s2=1029532985.e0;
    const double     s3=9494680.718e0;
    const double     s4=59272.64853e0;
    const double     s5=267.8532712e0;
    const double     s6=1.e0;   
    
	if (fabs(x) < 8){
		double y=x*x;
		return (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))));        
	}
	else
	{
		double ax=fabs(x);
		double z=8/ax;
		double y=z*z;
		double xx=ax-0.785398164;
		return sqrt(0.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-
									 z*sin(xx)* (q1+y*(q2+y*(q3+y*(q4+y*q5)))));       
	}
}


double BesselJ::BesselJ1(double x)
{
    const double    r1=72362614232.e0;
    const double    r2=-7895059235.e0;
    const double    r3=242396853.1e0;
    const double    r4=-2972611.439e0;
    const double    r5=15704.48260e0;
    const double    r6=-30.16036606e0;
    const double    s1=144725228442.e0;
    const double    s2=2300535178.e0;
    const double    s3=18583304.74e0;
    const double    s4=99447.43394e0;
    const double    s5=376.9991397e0;
    const double    s6=1.e0;
    const double    p1=1.e0;
    const double    p2=.183105e-2;
    const double    p3=-.3516396496e-4;
    const double    p4=.2457520174e-5;
    const double    p5=-.240337019e-6;
    const double    q1=.04687499995e0;
    const double    q2=-.2002690873e-3;
    const double    q3=.8449199096e-5;
    const double    q4=-.88228987e-6;
    const double    q5=.105787412e-6;
	
	if (fabs(x) < 8){
		double y=x*x;
		return x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))));        
	}
	else
	{
		double ax=fabs(x);
		double z=8/ax;
		double y=z*z;
		double xx=ax-2.356194491;
		return sqrt(0.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-
									 z*sin(xx)* (q1+y*(q2+y*(q3+y*(q4+y*q5)))));       
	}
}





double BesselJ::value(int n, double x)
{

	
	double bessj=0;
	if (n==0){
  	  return this->BesselJ0(x);
	}
	if (n==1) {
	  return this->BesselJ1(x);
	}
	
	if (x==0) {
		return 0;
	}
	
    double ax = fabs(x);
	
	if (ax > static_cast<double> (n)){
		double tox=2/ax;
		double bjm=this->BesselJ0(ax);
		double bj =this->BesselJ1(ax);
		for (int j=1; j<n; j++) {
			double bjp=j*tox*bj-bjm;
			bjm=bj;
			bj=bjp;
		}
	    bessj=bj;
	} else {
		double tox=2./ax;
		int m=2*(n+floor(sqrt(static_cast<double>(40*n)))/2);
	    bessj=0;
		int jsum=0;
		double sum=0;
		double bjp=0;
		double bj=1;
		for (int j=m; j>0; j--) {
			double bjm=j*tox*bj-bjp;
			bjp=bj;
			bj=bjm;
			if (fabs(bj)>1e10){
				bj*=1e-10;
				bjp*=1e-10;
				bessj*=1e-10;
				sum*=1e-10;
			}
			if (jsum!=0){sum+=bj;}
			jsum=1-jsum;
			if (j==n) {bessj=bjp;}
		}
		sum=2*sum-bj;
		bessj=bessj/sum;
				 
	}
	
    if ((x<0)&& ((n%2)==1 )) {bessj=-bessj;}
		
	return bessj;
}