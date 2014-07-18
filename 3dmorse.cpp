// 3dmorse
// ==============

// Oleg Devinyak (a,*), Dmytro Havrylyuk (b) and Roman Lesyk (b)

// a)	Department of Pharmaceutical Disciplines, Uzhgorod National University, Uzhgorod 88000, Ukraine
	// e-mail: o.devinyak@gmail.com
// b)	Department of Pharmaceutical, Organic and Bioorganic Chemistry, Danylo Halytsky Lviv National Medical University, Lviv 79010, Ukraine

// http://github.com/devinyak/3dmorse
// --------------

// 3dmorse is an open-source small program to calculate 3D-MoRSE molecular descriptors. 
// Currently it supports only MOPAC2012 output files (*.out) as input. 
// The descriptors produced are 3D-MoRSE weighted with atomic mass, van der Waals volume, electronegativity, polarizability, atomic partial charge and unweighted descriptors. 
// The naming convention is consistent with DRAGON 6 (powerful but only commercially available program for molecular descriptors calculation). 
// That is, the numeration of 3D-MoRSE descriptors starts from 1, so, for example, Mor01u denotes unweighted descriptor with scattering parameter s=0 (since scattering parameter starts from zero), Mor02u denotes descriptor with  s=1 and so on. 
// There is a possibility to obtain a table of 3D-MoRSE terms that correspond to each atomic pair in the molecular structure (for all descriptors at once). 
// This table makes interpretation of 3D-MoRSE descriptors in a QSAR model much easier.

// The typical usage of program is:
// 3dmorse path_to_input_file path_to_output_file <terms flag>
// Terms flag is an optional argument, valid values are 0 (do not return 3D-MoRSE terms) or 1 (return 3D-MoRSE terms). 
// The terms are not returned by default.

// The output is comma separated values file with descriptors in columns. 
// The output of terms has additional fragment "terms" in output file name and its columns are: N - serial number, firstAtom and secondAtom - correspond to atomic pair, s - scattering parameter, Distance - interatomic distance, term - corresponding summand value, weight - weighting scheme.

// The program and its source are distributed under the GNU GPLv3 license.

// For reference or citation use
// Devinyak, O.; Havrylyuk, D.; Lesyk, R. 3D-MoRSE descriptors explained. Submitted to J. Chem. Inf. Model., 2014.


#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iostream>
#include <tchar.h>
#include <string>
//Some constants
const int ATOMS_NUMBER=29;
// Data is taken from CRC Handbook of Chemistry and Physics by D.R. Lide (editor), CRC press 2009-2010, 90th edition. (consistent with DRAGON 6)
const float atomicMasses[ATOMS_NUMBER] = {1.01,6.941,12.01,14.01,16,19,22.991,24.305,28.09,30.97,32.07,35.45,39.098,58.69,63.55,65.39,69.72,74.92,78.96,79.9,107.87,112.41,114.82,118.71,127.6,126.9,200.59,204.38,207.2};
const std::string atomicSymbols[ATOMS_NUMBER] = {"H","Li","C","N","O","F","Na","Mg","Si","P","S","Cl","K","Ni","Cu","Zn","Ga","As","Se","Br","Ag","Cd","In","Sn","Te","I","Hg","Tl","Pb"};
const float atomicVolumes[ATOMS_NUMBER] = {5.42,25.25,20.58,15.6,14.71,13.31,49,21.69,38.79,24.43,24.43,22.45,87.11,18.14,11.49,11.25,27.39,26.52,28.73,26.52,21.31,16.52,30.11,42.8,36.62,32.52,15.6,31.54,34.53};
const float atomicPolarizabilities[ATOMS_NUMBER] = {0.67,24.3,1.76,1.1,0.8,0.56,23.6,10.6,5.38,3.63,2.9,2.18,43.4,6.8,6.1,7.1,8.12,4.31,3.77,3.05,7.2,7.2,10.2,7.7,5.5,5.35,5.7,7.6,6.8};
const float atomicElectronegativities[ATOMS_NUMBER] = {2.59,0.89,2.75,3.19,3.65,4,0.56,1.32,2.14,2.52,2.96,3.48,0.45,1.94,1.98,2.23,2.42,2.82,3.01,3.22,1.83,1.98,2.14,2.3,2.62,2.78,2.2,2.25,2.29};
const int carbonPosition=2;


float atomDistance(float *x, float *y)
{
	float sumofSq=(x[0]-y[0])*(x[0]-y[0])+(x[1]-y[1])*(x[1]-y[1])+(x[2]-y[2])*(x[2]-y[2]);
	float dist=sqrt(sumofSq);
	return(dist);
}


float* prepareWeight(std::string weightType, std::string* Atoms, int nAtoms)
{
	float* weights;
	weights= new float[nAtoms];
	int i=0;

	if (weightType=="mass")
	{
	for (int i=0;i<nAtoms;i++)
	{
	
		for (int j=0;j<ATOMS_NUMBER;j++)
		{
			if (Atoms[i]==atomicSymbols[j])
			{
			weights[i]=atomicMasses[j]/atomicMasses[carbonPosition];
			break;
			}
		}
	}
	}
	if (weightType=="volume")
	{
	for (int i=0;i<nAtoms;i++)
	{
	
		for (int j=0;j<ATOMS_NUMBER;j++)
		{
			if (Atoms[i]==atomicSymbols[j])
			{
			weights[i]=atomicVolumes[j]/atomicVolumes[carbonPosition];
			break;
			}
		}
		}
	}
	if (weightType=="polarizability")
	{
	for (int i=0;i<nAtoms;i++)
	{
	
		for (int j=0;j<ATOMS_NUMBER;j++)
		{
			if (Atoms[i]==atomicSymbols[j])
			{
			weights[i]=atomicPolarizabilities[j]/atomicPolarizabilities[carbonPosition];
			break;
			}
		}
		}
	}
	if (weightType=="electronegativity")
	{
	for (int i=0;i<nAtoms;i++)
	{
	
		for (int j=0;j<ATOMS_NUMBER;j++)
		{
			if (Atoms[i]==atomicSymbols[j])
			{
			weights[i]=atomicElectronegativities[j]/atomicElectronegativities[carbonPosition];
			break;
			}
		}
	
	}
	}
	return(weights);
}

float* prepareWeight(std::string weightType, std::string* Atoms, float * Charge, int nAtoms)
{
	float* weights;
	weights= new float[nAtoms];
	int i=0;
	for (int i=0;i<nAtoms;i++)
	{
		weights[i]=Charge[i];
	}
	return(weights);
}


int main(int argc, char* argv[]) {
	
	    if (argc < 3) 
		{ 
        std::cout << "Usage is 3dmorse input_file output_file <terms flag>\n";
        std::cin.get();
        exit(0);
    } 
		else 
		{ 
		using namespace std; 
        char* myFile;
		char* myOutput;
		int termFlag=0;
		string line;

        
         
                myFile = argv[1];
                myOutput = argv[2];
				if (argc==4) try{termFlag=atoi(argv[3]);} catch(...){std::cout << "Wrong terms flag, should be 0 or 1\n";};
				
				
		std::ifstream moleculeFile(myFile);
		bool isSecondPart=false;
		bool isCoords=false;
		bool isCharges=false;
		bool isFirstLine=true;
		string charAtoms[500];
		float floatCoords[500][3];
		float floatCharges[500];
		int iterator=0;
		int nAtoms;

while(getline(moleculeFile,line))	
{
	if (!isSecondPart){
	if(line.find("----------")!=string::npos){
		isSecondPart=true;};
	continue;
	};

	if (isSecondPart & !isCoords & !isCharges){
		if(line.find("CARTESIAN COORDINATES")!=string::npos){
			isCoords=true;};
			if(line.find("ATOM NO")!=string::npos){
			isCharges=true;}
			continue;
		}

		if (isSecondPart & isCoords & !isCharges){
			int N;
			if (line.length()>5){
		stringstream ss(line);
		ss >> N >> charAtoms[iterator] >> floatCoords[iterator][0] >> floatCoords[iterator][1] >> floatCoords[iterator][2];
		iterator++;
			}
			else {
				if(isFirstLine) {isFirstLine=false;continue;}
				isCoords=false; nAtoms=iterator; iterator=0; 			}
			
		continue;
		}


		if (isSecondPart & !isCoords & isCharges){
			int N;
			char C;
			if (line.find("DIPOLE")==string::npos)
			{
		stringstream ss(line);
		ss >> N >> C >> floatCharges[iterator];
		iterator++;
			}else {
				isCharges=false;
				isSecondPart=false;
			}
			
		continue;
						
        }
				
}
//********************************************************
//END-OF-INPUT; START CALCULATIONS

float morU[32]={0};
float morM[32]={0};
float morV[32]={0};
float morP[32]={0};
float morE[32]={0};
float morC[32]={0};

	float* weightM=prepareWeight("mass",charAtoms,nAtoms);
	float* weightV=prepareWeight("volume",charAtoms,nAtoms);
	float* weightP=prepareWeight("polarizability",charAtoms,nAtoms);
	float* weightE=prepareWeight("electronegativity",charAtoms,nAtoms);
	float* weightC=prepareWeight("pcharge",charAtoms,floatCharges,nAtoms);
	float thisDistance;
	float member;
	int rownumerator=0;
	ofstream table;
	if(termFlag==1)	{
		
		std::string helpbuffer=string(myOutput);
		helpbuffer.insert(helpbuffer.length()-4,"terms");
		table.open(helpbuffer);
		table << "N" << "," << "firstAtom" << "," << "secondAtom" << ","<< "s" << "," << "Distance" << "," << "term" << "," << "weight" << "\n" ;
					}

	for (int i=0;i<(nAtoms-1);i++)
	for (int j=i+1;j<nAtoms;j++)
	for (int s=0;s<32;s++){

		thisDistance=atomDistance(floatCoords[i],floatCoords[j]);
		if(s==0) member=1; else 
		{	
				member=sin(s*thisDistance)/(s*thisDistance);
		}
		morU[s]+=member;
		morM[s]+=weightM[i]*weightM[j]*member;
		morV[s]+=weightV[i]*weightV[j]*member;	
		morP[s]+=weightP[i]*weightP[j]*member;
		morE[s]+=weightE[i]*weightE[j]*member;
		morC[s]+=weightC[i]*weightC[j]*member;

		if(termFlag==1)	{
		table << rownumerator++ << "," << charAtoms[i] << "," << charAtoms[j] << "," << s << "," << thisDistance << "," << member << "," << "unweighted" << "\n" ;
		table << rownumerator++ << "," << charAtoms[i] << "," << charAtoms[j] << ","<< s << "," << thisDistance << "," << weightM[i]*weightM[j]*member << "," << "mass" << "\n" ;
		table << rownumerator++ << "," << charAtoms[i] << "," << charAtoms[j] << ","<< s << "," << thisDistance << "," << weightV[i]*weightV[j]*member << "," << "volume" << "\n" ;
		table << rownumerator++ << "," << charAtoms[i] << "," << charAtoms[j] << ","<< s << "," << thisDistance << "," << weightP[i]*weightP[j]*member << "," << "polarizability" << "\n" ;
		table << rownumerator++ << "," << charAtoms[i] << "," << charAtoms[j] << ","<< s << "," << thisDistance << "," << weightE[i]*weightE[j]*member << "," << "electronegativity" << "\n" ;
		table << rownumerator++ << "," << charAtoms[i] << "," << charAtoms[j] << ","<< s << "," << thisDistance << "," << weightC[i]*weightC[j]*member << "," << "charge" << "\n" ;
						}
	}
	
	
	if(termFlag==1) table.close();
	delete [] weightM;
	delete [] weightV;
	delete [] weightP;
	delete [] weightE;
	delete [] weightC;
// -----OUTPUT-BLOCK-------	
      ofstream output;
		output.open(myOutput);
		
		for (int i=0;i<32;i++)
		{output << "Mor";
		if(i<9) output << "0";
		output << i+1 << "u" << ",";
		}
		for (int i=0;i<32;i++)
		{output << "Mor";
		if(i<9) output << "0";
		output << i+1 << "m" << ",";
		}
		for (int i=0;i<32;i++)
		{output << "Mor";
		if(i<9) output << "0";
		output << i+1 << "v" << ",";
		}
		for (int i=0;i<32;i++)
		{output << "Mor";
		if(i<9) output << "0";
		output << i+1 << "p" << ",";
		}
		for (int i=0;i<32;i++)
		{output << "Mor";
		if(i<9) output << "0";
		output << i+1 << "e" << ",";
		}
		for (int i=0;i<32;i++)
		{output << "Mor";
		if(i<9) output << "0";
		if (i!=31) output << i+1 << "c" << ",";
		else output << i+1 << "c" << "\n";
		}
		
		for (int i=0;i<32;i++)
		{output << morU[i] << ",";
		}
		for (int i=0;i<32;i++)
		{output << morM[i] << ",";
		}
		for (int i=0;i<32;i++)
		{output << morV[i] << ",";
		}
		for (int i=0;i<32;i++)
		{output << morP[i] << ",";
		}
		for (int i=0;i<32;i++)
		{output << morE[i] << ",";
		}
		for (int i=0;i<32;i++)
		{output << morC[i];
		if (i!=31) output << ","; else output << "\n";
		}
		
		output.close();
		
}
	return 0;

		}
