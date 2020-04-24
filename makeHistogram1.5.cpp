#include<stdio.h>
#include <stdlib.h>
#include<iostream>
#include<fstream>
#include<cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include <string.h>
using namespace std;
ifstream infile;
//1.5 Bugfix, removed junk code, increased string buffer from 1024 to 4096 ---WARNING if your dihedral line is greater than 4096 bits, you must increase this buffer otherwise floating pt exception will occur or other overflow errors...
//version 1.4 fix bugs in last column numerical value, use stable conversion methods 
//version 1.2 set the number of coordinates, replace natoms*ndihedrals with just ncoords
//bugfix: last column histogram

int main(int argc, char* argv[])
{
if(argc!=2) 
{
cout<<"No input file. Use: makeHistogram [dihedral coordinate file]"<<endl;
return 0;
}
string input="";	
cout<<"Enter the number of coordinates being analyzed"<<endl;
cin>>input;

int ncoords=1,charcount=0,linecount=0,index=0,ind2=0;

ncoords=atoi(input.c_str());
cout<<"You entered "<<ncoords<<" coordinates. Begin histogram"<<endl;

double * max=new double[ncoords];
double * min=new double[ncoords];
double coord=0;
double * coordinates = new double[ncoords]; //dynamic bruh.

for (int i=0; i<(ncoords); i++) 
{
coordinates[i]=0;
min[i]=0;
max[i]=0;

}



char c[4096];
char * s_text;
string number="";
infile.open(argv[1], std::ifstream::in); 	
bool beginAtom=false,record=false;



while(!infile.fail()&&!infile.eof())
{

infile.getline(c,4096);	
s_text=strtok(c," ,{}");//delimiters & set the search at beginning of the line
while(s_text!=NULL) //convert the line to numbers
{
	coord=atof(s_text);	
	if(linecount==0)
	{
		max[index]=coord;
		min[index]=coord;	
	}
	else
	{
		if(max[index]<coord) {max[index]=coord;}
		else if(min[index]>coord) {min[index]=coord;}
	}	
		
	s_text=strtok(NULL," ,{}");
	index++;
}

linecount++;
index=0;	
}	


linecount--;//correction for extra line read on failbit
cout<<"max values for dihedrals"<<endl;
for (int i=0; i<ncoords; i++)
{
	cout<<max[i]<<", ";
//cout<<"res "<<i<<" phi min max: "<<min[2*i]<<' '<<max[2*i]<<endl;
//cout<<"res "<<i<<" psi min max: "<<min[2*i+1]<<' '<<max[2*i+1]<<endl;
}
cout<<endl;
cout<<"minvalues"<<endl;

for (int i=0; i<ncoords; i++)
{
	cout<<min[i]<<", ";
}

cout<<'\n'<<linecount<<" lines analysed, begin histogram for 100 bars at 3.6 degree increments..."<<endl;



//BEGIN HISTOGRAM********************
int totallines=linecount;
infile.clear();//clear failbits/eof
infile.seekg(0,infile.beg);//rewind
linecount=0;
int ** histogram = new int*[100];//natoms columns and 100 rows

for (int j=0; j<100; j++)
{
histogram[j]=new int[ncoords];	
for (int i=0;i<ncoords;i++)
{
histogram[j][i]=0;
}

}//make array 0s

double * hrange=new double[ncoords];

cout<<"histogram differentials: "<<endl;

for(int i=0; i<ncoords; i++)
{
hrange[i]=(max[i]-min[i])/100;//size of the histogram groupings    ********NOT USED NOW
//cout<<hrange[i]<<" "; DEBUG
}

cout<<"percent complete: "<<endl;
index=0;
number="";
while(!infile.fail()&&!infile.eof())
{

printf("%.2i\b\b",(100*linecount/totallines));

infile.getline(c,4096);		
s_text=strtok(c," ,{}");//delimiters & set the search at beginning of the line

while(s_text!=NULL) //convert the line to numbers
{
coordinates[index]=atof(s_text);	
s_text=strtok(NULL," ,{}"); //next number
index++;//next index 	
}	
index=0;
linecount++;	
	for(int i=0; i<ncoords; i++) //coord ary to histogram
	{
		for(int j=0; j<100; j++)
		{
		//ALTERNATIVE    if(coordinates[i]>=(hrange[i]*j)+min[i]&&coordinates[i]<hrange[i]*(j+1)+min[i])//go from smallest value to largest value
		if(coordinates[i]>=(3.6*j)-180&&coordinates[i]<3.6*(j+1)-180)
		{
		histogram[j][i]++;
		break;//for speed
		}
		}
	}

	
}


infile.close();//close the stream


ofstream outputf;
outputf.open("histogramOutput.csv",std::ofstream::trunc);//discard file if it already exists
std::ostringstream oss;


	for(int j=0; j<100; j++)
		{
	for(int i=0; i<ncoords; i++)
	{		
		oss<<histogram[j][i];
		if(i==ncoords-1){oss << "\n";}
		else {oss <<", ";}
		
		}
	}



outputf<<oss.str();

cout<<"\nComplete! Exiting normally"<<endl;

return 0;

}
