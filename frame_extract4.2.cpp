#include<stdio.h>
#include<iostream>
#include<fstream>
#include<cmath>
#include <iomanip>
#include <string>
#include <sstream>
#include <cstring>
#include <stdlib.h>


using namespace std;

ifstream traj,deviation,probfile;
ofstream outputf,topclustersout;
//input arguments: frame_extract1 [traj_dihedrals] [deviation file] [probfile] [epsilon] [#of top clusters] [% (optional)]
//version 2: realized that the correlation function is required to choose structures. Skip top coordinate file and just compare each frame in the trajectory to 
//the probability file with deviation. Record the indexes at each frame and cluster. the probability score will be based on the number of frames corresponding to the specific Pcoord combination.
//this will take into account the correlation and drastically reduce computational time.
//2.1 Added a third probability/std
//2.2 print traj sample eery 1000 frames, and changed natoms and ncoords to just ncoords -- for user input too.
//2.3 variable number of top clusters to display and traj fraction accounted for.
//2.4 peak probability cut off -- using the probability of the peak, dont bother clustering coordinate if the peak pvalue is less than 5% ~ arbitrary for now
//also just read the lines in the traj to get frame numbers. 
//2.5 added a cluster relaxation method
//2.6 read the numframes and ncoords from the p input
//2.7 save the output of the top clusers
//2.8 make the relaxation default and sigma values, topclusters to extract arguments. If a % follows the top clusters, extract until that %traj occurs
//2.9 made the relaxation of the populations a user enterable value
//3.1 -- Changed comparisons between FT, changed output file format for top coords, changed some relaxation criteria.
//3.2 added option to exclude clusters containing an unmatched coordinate from the top clusters 
//3.3 write the transformed trajectory to a file transformed_traj.dat
//4.1 added markov state model calculation
//4.2 added the distribution probabilities 
double average(double * array,int size)
{
	double total=0;
	for(int k=0; k<size; k++)
	{
		total+=array[k]/size;
		
	}
	
	return total;
	
}

void getattributes(ifstream & pfile, ifstream & dfile, int & lcount, int & ccount)
{
	char c[1024];
	
	while(!pfile.fail()&&!pfile.eof())
{
	pfile.getline(c,1024);
	if(c[0]!='\0') lcount++;
}
pfile.clear(); //clear flags
pfile.seekg(0,ios::beg);//rewind streambuf


	while(!dfile.fail()&&!dfile.eof())
{
	dfile.getline(c,1024);
	if(c[0]!='\0') ccount++;
}

dfile.clear(); //clear flags
dfile.seekg(0,ios::beg);//rewind streambuf



}

int main(int argc, char* argv[])
{
double epsilon=1;//scale the stdev error
double pcutoff=0.01; //probability cutoff for a peak. If its less than 10%, dont bother clustering it.
string input="";	
int ncoords=0;	
int numframes=0,index1=0;

	if(argc<6||argc>7)
	{
		cout<<"Frame_extract [coordinate trajectory] [deviation file] [probfile] [epsilon] [topclusters] [%to extract (optional)]"<<endl;
		return 1;
	}

traj.open(argv[1], std::ifstream::in);	//open the trajectory file of dihedral angles	
deviation.open(argv[2], std::ifstream::in);	
probfile.open(argv[3], std::ifstream::in);	
	
getattributes(traj,probfile,numframes,ncoords);

printf("Performing CATS using %i lines with %i coordinates\n\n",numframes,ncoords);

//cout<<"Enter epsillon value"<<endl;
//cin>>input;
epsilon=atof(argv[4]);
input="";

	
//the parray has 3 p values per coordinate and 2 avg values per coordinate. 4*number of coordinates total
double * dev_ary= new double[3*ncoords];
double * prob_ary=new double[ncoords*6];

//double * topcoords_ary=new double[ntopcoords*ncoords+ntopcoords]; 

for(int i=0; i<6*ncoords;i++) //initialize to 0... watch out if the topcoords are less than 4!
{
	if(i<3*ncoords) dev_ary[i]=0;
	prob_ary[i]=0;
	
}





char c[1024];
bool record=false; 
string number="";

 
 //retrieve the values from the deviation file
 cout<<"retrieving deviations..."<<endl;
 
while(!deviation.fail()&&!deviation.eof())
{
	deviation.getline(c,1024);
	number="";
	for (int i=0; i<1024; i++)
	{
	if(c[i]=='\0')
	{
		break;
	}
	else if ((c[i]-'0')>-4&&(c[i]-'0')<10)
	{
		record=true;
		number+=c[i]; //assign char to temp array
		

	}
	else if ((c[i]==' '||c[i]==',')&& record) //at the end of a number
	{
		
		record=false;
		dev_ary[index1]=atof(number.c_str());
		printf("%f ",dev_ary[index1]);
		number=""; //clear string
		index1++;	
	}
	
	}
	printf("\n");
	
}

cout<<"retrieving prob values..."<<endl;

//reset values
int index=0;
record=false;
index1=0;
number="";
bool * skip_coords=new bool[ncoords];
for(int i=0; i<ncoords; i++)//based on the # sign in the p file. This array will allow for coordinates to not control clustering. 
{
	skip_coords[i]=false;
}

//get probs

while(!probfile.fail()&&!probfile.eof())
{
	probfile.getline(c,1024);
	if(c[0]=='#') skip_coords[index]=true;
	
	for (int i=0; i<1024; i++)
	{
	if(c[i]=='\0'&&i>0)//if the character is a null at the end of the line
	{
		record=false;
		prob_ary[index1]=atof(number.c_str());
		printf("%f ",prob_ary[index1]);
		number="";
		break;
	}
	else if ((c[i]-'0')>-4&&(c[i]-'0')<10)
	{
		record=true;
		number+=c[i]; //assign char to temp array

	}
	else if ((c[i]==' ' || c[i]==',') && record) //at the end of a number
	{
		
		record=false;
		prob_ary[index1]=atof(number.c_str());
		printf("%f ",prob_ary[index1]);
		number=""; //clear string
		index1++;
		
	}
	else if(c[i]=='\0') break;  
	
	}
	printf("\n");
	index++;//indexer for the skip coords array
}

record=false;
index1=0;
number="";
index=0;
//get the indexes of the top coords compared to the stddev and the prob file

//int *topcrd_indexer=NULL;
//topcrd_indexer=new int[ncoords*ntopcoords];
//double *topcoords_devs=new double[ncoords*ntopcoords]; //do a 1 to 1 relation between devs and the top coord


int traj_framenumber=0;
bool frame_captured=false;
ofstream outputf;
outputf.open("frame_extract_results.csv",std::ofstream::trunc);//discard file if it already exists
std::ostringstream oss;

//need dynamic allocation of memory with these large arrays.

int * matched_frames=NULL;
matched_frames=new int[numframes];//the line number corresponds to the cluster number, and the line contains the frames matched

int * matched_frames_inexer=new int[numframes]; //keeps track of how many clusters in each line... each element corresponds to a line, and the number is the number of frames inside the cluster

double *lowest_difference=new double[numframes]; //this records the difference value from the average -- [frame number, off center distance]

double avgstd=0;
int * traj_ary=NULL;
traj_ary=new int[numframes*ncoords];

for(int i=0; i<numframes*ncoords; i++) //assign zeros
{
		if(i<numframes) matched_frames[i]=0;
		traj_ary[i]=0;
}


cout<<"Begining comparison..."<<endl;

printf("Reading trajectory     ");



double *traj_line=new double[ncoords];

while(!traj.fail()&&!traj.eof()) //get the trajectory from file. Don't even bother saving it to memory, just compare the prob file
{
	
	printf("\b\b\b\b%3i%%",100*traj_framenumber/numframes);
	
	traj.getline(c,1024); //get each line of the traj

	for (int i=0; i<1024; i++)
	{
		
	if(c[i]=='\0'&&i>0)//if the character null and its not the last line
	{
		record=false;
		traj_line[index1]=atof(number.c_str());
		//printf("%f ",prob_ary[index1]);
		number="";
		break;
	}
	
	else if ((c[i]-'0')>-4&&(c[i]-'0')<10)
	{
		record=true;
		number+=c[i]; //assign char to string

	}
	else if ((c[i]==' '||c[i]==',') && record) //at the end of a number
	{
		
		record=false;
		traj_line[index1]=atof(number.c_str());//get the line
		//traj_ary[traj_framenumber*ncoords+index1]=atof(number.c_str());
		//printf("%f ",traj_ary[traj_framenumber*ncoords+index1]);
		number=""; //clear string
		index1++;
		
	}
	else if(c[i]=='\0'&&i==0) break; //eof.
	
	
	} //end of the line read
	lowest_difference[traj_framenumber]=0; //easier than doing another loop to make zeros
	
	for (int i=0; i<ncoords; i++) //translate the line into the probfile index 0 or 1 aka first peak or second peak
	{
		//if the probability of the peak is less than 5%, just ignore it. 
		if(180-abs(abs(traj_line[i]-prob_ary[6*i])-180)-epsilon*dev_ary[3*i]<0&&prob_ary[6*i+1]>pcutoff) //see if coordinate 1 is within tolerance of first peak if not then it must be in the 2nd peak
		{
			traj_ary[ncoords*traj_framenumber +i]=0;
			lowest_difference[traj_framenumber]+=(180-abs(abs(traj_line[i]-prob_ary[6*i])-180))/ncoords; //the average of the differences
		}
		else if(180-abs(abs(traj_line[i]-prob_ary[6*i+2])-180)-epsilon*dev_ary[3*i+1]<0&&prob_ary[6*i+3]>pcutoff)
		{
			traj_ary[ncoords*traj_framenumber +i]=1;
			lowest_difference[traj_framenumber]+=(180-abs(abs(traj_line[i]-prob_ary[6*i+2])-180))/ncoords;
		}
		else if(180-abs(abs(traj_line[i]-prob_ary[6*i+4])-180)-epsilon*dev_ary[3*i+2]<0&&prob_ary[6*i+5]>pcutoff)
		{
			traj_ary[ncoords*traj_framenumber +i]=2;
			lowest_difference[traj_framenumber]+=(180-abs(abs(traj_line[i]-prob_ary[6*i+4])-180))/ncoords;
		}
		else
		{
		//cout<<"Frame "<<traj_framenumber<<" does not match any probs on coordinate "<<i<<endl;//might actually show up for 3 peak distributions
		traj_ary[ncoords*traj_framenumber +i]=3; //not matched
		lowest_difference[traj_framenumber]+=180;//basically make it larger than anything.
		}
		//lowest_difference[2*traj_framenumber]=traj_framenumber;// frame number is the first column
		//printf("%i,",traj_ary[ncoords*traj_framenumber +i]);
	}
	//printf("\n");
	
	
	index1=0;//reset index1
	traj_framenumber++;//increase frame number

	
} //end of traj input

ofstream outputtraj;
outputtraj.open("transformed_trajectory.dat",std::ofstream::trunc);//discard file if it already exists
std::ostringstream traj_stream;
for(int i=0; i<numframes;i++)
{
	for(int j=0; j<ncoords; j++)
	{
		traj_stream<<traj_ary[ncoords*i+j]<<' ';
		
	}
	traj_stream<<'\n';
	
}
outputtraj<<traj_stream.str();
outputtraj.close();



for(int i=0; i<numframes; i+=1000)//print out traj every thousand frames
{
	for(int j=0; j<ncoords; j++)//print out traj every thousand frames
	{
	
	 printf("%i, ",traj_ary[ncoords*i +j]);
	
	}
	printf("\n");
}


printf("\nCounted %i frames. Begin clustering\n\n",traj_framenumber);
//numframes=traj_framenumber-1; //in case the user count was off

index1=1; //reset counter of matches

//printf("  0%%");

//cluster the traj_ary results 

int * clusters =new int[numframes*ncoords]; //too big of an array but safe.
int nclusters=1, n_ignorable_coords=round(ncoords*.1), n_unclustered=0, n_lowestpop=numframes*.0001;//from the first frame assignment, absorb cluster if unclusterable count under 10%, count the number of unclusterable coords, lowest population to absorb cluster.
bool unclustered=false;
int * largest_deviations=new int[ncoords]; //record what residue/coordinates produce the most deviations in the trajectory and create new clusters
matched_frames[0]=0;//cluster 0 has frame 0 in it

printf("Enter the population cutoff for clusters in relaxation\n");
cin>>n_lowestpop;
printf("Enter the number of ignorable coordinates\n");
cin>>n_ignorable_coords;
printf("The relaxation population cutoff will be %i...and the number of ignorable coords will be %i Beginning relaxation...\n",n_lowestpop,n_ignorable_coords);

	for (int i=0; i<ncoords; i++) //the first frame is the first cluster...
	{
		clusters[i]=traj_ary[i];
		largest_deviations[i]=0;
		//matched_frames_inexer[0]=1;//one frame in the cluster
		
	}
	matched_frames[0]=0; //pair cluster number,traj frame number.
	printf("Completed:         "); 
	
	for (int i=1; i<numframes; i++) //starting from the next line
	{
		printf("\b\b\b\b%3i%%",100*i/numframes);
		for(int k=0; k<nclusters; k++) 
		{						
			n_unclustered=0; //reset 
			for(int j=0;j<ncoords;j++)//compare each coordinate
			{
				if(skip_coords[j]) frame_captured=true;
				
				//see what happens if we dont let the unclustered ones slide...
				/*
				else if(traj_ary[ncoords*i+j]==3) 
					{//this does not cluster that frame.
						//unclustered=true;
						//frame_captured=false;
						//break;
						frame_captured=true; //this ignores that unclusterable coordinate, up to a point
						n_unclustered++;
						if(n_unclustered>n_ignorable_coords)//only uncluster it if we need to ignore 10% of the coords. 
						{
						unclustered=true;
						frame_captured=false;
						break;	
						}
					}
					*/
				
				else
				{
					if(clusters[ncoords*k+j]==traj_ary[ncoords*i+j]) frame_captured=true;
					else
					{
						//largest_deviations[j]++;
						frame_captured=false;
						break;
					}
				}
				
			}
				
			if(frame_captured) //if the frame is already represented by a cluster
			{
			
			matched_frames[i]=k;//traj frame i has cluster k associated with it.
			frame_captured=false;
			index1++;//increase the match count
			break; //stop going through the k clusters.
			}
			
			/*
			else if(unclustered)//skip this frame and go to the next.
			{
				unclustered=false; 
				matched_frames[i]=-1;//set the cluster to -1
				break;
			}
			*/
			
			
			else if(!frame_captured&&k==(nclusters-1))//looked through all of the clusters with no match, create new cluster.
			{
				
				for (int j=0;j<ncoords;j++)//assign each coordinate
				{	
				clusters[ncoords*nclusters+j]=traj_ary[ncoords*i+j];
				}
				matched_frames[i]=nclusters;//record the new cluster/frame
				index1++; //increase match count
				nclusters++;//keep count.
				break;
			}

			
		}

		
	}//end clustering 
	
	n_unclustered=0;
	

	cout<<"\n\nFinished with "<<nclusters<<" clusters... refining clusters based on "<<n_ignorable_coords<<" ignorable coordinates in clusters with populations less than "<<(int)n_lowestpop<<endl;
	int * cluster_counter=new int[nclusters];//how many frames are in this cluster
	int * cluster_addons =new int[3*nclusters];//add on to a cluster number a maximum of nlowestpop frames. [old cluster#, new cluster#, frame to add]
	int ** cluster_frame_tracker = new int*[nclusters];//keep track of the cluster frame numbers in each cluster
	
	for(int i=0; i<nclusters; i++)
	{
		cluster_counter[i]=0;
	}
	
	
	
	for(int i=0; i<nclusters; i++)//number of frames in each cluster
	{

		for(int j=0; j<numframes; j++)
		{
			if(matched_frames[j]==i) 
			{
				cluster_counter[i]++;
			}
		}
	}
	
	
	
	index=0;
	int index2=0;
	int ccounter=0;
	//***********BEING REFINING CLUSTERS BASED ON RELAXATION*****************
	
	

	cout<<"Beginngin relaxation of clusters..."<<endl;
	printf("completed:     ");
for(int i=0; i<nclusters; i++)
{
		printf("\b\b\b\b%3i%%",(int)(100*i/nclusters));	
		cluster_frame_tracker[i]=new int[cluster_counter[i]+1];//track which frames are in what cluster. oversized matrix.
		
		for(int j=0; j<numframes; j++)
		{
			
			if(matched_frames[j]==i) 
			{
			cluster_frame_tracker[i][index]=j;//cluster i has the following frames attached to it.
			index++;
			}
			
		}
		index=0;//reset
		
if(cluster_counter[i]<n_lowestpop)//if the population is too low in the cluster (~20 frames) then try absorbing them into other clusters. 
{
	ccounter=cluster_counter[i];
	for(int l=0; l<ccounter; l++)//for each frame in the low populated cluster try to match them with relaxation to the existing clusters
	{
		for(int k=0; k<nclusters; k++) //for each cluster, not it's native cluster.
		{		
			
			n_unclustered=0; //reset 
			for(int j=0;j<ncoords;j++)//compare each coordinate
			{
				if(skip_coords[j]) frame_captured=true; //if theres a # sign ignoring corrdinate
				
				else
				{
					if(clusters[ncoords*k+j]==traj_ary[ncoords*cluster_frame_tracker[i][l]+j]) frame_captured=true;
					else
					{
						//largest_deviations[j]++;
						n_unclustered++;
						if(n_unclustered>n_ignorable_coords)//only uncluster it if we need to ignore 10% of the coords. 
						{
						frame_captured=false;
						break;
						}
					}
				}
				
			}
				
			if(frame_captured) //if the frame represented by a cluster
			{
				
			//cluster_addons[3*index2]=i; //old cluster i
			//cluster_addons[3*index2+1]=k; //new cluster k &
		//	cluster_addons[3*index2+2]=cluster_frame_tracker[i][l]; //this frame added to it.
			matched_frames[cluster_frame_tracker[i][l]]=k;//change the frame number to the cluster k.			
			
			
			cluster_counter[i]--;//take a frame away from cluster counter
			cluster_counter[k]++;

			frame_captured=false;
			//index2++;//increase the index
			break; //stop going through the k clusters.
			}
			

			if(k==i-1) k++; //this will make the loop skip it's own cluster so it cant match to itself.
		}//end for each cluster

		
	}//end for each frame in the cluster
			
				
				
				
			
			
			
			
	}//end if the population is too low
			
	}//end for each cluster.
			

	
		input="";
	//sorting clusters
	int ntopclusters=1000; //default
	int nclusters2=0;//count the number of non zero clusters after refinement. 
	bool topclusterpercent=false;
	double topclusterpercent_sum=0,topclusterspercent_val=0;
	if(argv[argc-1]!= "%") //if the % argument is left out
	{
	ntopclusters=atoi(argv[argc-1]);
	}
	else if(argv[argc-1]=="%")//if we add up to a percent then sum the top clusters until that percent is reached 
	{
		topclusterspercent_val=atof(argv[argc-2]);//convert the percent to a float
		cout<<"top "<<topclusterspercent_val<<'%'<<endl;
		topclusterpercent=true;
		ntopclusters=1000; //go 1000 frames and stop when the % is reached. 
	}
	else
	{
	cout<<"Enter the number of top clusters to show:"<<endl;
	cin>>ntopclusters;	

	}
	
	cout<<"Exclude top ranked clusters with unsortable coordinates in them? 0=no 1=yes"<<endl;
	bool exclude_clusters=false,bad_coordinate=false;
	cin>>exclude_clusters;
	if(exclude_clusters) cout<<"Excluding clusters with unsortable codes in them"<<endl;
	else if(!exclude_clusters) cout<<"Not excluding any clusters"<<endl;
	else cout<<"something went wrong..."<<endl;
	
	double * cluster_distprob=new double[nclusters];//distribution based probability for each cluster. 
	int *topclusters=new int[2*ntopclusters]; //top 10 clusters, [clusternumber,size]
		for(int j=0; j<2*ntopclusters; j++) 
			{
				topclusters[j]=0; //set to zeros
			}
		
		//assign defaults:
		topclusters[0]=0;
		topclusters[1]=cluster_counter[0];
		
		printf("sort by distribution probability? (y/n) \n");
		cin>>input;
		if(input=="y")
		{
			cout<<"sorting by distribution probability... "<<endl;
		for(int i=0; i<nclusters; i++)//calculate all dist probabilities. 
		{
			cluster_distprob[i]=1; //initialize
			
			for(int j=0; j<ncoords; j++)
			{
				if(!skip_coords[j]) cluster_distprob[i]=cluster_distprob[i]*prob_ary[6*j+2*clusters[ncoords*i+j]+1];
				
				
			}
			
			
		}
			
		
		for(int i=1; i<nclusters; i++)//compare one cluster to all the others 
		{
			
			
			if(cluster_counter[i]!=0) nclusters2++;//count the number of non zero clusters after relaxation.
								
					if(exclude_clusters)
					{
						bad_coordinate=false;
						for(int k=0; k<ncoords;k++)
						{
							if(!skip_coords[k]&&traj_ary[ncoords*cluster_frame_tracker[i][0]+k]==3) 
							{
								bad_coordinate=true;
								break;
							}
							
						}
						
					}
					
					
		if((!bad_coordinate&&exclude_clusters)||!exclude_clusters)		
		{
			for(int j=0; j<ntopclusters; j++) 
			{
				
				if(cluster_distprob[i]>topclusters[2*j+1])//population bigger than on the list
				{


					if(j==ntopclusters-1) //at the bottom, we cant switch
					{
						topclusters[2*j+1]=cluster_counter[i];//pop
						topclusters[2*j]=i;//cluster number
					}
					else
					{
						for(int k=ntopclusters-1; k>j; k--)//move the list down by one
						{
							if(k!=0)
							{
							topclusters[2*k+1]=topclusters[2*k-1];//pop
							topclusters[2*k]=topclusters[2*k-2];//cluster number	
							}
						}
						topclusters[2*j+1]=cluster_counter[i];//pop
						topclusters[2*j]=i;//cluster number	
						break;
					}
					
					
				}
				
				
			}
		}
		//printf("C%i:%i, ",i,cluster_counter[i]);
		
		}
			
		}

		else
		{
			cout<<"sorting by population"<<endl;
			
		
		for(int i=1; i<nclusters; i++)//compare one cluster to all the others 
		{
			
			
			if(cluster_counter[i]!=0) nclusters2++;//count the number of non zero clusters after relaxation.
								
					if(exclude_clusters)
					{
						bad_coordinate=false;
						for(int k=0; k<ncoords;k++)
						{
							if(!skip_coords[k]&&traj_ary[ncoords*cluster_frame_tracker[i][0]+k]==3) 
							{
								bad_coordinate=true;
								break;
							}
							
						}
						
					}
					
					
		if((!bad_coordinate&&exclude_clusters)||!exclude_clusters)		
		{
			for(int j=0; j<ntopclusters; j++) 
			{
				
				if(cluster_counter[i]>topclusters[2*j+1])//population bigger than on the list
				{


					if(j==ntopclusters-1) //at the bottom, we cant switch
					{
						topclusters[2*j+1]=cluster_counter[i];//pop
						topclusters[2*j]=i;//cluster number
					}
					else
					{
						for(int k=ntopclusters-1; k>j; k--)//move the list down by one
						{
							if(k!=0)
							{
							topclusters[2*k+1]=topclusters[2*k-1];//pop
							topclusters[2*k]=topclusters[2*k-2];//cluster number	
							}
						}
						topclusters[2*j+1]=cluster_counter[i];//pop
						topclusters[2*j]=i;//cluster number	
						break;
					}
					
					
				}
				
				
			}
		}
		//printf("C%i:%i, ",i,cluster_counter[i]);
		
		}
		
		}
			
		for(int i=1; i<numframes; i++)//scroll through the trajectory and find what residues/coordinates change alot.
		{
			for(int j=0;j<ncoords; j++)
			{
				if(traj_ary[ncoords*i+j]!=traj_ary[ncoords*(i-1)+j]&&!skip_coords[j]) largest_deviations[j]++;
				
			}
			
			
		}
		
		if(topclusterpercent)
		{
			
			for (int i=0; i<nclusters; i++)
			{
				topclusterpercent_sum+=(double)topclusters[2*i+1]/numframes;
				if(topclusterpercent_sum>=(double)topclusterspercent_val/100)
				{
					nclusters2=i;
					ntopclusters=i;
				}
			}
				
		}
		
		
		int * topcluster_frames=new int[ntopclusters];//record the frames with the lowest avg deviation from avg value in prob table.
		double lowerbound=0;//default initial.
		
		for(int i=0; i<ntopclusters; i++)
		{
			for(int j=0;j<numframes; j++)
			{
				if(topclusters[2*i]==matched_frames[j])		
				{
					if(lowerbound==0)//default initial
					{
					topcluster_frames[i]=j;
					lowerbound=lowest_difference[j];
					}
					else if(lowest_difference[j]<lowerbound)
					{
					topcluster_frames[i]=j;
					lowerbound=lowest_difference[j];	
					}
				}
			}
			lowerbound=0;//reset
			
		}
		
//create a map of what each top cluster's transform is based on the best frame		
cout<<"The top "<<ntopclusters<<" clusters have the form:"<<endl;

for(int i=0;i<ntopclusters;i++)
{
printf("cluster %i: ",i);	
for(int j=0;j<ncoords;j++)
{
	if(skip_coords[j]) printf("x ");
	else
	{
printf("%i, ", traj_ary[ncoords*topcluster_frames[i]+j]);
	}
	
}
printf("\n"); 
	
}

		

cout<<"Total of "<<index1<<" matches with epsilon="<<epsilon<<" accounting for "<<100*index1/numframes<<"% of the trajectory in "<<nclusters2<<" with relaxation and "<<nclusters<<" clusters without relaxation"<<endl;

printf("The top %i clusters/population:\n",ntopclusters);
double total=0,pop=0;//totals the %s of the clusters

topclustersout.open("topclusters.out",std::ofstream::trunc);




for(int i=0; i<ntopclusters; i++)
{
	probability=1;//reset
	
	for(int j=0; j<ncoords; j++)
	{
		if(!skip_coords[j])
		{
			probability=probability*prob_ary[6*j+2*clusters[ncoords*topclusters[2*i]+j]+1]; //use cluster array to take peak probability for conditional propbability. 
		}
	}
	
	
pop=(double)100*topclusters[2*i+1]/numframes;//percent population of a cluster.

topclustersout<<"Cluster "<<topclusters[2*i]<<" Population: "<<topclusters[2*i+1]<<" -> "<<pop<<'%'<<" distribution prob: "<<probability<<" bestframe: "<<topcluster_frames[i]<<endl;
printf("Cluster %i: %i - %.3f%% of total trajectory & Pdist= %f . Best frame: %i  \n",topclusters[2*i],topclusters[2*i+1],pop,probability,topcluster_frames[i]);
total+=	pop;
}

cout<<"Totaling "<<total<<"% of the trajectory"<<endl;
topclustersout<<"Totaling "<<total<<"% of the trajectory"<<endl;
cout<<"\n\nThe largest coordinate deviations:"<<endl;

for(int i=0; i<ncoords; i++)//print out the deviations of clusters
{
	printf("%i, ",largest_deviations[i]);
}

for (int i=0; i<nclusters; i++)
{
	
	oss<<"Cluster "<<i<<": ";
	for(int j=0; j<numframes; j++)
	{
		if(matched_frames[j]==i)
		{
	oss<<j<<' ';
		}
		
	}
	
	oss<<'\n';
	
}


outputf<<oss.str();
topclustersout.close();
//delete [] matched_frames;
//delete [] traj_ary;
//delete [] dev_ary;
//delete [] topcrd_indexer;
//delete [] prob_ary;

//garbage disposal
outputf.close();



cout<<"Generate MSM? 1=yes 0=no"<<endl;
cin>>input;
if(input=="1")
{
	int lagtime=0;
	oss.str("");
	outputf.open("msm.dat",std::ofstream::trunc);
	cout<<"Enter max lagtime"<<endl;
	cin>>lagtime;
	cout<<"Using lagtime of "<<lagtime<<" seconds"<<endl;
	double * msmm=new double[lagtime*nclusters*nclusters];
	
	for(int i=0;i<lagtime;i++)//lagtime starts at 1
	{
		for(int j=0;j+i+1<numframes;j++)
		{	
		msmm[i*nclusters*nclusters+nclusters*matched_frames[j]+matched_frames[j+i+1]]++;
			
			
		}
		
		
	}
	index=0;
	index1=0;
	//print out a matrix of the top clusters only
	for(int i=0;i<lagtime;i++)//lagtime starts at 1
	{
		oss<<"T="<<i+1<<endl;
		
		for(int j=0;j<ntopclusters;j++)
		{
			oss<<'c'<<topclusters[2*j]<<" ";
		for(int k=0;k<ntopclusters;k++)
		{	
			oss<<msmm[i*nclusters*nclusters+topclusters[2*j]*nclusters+topclusters[2*k]]<<' ';
			
		}
		oss<<'\n';
		}
		
printf("Completed T=%i\n",i);
	}		
	outputf<<oss.str();
	outputf.close();
	
}

cout<<"\nComplete"<<endl;


return 0;	
}
