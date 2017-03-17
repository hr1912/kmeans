
#include <stdio.h>
#include <math.h>
#include<unistd.h>
#include <stdlib.h>

typedef struct subgraph
{
	int subID;
	int nodeNum;
	int avgLenOrig;
	int avgCvgOrig;
	int avgGCOrig;
	int avgLen;
	int avgCvg;
	int avgGC;
	int clusterID;
}SubGraph;

static SubGraph *subGarray;
static int subGnum;

void getData(char *filename)
{
	FILE *fp;
	subGarray = NULL;

	if ((fp = fopen(filename, "r")) == NULL) {
		printf("Unable to open file: %s, now exit to system\n",filename);
		exit(-1);
	}
	
	char line[1024];
	subGnum = 0;
	// count how namy lines in the file
	while(fgets(line, 1024, fp)){
		subGnum++;
	}
	printf("......There're %d subgraphs in file %s\n",subGnum,filename);
	if(subGnum<1)
		exit(-1);
	
	//allocate array for subgraphs
	subGarray = (SubGraph *)calloc(subGnum,sizeof(SubGraph));
	
	rewind(fp);
	subGnum = 0;
	int minLen=100000;
	int maxLen = 10;
	int minCvg = 100000;
	int maxCvg = 1;
	int minGC = 100;
	int maxGC = 0;

	int subID,nodeNum;
	double avgLen,avgCvg,avgGC;

	while(fgets(line, 1024, fp)){
		sscanf(line+1,"%d node %d avgLen %lf avgCvg %lf avgGC %lf\n",
			&subID,&nodeNum,&avgLen,&avgCvg,&avgGC);
		subGarray[subGnum].subID = subID;
		subGarray[subGnum].nodeNum = nodeNum;
		subGarray[subGnum].avgLen = (int)avgLen;
		subGarray[subGnum].avgCvg = (int)avgCvg;
		subGarray[subGnum].avgGC = (int)avgGC;
		subGarray[subGnum].avgLenOrig = (int)avgLen;
		subGarray[subGnum].avgCvgOrig = (int)avgCvg;
		subGarray[subGnum].avgGCOrig = (int)avgGC;

		minLen = minLen < subGarray[subGnum].avgLen ? minLen : subGarray[subGnum].avgLen;
		minCvg = minCvg < subGarray[subGnum].avgCvg ? minCvg : subGarray[subGnum].avgCvg;
		minGC = minGC < subGarray[subGnum].avgGC ? minGC : subGarray[subGnum].avgGC;
		maxLen = maxLen > subGarray[subGnum].avgLen ? maxLen : subGarray[subGnum].avgLen;
		maxCvg = maxCvg > subGarray[subGnum].avgCvg ? maxCvg : subGarray[subGnum].avgCvg;
		maxGC = maxGC > subGarray[subGnum].avgGC ? maxGC : subGarray[subGnum].avgGC;
		subGnum++;
	}
	printf("*************Len(%d_%d) Cvg(%d_%d) GC(%d_%d)\n",minLen,maxLen,minCvg,maxCvg,minGC,maxGC);
	// normalize the data
	//
	int i;
	for(i=0;i<subGnum;i++)
	{
		if(maxLen-minLen>0){
			subGarray[i].avgLen -= minLen;
			subGarray[i].avgLen *= 100;
			subGarray[i].avgLen /= (maxLen-minLen);
		}
		if(maxCvg-minCvg>0){
			subGarray[i].avgCvg -= minCvg;
			subGarray[i].avgCvg *= 100;
			subGarray[i].avgCvg /= (maxCvg-minCvg);
		}
		if(maxGC-minGC>0){
			subGarray[i].avgGC -= minGC;
			subGarray[i].avgGC *= 100;
			subGarray[i].avgGC /= (maxGC-minGC);
		}
		fprintf(stderr,"C%d\tnode%d\tavgLen%d\tavgCvg%d\tavgGC%d\n",subGarray[i].subID,subGarray[i].nodeNum,
				subGarray[i].avgLen,subGarray[i].avgCvg,subGarray[i].avgGC);
	}
	fclose(fp);
}

static void display_usage()
{
	printf("=============  K-means cluster generator usage  =================\n");
	printf("\t./Kcluster -i fileName [ -k #cluster ]\n\n");
}

long long distance(SubGraph subG,SubGraph cent)
{
	long long dist = 0;
	dist = (subG.avgLen - cent.avgLen)*(subG.avgLen-cent.avgLen);
	dist += (subG.avgCvg - cent.avgCvg)*(subG.avgCvg-cent.avgCvg);
	dist += (subG.avgGC - cent.avgGC)*(subG.avgGC-cent.avgGC);
	
	return dist;
}

void group(SubGraph *seed, int K)
{
	int i,j,index;
	long long  minDis,dist;
	for(i=0;i<subGnum;i++){
		minDis = distance(subGarray[i],seed[0]); index=0;
		for(j=1;j<K;j++){
			dist = distance(subGarray[i],seed[j]);
			if(minDis>dist){
				minDis = dist;
				index = j;
			}
		}
		subGarray[i].clusterID = index;
	}
}

void center(SubGraph *seedNew,int K)
{
	int i;
	for(i=0;i<K;i++){
		seedNew[i].avgLen = 0;
		seedNew[i].avgCvg = 0;
		seedNew[i].avgGC = 0;
		seedNew[i].clusterID = 0;
	}
	int cID;
	for(i=0;i<subGnum;i++){
		cID = subGarray[i].clusterID;
		seedNew[cID].clusterID++;
		seedNew[cID].avgLen += subGarray[i].avgLen;
		seedNew[cID].avgCvg += subGarray[i].avgCvg;
		seedNew[cID].avgGC += subGarray[i].avgGC;
	}
	for(i=0;i<K;i++){
		if(seedNew[i].clusterID>0){
			seedNew[i].avgLen /= seedNew[i].clusterID;
			seedNew[i].avgCvg /= seedNew[i].clusterID;
			seedNew[i].avgGC /= seedNew[i].clusterID;
		}
	}
}
#define  MinDiff 0

unsigned char equal(SubGraph *seed,SubGraph *seedNew,int K)
{
	int i;
	for(i=0;i<K;i++){
		if(seed[i].avgLen-seedNew[i].avgLen>MinDiff)
			return 0;
		if(seed[i].avgCvg-seedNew[i].avgCvg>MinDiff)
			return 0;
		if(seed[i].avgGC-seedNew[i].avgGC>MinDiff)
			return 0;
	}
	return 1;
}

void cpSeed(SubGraph *seed,SubGraph *seedNew,int K)
{
	int i;
	for(i=0;i<K;i++)
		seed[i] = seedNew[i];
}

void cluster(int K)
{
	SubGraph *seed = (SubGraph *)calloc(K,sizeof(SubGraph));
	SubGraph *seedNew = (SubGraph *)calloc(K,sizeof(SubGraph));
	int i;
	//select K seed
	for(i=0;i<K;i++)
		seed[i] = subGarray[i];
	i = 0;
	while(1){
		group(seed,K);
		center(seedNew,K);
		i++;
		printf(".......Completed %d clusterring\n",i);
		if(!equal(seed,seedNew,K))
			cpSeed(seed,seedNew,K);
		else
			break;
	
	}	

	free((void *)seed);
	free((void *)seedNew);
}

static int  cmp_subGraph(const void *e1,const void *e2)
{
	SubGraph *p1,*p2;
	p1 = (SubGraph*)e1;
	p2 = (SubGraph*)e2;
	
	if(p1->clusterID==p2->clusterID){
		if(p1->nodeNum>p2->nodeNum)
			return -1;
		else if(p1->nodeNum<p2->nodeNum)
			return 1;
		else
			return 0;

	}else if(p1->clusterID<p2->clusterID){
		return -1;
	}else
		return +1;
}

void output_cluster()
{
	int i;
	int prevCluster = -1;
	for(i=0;i<subGnum;i++){
		if(subGarray[i].clusterID!=prevCluster){
			prevCluster = subGarray[i].clusterID;
			printf("==============cluster %d ==================\n",subGarray[i].clusterID+1);
		}
		printf("C%d,#node %d, avgLen %d, avgCvg %d, avgGC %d,cluster %d\n",subGarray[i].subID,
			subGarray[i].nodeNum,subGarray[i].avgLenOrig,subGarray[i].avgCvgOrig,subGarray[i].avgGCOrig,subGarray[i].clusterID+1);
		
	}
}

int main(int argc, char *argv[])
{
	int c,K=2;
	unsigned char inputFN = 0;
	char fileName[1024];
        while ((c = getopt(argc, argv, "i:k:")) >= 0) {
                switch (c) {
                case 'k': K = atoi(optarg); break;
                case 'i': inputFN = 1; sscanf(optarg,"%s", fileName);break;
		default:
                	if (inputFN == 0) {       // 
				display_usage();
				exit(-1);
			}
                }
        }

	if (inputFN == 0) {        
		display_usage();
		exit(-1);
	}
	printf("........ to put data in %s into %d clusters  ..........\n",fileName,K);

	getData(fileName);

	cluster(K);

	qsort(subGarray,subGnum,sizeof(SubGraph),cmp_subGraph);

	output_cluster();

	if(subGarray)
		free((void *)subGarray);
}



