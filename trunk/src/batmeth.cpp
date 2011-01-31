#include "config.h"
#include <string>
#include <stdio.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <signal.h> 
//#include <dvec.h>
#include "zlib.h"
#include <queue>
#include <ctype.h>
#include "batlib.h"
#include <pthread.h>
#include "map.h"
#include "commandline.h"
#include <string>
#include <getopt.h>
using namespace std;

const int FILEBUFSIZE=60000;//Space for batman to write records..
const int MAX_REC_SIZE = 2000;//maximum size a record of batman will occupy..
const int INITIAL_PROGRESS_READS =1000;//progress bar initial count..
FILE* Output_File;
int MAXHITS = 50;
int MAX_MISMATCHES=2;

int getMismatch(char* s2);
inline void Print_List(bool & Header_Printed,int cntHit,int firstmismatch,int lastMismatch,string* tmpList,char* Header);
bool Get_Records(char* File_Buffer, string*  tmpList,int & firstmismatch,int & lastMismatch,int & cntHit,char Count);
void Build_Names(const char* Genome_Name,FMFILES & FM);//build FM index file names...
void Parse_Command_line(int argc, char* argv[],char* & GENOME,char* & CT_GENOME,char* & GA_GENOME,char* & INPUTFILE,char* & OUTPUTFILE,int & MAXHITS, int & MAXMISMATCH,unsigned & MAX_READS);
void ReplaceCtoT(char* Read);
inline void ReplaceGtoA(char* Read);
	
int main (int argc, char* argv[]) 
{

	time_t Start_Time,End_Time;
	time(&Start_Time);
	char *CT_GENOME,*GA_GENOME,*GENOME;
	char* INPUTFILE;
	char* OUTPUTFILE;
	unsigned MAX_READS=0;
	Parse_Command_line(argc, argv,GENOME,CT_GENOME,GA_GENOME,INPUTFILE,OUTPUTFILE,MAXHITS, MAX_MISMATCHES,MAX_READS);
	char Char_To_C[255];
	Char_To_C['A']='T';Char_To_C['a']='t';
	Char_To_C['C']='G';Char_To_C['c']='g';
	Char_To_C['G']='C';Char_To_C['g']='c';
	Char_To_C['T']='A';Char_To_C['t']='a';
	Char_To_C['N']='N';Char_To_C['n']='n';

	string tmpList1[100];
	string tmpList2[100];
	string tmpList3[100];
	string tmpList4[100];
	char Header[MAX_REC_SIZE];
//{-------------------------------- FM File Stuff ---------------------------------------------------
	FMFILES CT,GA;
	MEMLOOK MLookCT,MLookGA;
	INFILE F;
	MMPool *mmPool;
	LEN L;L.IGNOREHEAD=0;

	string CTS=GENOME,GAS=GENOME;
	if (GENOME)
	{
		CTS+="-CtoT";GAS+="-GtoA";
		CT_GENOME=(char*)CTS.c_str();GA_GENOME=(char*)GAS.c_str();
	}

	Build_Names(CT_GENOME,CT);Build_Names(GA_GENOME,GA);
	BWT *fwfmiCT,*revfmiCT,*fwfmiGA,*revfmiGA;
	Offset_Record Genome_OffsetsCT[80],Genome_OffsetsGA[80];
	unsigned Location_ArrayCT[80],Location_ArrayGA[80];

//------------ Analyse read file -------------
	F.Input_File=File_Open(INPUTFILE,"r");
	Output_File=File_Open(OUTPUTFILE,"w");
	F.Buffer=(char*) malloc(60000);
	char* Buffer_CT=(char*) malloc(FILEBUFSIZE);
	char* Buffer_RevCT=(char*) malloc(FILEBUFSIZE);
	char* Buffer_GA=(char*) malloc(FILEBUFSIZE);
	char* Buffer_RevGA=(char*) malloc(FILEBUFSIZE);

	Analyze_File(F,L);//Get File type, read length etc..
	Split_Read(L.STRINGLENGTH_ORG,L);//calculate read portions for Batman algo..
	printf("Loading C->T Index..");Load_FM_and_Sample(fwfmiCT,revfmiCT,mmPool,CT);printf("[Done]\n");//load FM index structures..
	printf("Loading A->G Index..");Load_FM_and_Sample(fwfmiGA,revfmiGA,mmPool,GA);printf("[Done]\n");
	
	int Genome_CountCT= Load_Location(CT.LOCATIONFILE,Genome_OffsetsCT,Location_ArrayCT);//Load informatin to decode locations..
	int Genome_CountGA= Load_Location(GA.LOCATIONFILE,Genome_OffsetsGA,Location_ArrayGA);
	Init(F.SOLID,0);//Set SOLiD /Solexa algo for batman..
	int LOOKUPSIZE=Get_Lookup_Size(MAX_MISMATCHES,L.STRINGLENGTH);//initialize lookup tables for batman..;
	MLookCT.Lookupsize=MLookGA.Lookupsize=LOOKUPSIZE;
	Build_Tables(fwfmiCT,revfmiCT,MLookCT);
	Build_Tables(fwfmiGA,revfmiGA,MLookGA);

//}-------------------------------- FM File Stuff ---------------------------------------------------

	unsigned Actual_Tag=0;
	unsigned Mapped=0;
	READ R,R_CT,R_GA;//+ strand reads..
	READ R_Rev_CT,R_Rev_GA;//Converted reverse strands..
	BATREAD B;
	GUESS G;
	MEMX MF,MC,MF_Rev,MF_GA,MC_GA;
	OUTPUT O_CT,O_GA,O_RevCT,O_RevGA;

	O_CT.SAM=0;
	O_CT.PLUSSTRAND=0;
	O_CT.MaxHits=MAXHITS;
	O_CT.Offset=0;
	O_CT.Location_Array=Location_ArrayCT;
	O_CT.Genome_Offsets=Genome_OffsetsCT;
	O_CT.Genome_Count=Genome_CountCT;
	O_CT.FILETYPE=FQ;
	O_CT.Length_Array[0]=O_CT.Length_Array[1]=O_CT.Length_Array[2]=L.STRINGLENGTH;
	O_RevCT=O_CT;
	O_GA=O_CT;O_RevGA=O_CT;
	O_CT.Buffer=Buffer_CT;//Buffers to write batman output ...
	O_GA.Buffer=Buffer_GA;
	O_RevCT.Buffer=Buffer_RevCT;
	O_RevGA.Buffer=Buffer_RevGA;

	char ONEFMINDEX=FALSE;

//----------- initialise memory structures ------------------
	B.IGNOREHEAD=L.IGNOREHEAD;B.StringLength=L.STRINGLENGTH;
	MF.Lookupsize=MF_Rev.Lookupsize=LOOKUPSIZE;
	MF_GA.Lookupsize=MC_GA.Lookupsize=LOOKUPSIZE;
	Copy_MEM(MLookCT,MF,MC,MAX_MISMATCHES);
	Copy_MEM(MLookCT,MF_Rev,MC,MAX_MISMATCHES);
	Copy_MEM(MLookGA,MF_GA,MC_GA,MAX_MISMATCHES);

//----------- Mapping Loop  ------------------
	int Progress=0;unsigned Number_of_Tags=INITIAL_PROGRESS_READS;
	Init_Progress();

	while (Read_Tag(F.Input_File,F.FILETYPE,R))
	{
		R.Tag_Number=1;R.Read_Number=Actual_Tag;
		if(++Actual_Tag == MAX_READS) break;
		Progress++;
		if (Progress>Number_of_Tags) 
		{
			off64_t Current_Pos=ftello64(F.Input_File);
			off64_t Average_Length=Current_Pos/Actual_Tag+1;//+1 avoids divide by zero..
			Number_of_Tags=(F.File_Size/Average_Length)/20;
			Progress=0;
			Show_Progress(Current_Pos*100/F.File_Size);
		}
		int Hits;
		bool Header_Printed = false;
		R_CT=R_Rev_CT=R_Rev_GA=R_GA=R;
		for (int i=0;i<=B.StringLength-1;i++){R_Rev_GA.Tag_Copy[B.StringLength-1-i]=R_Rev_CT.Tag_Copy[B.StringLength-1-i]=Char_To_C[R.Tag_Copy[i]];};//Reverse complement reads 
		int cntHit1 = 0, cntHit2 = 0, cntHit3 = 0, cntHit4 = 0;
		int firstmismatch1 = 0, firstmismatch2 = 0, firstmismatch3 = 0, firstmismatch4 = 0;
		int lastMismatch1 = -1, lastMismatch2 = -1, lastMismatch3 = -1, lastMismatch4 = -1;

//------------- C->T scan Start ----------------------------------------------
		Hits=0;
		ReplaceCtoT(R_CT.Tag_Copy);
		Process_Read(R_CT,B,MF,MC);
		if(Hits=Map_Strand(MAX_MISMATCHES,MAXHITS,L,fwfmiCT,revfmiCT,MF)){Mapped++;}
		Print_Hits(MF,L,Output_File,ONEFMINDEX,revfmiCT,fwfmiCT,O_CT);
		*O_CT.Buffer_End='&';*(O_CT.Buffer_End+1)='\n';

		char Banner[2000];
		int D;for(D=0;R.Description[D]!='\n' && R.Description[D]!='\r';D++);R.Description[D]=0;//Make read to C string
		R.Tag_Copy[L.STRINGLENGTH]=0;R.Quality[L.STRINGLENGTH]=0;
		sprintf(Banner,"%s\t%s\t%s\t%d:%d:%d:%d:%d:%d:%d\n",R.Description,R.Tag_Copy,R.Quality,MF.Stats[0],MF.Stats[1],MF.Stats[2],MF.Stats[3],MF.Stats[4],MF.Stats[5],MF.Stats[6]);
		Get_Records(O_CT.Buffer, tmpList1,firstmismatch1,lastMismatch1,cntHit1,'1');
		Print_List(Header_Printed,cntHit1,firstmismatch1,lastMismatch1,tmpList1,Banner);
//------------- G->A scan Start ----------------------------------------------
		Hits=0;
		ReplaceGtoA(R_GA.Tag_Copy);
		Process_Read(R_GA,B,MF_GA,MC);
		if(Hits=Map_Strand(MAX_MISMATCHES,MAXHITS,L,fwfmiGA,revfmiGA,MF_GA)){Mapped++;}
		Print_Hits(MF_GA,L,Output_File,ONEFMINDEX,revfmiGA,fwfmiGA,O_GA);
		*O_GA.Buffer_End='&';*(O_GA.Buffer_End+1)='\n';
		Get_Records(O_GA.Buffer, tmpList2,firstmismatch2,lastMismatch2,cntHit2,'2');
		Print_List(Header_Printed,cntHit2,firstmismatch2,lastMismatch2,tmpList2,Banner);
//------------- C->T reversed scan Start ----------------------------------------------
		Hits=0;
		ReplaceCtoT(R_Rev_CT.Tag_Copy);
		Process_Read(R_Rev_CT,B,MF,MC);
		if(Hits=Map_Strand(MAX_MISMATCHES,MAXHITS,L,fwfmiCT,revfmiCT,MF)){Mapped++;}
		Print_Hits(MF,L,Output_File,ONEFMINDEX,revfmiCT,fwfmiCT,O_RevCT);
		*O_RevCT.Buffer_End='&';*(O_RevCT.Buffer_End+1)='\n';
		Get_Records(O_RevCT.Buffer, tmpList3,firstmismatch3,lastMismatch3,cntHit3,'3');
		Print_List(Header_Printed,cntHit3,firstmismatch3,lastMismatch3,tmpList3,Banner);
//------------- G->A reversed scan Start ----------------------------------------------
		Hits=0;
		ReplaceGtoA(R_Rev_GA.Tag_Copy);
		Process_Read(R_Rev_GA,B,MF_GA,MC);
		if(Hits=Map_Strand(MAX_MISMATCHES,MAXHITS,L,fwfmiGA,revfmiGA,MF_GA)){Mapped++;}
		Print_Hits(MF_GA,L,Output_File,ONEFMINDEX,revfmiGA,fwfmiGA,O_RevGA);
		*O_RevGA.Buffer_End='&';*(O_RevGA.Buffer_End+1)='\n';
		Get_Records(O_RevGA.Buffer, tmpList4,firstmismatch4,lastMismatch4,cntHit4,'4');
		Print_List(Header_Printed,cntHit4,firstmismatch4,lastMismatch4,tmpList4,Banner);
//------------- G->A reversed scan Finish ----------------------------------------------



	}
	Done_Progress();
	time(&End_Time);printf("Total Reads - %u\n Time Taken  - %.0lf Seconds ..\n ",Actual_Tag,difftime(End_Time,Start_Time));

}
/* Given a base name, will find file names of related FM Index */
void Build_Names(const char* Genome_Name,FMFILES & FM)//build FM index file names...
{
	int Last_Dash=0;
	char* Name=(char*)Genome_Name;
	for(int i=0;Name[0]!=0;Name++,i++)
	{
		if (Name[0]=='/') 
		{
			Last_Dash=i;
		}
	}
	Name=(char*)Genome_Name+Last_Dash;Last_Dash++;

	char* Command_Line_Buffer=(char*)malloc(5000);
	FM.REVBWTINDEX = (char*)Command_Line_Buffer;
	//if(Last_Dash) Last_Dash=Genome_Name-N+1; else Genome_Name--;
	if(!Last_Dash) Name--;
	strncpy(FM.REVBWTINDEX,Genome_Name,Last_Dash);
	FM.REVBWTINDEX[Last_Dash+0]='r';FM.REVBWTINDEX[Last_Dash+1]='e';FM.REVBWTINDEX[Last_Dash+2]='v';
	strcpy(FM.REVBWTINDEX+Last_Dash+3, Name+1);
	strcat(FM.REVBWTINDEX+Last_Dash+3,".bwt"); 

	FM.BWTFILE=FM.REVBWTINDEX+500;
	strncpy(FM.BWTFILE,Genome_Name,Last_Dash);
	strcpy(FM.BWTFILE+Last_Dash,Name+1);
	strcat(FM.BWTFILE+Last_Dash,".bwt"); 


	FM.REVOCCFILE = FM.BWTFILE+500;
	strncpy(FM.REVOCCFILE,Genome_Name,Last_Dash);
	FM.REVOCCFILE[Last_Dash+0]='r';FM.REVOCCFILE[Last_Dash+1]='e';FM.REVOCCFILE[Last_Dash+2]='v';
	strcpy(FM.REVOCCFILE+Last_Dash+3,Name+1);
	strcat(FM.REVOCCFILE+Last_Dash+3,".fmv"); 


	FM.OCCFILE=FM.REVOCCFILE+500;			
	strncpy(FM.OCCFILE,Genome_Name,Last_Dash);
	strcpy(FM.OCCFILE+Last_Dash,Name+1);
	strcat(FM.OCCFILE+Last_Dash,".fmv"); 

	FM.SAFILE=FM.OCCFILE+500;			
	strncpy(FM.SAFILE,Genome_Name,Last_Dash);
	strcpy(FM.SAFILE+Last_Dash,Name+1);
	strcat(FM.SAFILE+Last_Dash,".sa");

	FM.REVSAFILE = FM.SAFILE+500;
	strncpy(FM.REVSAFILE,Genome_Name,Last_Dash);
	FM.REVSAFILE[Last_Dash+0]='r';FM.REVSAFILE[Last_Dash+1]='e';FM.REVSAFILE[Last_Dash+2]='v';
	strcpy(FM.REVSAFILE+Last_Dash+3,Name+1);
	strcat(FM.REVSAFILE+Last_Dash+3,".sa"); 

	FM.BINFILE=FM.REVSAFILE+500;			
	strncpy(FM.BINFILE,Genome_Name,Last_Dash);
	strcpy(FM.BINFILE+Last_Dash,Name+1);
	strcat(FM.BINFILE+Last_Dash,".bin");

	FM.LOCATIONFILE=FM.BINFILE+500;			
	strncpy(FM.LOCATIONFILE,Genome_Name,Last_Dash);
	strcpy(FM.LOCATIONFILE+Last_Dash,Name+1);
	strcat(FM.LOCATIONFILE+Last_Dash,".ann.location");

	FM.PACFILE=FM.LOCATIONFILE+500;			
	strncpy(FM.PACFILE,Genome_Name,Last_Dash);
	strcpy(FM.PACFILE+Last_Dash,Name+1);
	strcat(FM.PACFILE+Last_Dash,".pac");
}
	
int getMismatch(char* s2) 
{
	int tab=0;
	for (int i=0;s2[i]!='\n';i++)
	{
		if (s2[i]=='\t') { 
			tab++;
			//continue;
		}
		if (tab==4)
		{
			return atoi(s2+i+1);
		}
		else if (tab>4) break;
	}
	return 0;
}

bool Get_Records(char* File_Buffer, string*  tmpList,int & firstmismatch,int & lastMismatch,int & cntHit,char Count)
{
	int t=0,Skip=0;
	char Buffer[MAX_REC_SIZE];
	while(sscanf(File_Buffer,"%[^\n]%n",Buffer,&Skip)) 
	{
		if(Buffer[0] == '&') break;
		Buffer[0]=Count;//Buffer[Skip]='\n';
		tmpList[t++]=Buffer;
		if (cntHit==0) firstmismatch = getMismatch(Buffer);
		cntHit++;
		if (cntHit==MAXHITS)lastMismatch = getMismatch(Buffer);
		File_Buffer+=(Skip+1);
	}
	return true;
}

inline void Print_List(bool & Header_Printed,int cntHit,int firstmismatch,int lastMismatch,string* tmpList,char* Header)
{
	if ((cntHit<MAXHITS || firstmismatch!=lastMismatch)) 
	{
		if (!Header_Printed) {fprintf(Output_File,"@\n%s",Header);Header_Printed = true;}
		int currCnt=0;
		while(cntHit>0) 
		{
			fprintf(Output_File,"%s\n",tmpList[--cntHit].c_str());
			if(currCnt++==MAXHITS) break;
		}
	}
}

inline void ReplaceCtoT(char* Read)
{
        for (int i=0;Read[i]!='\n';i++)
        {
                if (Read[i] == 'C' || Read[i] == 'c') Read[i]='T';
        }
}

inline void ReplaceGtoA(char* Read)
{
        for (int i=0;Read[i]!='\n';i++)
        {
                if (Read[i] == 'G' || Read[i] == 'g') Read[i]='A';
        }
}

option Long_Options_Decode[]=
{
	{"help",0,NULL,'h'},
	{"inputfile",1,NULL,'i'},
	{"outputfile",1,NULL,'o'},
	{"CT",1,NULL,'C'},
	{"GA",1,NULL,'G'},
	{"genome",1,NULL,'g'},
	{"maxhits",1,NULL,'m'},
	{"maxtags",1,NULL,'t'},
	{0,0,0,0}
};


void Parse_Command_line(int argc, char* argv[],char* & GENOME,char* & CT_GENOME,char* & GA_GENOME,char* & INPUTFILE,char* & OUTPUTFILE,int & MAXHITS, int & MAXMISMATCH,unsigned & MAX_READS)
{
	int Current_Option=0;
	const char* Short_Options ="hi:C:G:g:o:n:m:t:";//allowed options....
	char* This_Program = argv[0];//Current program name....
	const char* Help_String=
"Parameters:\n"
" --help | -h\t\t\t\t Print help\n"
" --inputfile | -i <filename>\t\t Name of input file\n"
" --CT | -C | C->T genome...\n"
" --GA | -G | G->A genome...\n"
" --genome | -g | G->A genome...\n"
" --outputfile | -o <filename>\t\t Name of output file\n"
" --maxmis | -n <integer>\t\t Number of mismatches\n"
" --maxhits | -m <integer>\t\t Max hits per set...\n"
;
	char* Name;int Last_Dash;char* Genome_Name;
	char* Command_Line_Buffer=(char*)malloc(5000);
	GENOME=0;

	for(;;)	
	{
		Current_Option=getopt_long(argc, argv, Short_Options, Long_Options_Decode, NULL);
		if (Current_Option == -1 ) break;
		switch(Current_Option)
		{
			case 'h':
				printf("%s \n",Help_String);exit(0);
			case 'i':
				INPUTFILE=optarg;
				break;
			case 'o':
				OUTPUTFILE=optarg;
				break;
			case 'C':
				CT_GENOME=optarg;
				break;
			case 'g':
				GENOME=optarg;
				break;
			case 'G':
				GA_GENOME=optarg;
				break;
			case 'n':
				MAXMISMATCH=atol(optarg);
				break;
			case 'm':
				MAXHITS=atol(optarg);
				break;
			case 't':
				MAX_READS=atol(optarg);
				break;
			default:
				printf("%s \n",Help_String);
				exit(0);
		}
	}	
}

