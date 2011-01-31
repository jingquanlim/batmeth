//Notes : remove MAX_HITS_ALLOWED
//is chr labeling fixed?
//./split out2 ../hg19/hg19.fa 2 y out.split
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>


#define MAX_HITS_ALLOWED 251
#define CHROMSIZE 30
#define BATBUF 2000
#define MAXTAG 500 

struct Mismatch
{
	char Base;
	int Pos;
};
using namespace std;

//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------
FILE* File_Open(const char* File_Name,const char* Mode);
void Show_Progress(float Percentage);
#define Init_Progress()	printf("======================]\r[");//progress bar....
#define Done_Progress()	printf("\r[++++++++ 100%% ++++++++]\n");//progress bar....
inline unsigned char Hash(char* S);
int Get_String_Length(FILE* INFILE);
unsigned Unique_Count=0;
unsigned u=0;
unsigned Total_Reads;//total hits in the input file...
unsigned Tot_Unique_Org=0;//total unique hits obtained
unsigned Tot_Unique_Remdup=0;//total unique hits obtained after removing dups...
unsigned met_cnt=0;
unsigned non_met_cnt=0;

//}-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------
char Char2Comp[255];
const int INITIAL_PROGRESS_READS =1000;//progress bar initial count..

int main(int argc, char* argv[])
{
	time_t Start_Time,End_Time;
	printf("Split V 1.01\n");
	Char2Comp['A']=Char2Comp['a']='T';
	Char2Comp['C']=Char2Comp['c']='G';
	Char2Comp['G']=Char2Comp['g']='C';
	Char2Comp['T']=Char2Comp['t']='A';
	if (argc >1)
	{

		try
		{
			time(&Start_Time);
			char* Input_Name=argv[1];
			int UPPER_MAX_MISMATCH=atoi(argv[3]);
			bool Single_Output_File;if(*argv[4]=='y') Single_Output_File=true; else Single_Output_File=false;/*jq: set to default*/ 
			char* Output_Name=argv[5];
			string G=argv[2];G+=".bin";
			string L=argv[2];L+=".ann.location";

			FILE* BINFILE=File_Open(G.c_str(),"r");
			FILE* Location_File=File_Open(L.c_str(),"r");
			FILE* INFILE=File_Open(Input_Name,"r");
			FILE* OUTFILE=File_Open(Output_Name,"w");
			fseek(BINFILE, 0L, SEEK_END);off64_t Genome_Size=ftello64(BINFILE);rewind(BINFILE);//load original genome..
			char* Org_Genome=new char[Genome_Size];if(!Org_Genome) throw("Insufficient memory to load genome..\n"); 
			if(!fread(Org_Genome,Genome_Size,1,BINFILE)) throw ("Error reading file..\n");/*jq: format of genome*/
			char* Marked_Genome=new char[Genome_Size+1];if(!Marked_Genome) throw("Insufficient memory to Mark genome..\n"); 

			printf("Splitting Genome..\n");
			struct Offset_Record
			{
				char Genome[40];
				unsigned Offset;
			} Genome_Offsets[80];
			int Genome_Count=0;

			while (fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File)!=0 && Genome_Count<80)
			{
				Genome_Offsets[Genome_Count].Offset=atoi(Genome_Offsets[Genome_Count].Genome);
				fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File);
				for(int i=0;i<40;i++) 
				{
					if (Genome_Offsets[Genome_Count].Genome[i] == '\n' ||Genome_Offsets[Genome_Count].Genome[i] == '\r')
					{ 
						Genome_Offsets[Genome_Count].Genome[i]=0;
						break;
					} 
				}
				Genome_Count++;	
			}
			Genome_Count--;

			char* Split_Point=Org_Genome;//split and write...
			struct Gene_Hash
			{
				char* Genome;
				int Index;
			};
			Gene_Hash Genome_List[UCHAR_MAX];for (int i=0;i<UCHAR_MAX;Genome_List[i++].Index=CHAR_MAX)/*for collision detection*/;
			for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
			{
				printf("%s, ",Genome_Offsets[i].Genome);
				unsigned char H=Hash(Genome_Offsets[i].Genome);
				if(Genome_List[H].Index != CHAR_MAX) throw ("Collision in Hash: Please rename your chromosomes..\n");
				Genome_List[H].Genome=(Split_Point+=Genome_Offsets[i].Offset);
				Genome_List[H].Index=i;
			}
			printf("Loaded\n");


			/*unsigned Split_Point=0;//split and write...
			  for ( int i=0;i<Genome_Count;i++)
			  {
			  Split_Point+=Genome_Offsets[i].Offset;
			  string S=Genome_Offsets[i].Genome;
			  FILE *Chr_File=File_Open(S.c_str(),"w");
			  fwrite(&Org_Genome[Split_Point],1,Genome_Offsets[i+1].Offset,Chr_File);
			  fclose(Chr_File);
			  }*/


			int true_matches=0;unsigned pos;
			string readString="";
			int hitType=0;
			int n_mismatch[MAX_HITS_ALLOWED];
			string hits[MAX_HITS_ALLOWED];
			unsigned hitsToken[MAX_HITS_ALLOWED][5];
			char Comp_String[MAXTAG];for (int i=1;i<MAXTAG;Comp_String[i++]=0);
			//start to read batman hit file
			int jqread=0,Read_No;
			char Buf[BATBUF],s2t[BATBUF],Dummy[BATBUF],forReadString[BATBUF],rcReadString[BATBUF],Chrom[CHROMSIZE],Strand;

			int Progress=0;unsigned Number_of_Tags=INITIAL_PROGRESS_READS;
			fseek(INFILE, 0L, SEEK_END);off64_t File_Size=ftello64(INFILE);rewind(INFILE);//
			Init_Progress();

			int Read_Len=Get_String_Length(INFILE);
			fgets(Buf,BATBUF,INFILE);//read first header marker..
			while (!feof(INFILE)) 
			{
				Total_Reads++;
				Progress++;
				if (Progress>Number_of_Tags) 
				{
					off64_t Current_Pos=ftello64(INFILE);
					off64_t Average_Length=Current_Pos/Total_Reads+1;//+1 avoids divide by zero..
					Number_of_Tags=(File_Size/Average_Length)/20;
					Progress=0;
					Show_Progress(Current_Pos*100/File_Size);
				}
				fgets(s2t,BATBUF,INFILE);//read description..
				int cntHit=0;
				fgets(Buf,BATBUF,INFILE);
				while(!feof(INFILE) && Buf[0]!='@')//while in the record.. 
				{
					hits[cntHit++]=Buf; //read hit info..
					fgets(Buf,BATBUF,INFILE);
					assert(MAX_HITS_ALLOWED > cntHit);
				}

				sscanf(s2t,"%s%s%s",Dummy,forReadString,Dummy);
				for (int i=0;i<Read_Len;i++) {rcReadString[Read_Len-i-1]=Char2Comp[forReadString[i]];}

				char* Genome;
				for (int j=0; j<cntHit; j++) 
				{
					//we check for the mismatches information here and decide if we should carry them over for future computation

					assert(!hits[j].empty());
					sscanf(hits[j].c_str(),"%d%s%s%u%d%d%d",&hitType,Chrom,&Strand,&pos,&true_matches,&Read_Len,&Read_No);
					true_matches=0;
					char const *C=hits[j].c_str();char Mark,Mis_Letter;
					int Num,Skip,k=0,l=0;
					Mismatch Mismatches[MAX_HITS_ALLOWED];
					Mismatch N[MAX_HITS_ALLOWED];
					for (int Tab_Count=0;Tab_Count!=7;C++) if (*C=='\t') Tab_Count++; //seek mismatches..
					while(sscanf(C,"%d%c%c%n",&Num,&Mark,&Mis_Letter,&Skip)==3)	
					{
						if ('>' == Mark)
						{
							if (!Comp_String[Mis_Letter]) 
							{
								true_matches++;
								if (true_matches>UPPER_MAX_MISMATCH) {
									hits[j].clear();
									break;
								}
							}
							//Mismatches[k].Base=Mis_Letter;
							//Mismatches[k++].Pos=Num;
						}
						else
						{
							assert(Mark == ':');
							N[l].Base=Mis_Letter;
							Comp_String[Mis_Letter]=1;
							N[l++].Pos=Num;
							
						}
						C+=(Skip+1);
					}
					while(l)Comp_String[N[--l].Pos]=0;//restore Comp_String...

					if (hitType==3 || hitType==4) readString = rcReadString;
					else readString=forReadString;
					unsigned char H=Hash(Chrom);
					Genome=Genome_List[H].Genome;//load current genome..
					int Hash_Index=Genome_List[H].Index;//load current genome..
					hitsToken[j][0]=hitType;hitsToken[j][1]=H;hitsToken[j][2]=Strand;hitsToken[j][3]=pos;
					for (int k=0; k<Read_Len; k++) 
					{
						if (pos+k >= Genome_Offsets[Hash_Index+1].Offset) break;
						if ((hitType==1 || hitType==3) && (readString[k]=='C' && toupper(Genome[pos+k])=='T')) true_matches++;
						else if ((hitType==2 || hitType==4) && (readString[k]=='G' && toupper(Genome[pos+k])=='A')) true_matches++;
					}

					hitsToken[j][4] = true_matches;
				}
				for (int i=0; i<cntHit;i++) 
				{
					if (hits[i].empty()) continue;
					for (int j=i+1; j<cntHit;j++) 
					{
						if (hits[j].empty()) continue;
						bool same = true;
						for (int z=1;z<=4;z++) 
						{
							if (hitsToken[i][z]!=hitsToken[j][z]) 
							{
								same = false;
								break;
							}
						}
						if (same) hits[j].clear();		
					}
				}

				//making sure the hit is unique
				int min_Mis=UPPER_MAX_MISMATCH+1;
				int ind = -1;
				short gdHit = 0;
				for (int i=0; i<cntHit; i++) 
				{
					if (!hits[i].empty()) 
					{
						int tempNum = hitsToken[i][4];
						if (tempNum<min_Mis) 
						{
							min_Mis=tempNum;
							ind=i;
						}
					}
				}
				for (int i=0; i<cntHit; i++)
				{
					if (!hits[i].empty() && hitsToken[i][4]==min_Mis) gdHit++;
				}
				if (gdHit==1 && ind!=-1) 
				{

					Tot_Unique_Org++;
					string s3t = hits[ind];
					unsigned *HL = hitsToken[ind];char Flag=0;
					unsigned char H=HL[1];Genome=Genome_List[H].Genome;//load current genome..
					unsigned G_Skip=Genome-Org_Genome;
					if(HL[0] == '1' || HL[0] == '3' ) Flag=2; else Flag=4;
					char Mark=Marked_Genome[HL[3]+G_Skip];
					if (!Mark || !(Mark & Flag))//Remove duplicate hits..
					{

						Genome=Genome_List[H].Genome;//load current genome..
						int Hash_Index=Genome_List[H].Index;//load current genome..
						unsigned pos=HL[3];

						int hitType=HL[0];
						if (hitType==3 || hitType==4) readString = rcReadString;
						else readString=forReadString;
						for (int k=0; k<Read_Len; k++) 
						{
							if (pos+k >= Genome_Offsets[Hash_Index+1].Offset) break;
							if ((hitType==1 || hitType==3)&& toupper(Genome[pos+k+1])=='C') 
							{
								if (readString[k]=='C') met_cnt++;
								else if (readString[k]=='T' ) non_met_cnt++;
							}
							else if ((hitType==2 || hitType==4) && toupper(Genome[pos+k+1])=='G')
							{
								if (readString[k]=='G') met_cnt++;
								else if (readString[k]=='A' ) non_met_cnt++;
							}
						}

						fprintf(OUTFILE,"@\n%s%u\t%s\t%c\t%u\t%u\t%d\t%d\n",s2t,HL[0],Genome_Offsets[Genome_List[HL[1]].Index].Genome,HL[2],HL[3],min_Mis,Read_Len,Read_No);
						Tot_Unique_Remdup++;
					}
					Marked_Genome[HL[3]+G_Skip] |= Flag;
				}
			}
			Done_Progress();
		}
		catch(char* Err)
		{
			printf(Err);
			exit(-1);
		}
		printf("Met_loc / Non-Met_loc : %u/%u\n ", met_cnt,non_met_cnt);
		printf("Methylation Rate : %.01f \n",(float(non_met_cnt)*100)/float((non_met_cnt+met_cnt)));
		printf("Total Reads/Total mapped/Total uniqu mapped : %u/%u/%u \n",Total_Reads,Tot_Unique_Org,Tot_Unique_Remdup);

		time(&End_Time);printf("Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	}
	else {printf("Command Format : splitgenome Base_Name\n");exit(-1);}

}

inline unsigned char Hash(char* S)
{
	unsigned char C=0;
	for (int i=2;S[i];C+=i*S[i++]);
	return C;
}

int Get_String_Length(FILE* INFILE)
{
	char Buf[BATBUF],Dummy[BATBUF],Tag[BATBUF];int L;
	fgets(Buf,BATBUF,INFILE);
	fgets(Buf,BATBUF,INFILE);
	sscanf(Buf,"%s%s",Dummy,Tag);	
	for (L=0;Tag[L];L++);
	rewind (INFILE);
	return L;
}
//{----------------------------------- FILE HANDLING ---------------------------------------------------------

/* 
 * ===  FUNCTION  ======================================================================
 *         Name:  File_Open
 *  Description:  Open a file:
 *  Mode - "w" create/Write text from top    - "wb" Write/Binary  -"w+" open text read/write            -"w+b" same as prev/Binary
 *         "r" Read text from top            - "rb" Read/Binary   -"r+" open text read/write/nocreate   -"r+b" same as prev/binary
 *       - "a" text append/write                                  -"a+" text open file read/append      -"a+b" open binary read/append
 *
 * =====================================================================================
 */
FILE* File_Open(const char* File_Name,const char* Mode)
{
	FILE* Handle;
	Handle=fopen64(File_Name,Mode);
	if (Handle==NULL)
	{
		printf("File %s Cannot be opened ....",File_Name);
		exit(1);
	}
	else return Handle;
}

//}----------------------------------- FILE HANDLING ---------------------------------------------------------


void Show_Progress(float Percentage)
{
	if (Percentage >98) return;
	printf("+%.0f%\b\b\b",Percentage);
	fflush(stdout);
}
