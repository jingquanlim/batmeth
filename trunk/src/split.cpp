//Notes : remove MAX_HITS_ALLOWED
//is chr labeling fixed?
//./split out2 ../hg19/hg19.fa 2 y out.split
#include <string.h>
#include <cstdio>
#include <assert.h>
#include <cstdlib>
#include <string>
#include <math.h>
#include "limits.h"
#include <map>
#define MAX_HITS_ALLOWED 400
#define CHROMSIZE 30
#define BATBUF 2000
#define MAXTAG 500
#define ENTROPY_CUTOFF 0.25

struct Mismatch
{
	char Base;
	int Pos;
};
struct Gene_Hash
{
	char* Genome;
	int Index;
};
using namespace std;
bool Collision=false;
map <string,int> String_Hash;


//{-----------------------------  FUNCTION PRTOTYPES  -------------------------------------------------
FILE* File_Open(const char* File_Name,const char* Mode);
void Show_Progress(float Percentage);
#define Init_Progress()	printf("======================]\r[");//progress bar....
#define Done_Progress()	printf("\r[++++++++ 100%% ++++++++]\n");//progress bar....
inline unsigned char Hash(char* S);
int Get_String_Length(FILE* INFILE);
unsigned Unique_Count=0;
unsigned u=0;
unsigned Total_Reads=0;//total hits in the input file...
unsigned Progress_Reads=0;
unsigned Tot_Unique_Org=0;//total unique hits obtained
unsigned Tot_Unique_Remdup=0;//total unique hits obtained after removing dups...
unsigned met_cnt=0;
unsigned non_met_cnt=0;
bool REMOVE_DUP=true; //true to removeDup, false will not remove PCR-dup
unsigned Mismatch_Qual[255][255][255]; //[readLength][255][255]

//}-----------------------------   GLOBAL VARIABLES  -------------------------------------------------
char Char2Comp[255];
const int INITIAL_PROGRESS_READS =1000;//progress bar initial count..
float calc_Entropy(string readString, int L);
void Print_Mismatch_Quality(FILE* OUTFILE_MM, int L);

int main(int argc, char* argv[])
{
	time_t Start_Time,End_Time;
	printf("\nBatMeth: Split  v1.04\n");
	printf("Command Format :  split <FINAL_RESULT> GENOME <number of mismatches> <y for base-read, n for color-read> <all..temp(s)...from...batmeth>\n\n");
	Char2Comp['A']=Char2Comp['a']='T';
	Char2Comp['C']=Char2Comp['c']='G';
	Char2Comp['G']=Char2Comp['g']='C';
	Char2Comp['T']=Char2Comp['t']='A';
	Char2Comp['N']=Char2Comp['n']='N';
	int Genome_CountX=0;
	char* Output_Name=argv[1];
	
	if (argc >1 && !strcmp(argv[4], "y")) //this is branching into the solexa part
	{
		try
		{
			time(&Start_Time);
			char* Input_Name=argv[5];
			int UPPER_MAX_MISMATCH=atoi(argv[3]);
			//bool Single_Output_File;if(*argv[4]=='y') Single_Output_File=true; else Single_Output_File=false;/*jq: set to default*/ 
			string G=argv[2];G+=".bin";
			string L=argv[2];L+=".ann.location";

			FILE* BINFILE=File_Open(G.c_str(),"r");
			FILE* Location_File=File_Open(L.c_str(),"r");
			FILE* INFILE=File_Open(Input_Name,"r");
			FILE* OUTFILE=File_Open(Output_Name,"w");
			//FILE* OUTFILE_MM_QUALITY=File_Open(strcat(Output_Name,".mismatch"), "w");
			fseek(BINFILE, 0L, SEEK_END);off64_t Genome_Size=ftello64(BINFILE);rewind(BINFILE);//load original genome..
			char* Org_Genome=new char[Genome_Size];if(!Org_Genome) throw("Insufficient memory to load genome..\n"); 
			if(!fread(Org_Genome,Genome_Size,1,BINFILE)) throw ("Error reading file..\n");/*jq: format of genome*/
			char* Marked_Genome=new char[Genome_Size+1];if(!Marked_Genome) throw("Insufficient memory to Mark genome..\n"); 

			printf("Splitting Genome..\n");

			struct Offset_Record
			{
				char Genome[40];
				unsigned Offset;
			} Temp_OR; 

			while (fgets(Temp_OR.Genome,39,Location_File)!=0)//count genomes..
			{
				fgets(Temp_OR.Genome,39,Location_File);
				Genome_CountX++;	
			}
			rewind(Location_File);
			Offset_Record Genome_Offsets[Genome_CountX+1];

			int Genome_Count=0;
			while (fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File)!=0 && Genome_Count<=Genome_CountX)
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
			//Gene_Hash Genome_List[UCHAR_MAX];for (int i=0;i<UCHAR_MAX;Genome_List[i++].Index=CHAR_MAX)/*for collision detection*/;
			/*for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
			{
				printf("%s, ",Genome_Offsets[i].Genome);
				unsigned char H=Hash(Genome_Offsets[i].Genome);
				if(Genome_List[H].Index != CHAR_MAX)
				{
					Collision=true;
					break;
				       	//throw ("Collision in Hash: Please rename your chromosomes..\n");
				}
				Genome_List[H].Genome=(Split_Point+=Genome_Offsets[i].Offset);
				Genome_List[H].Index=i;
			}
			if(Collision)
			{
			}
			printf("Loaded\n");*/

			Gene_Hash Genome_List[Genome_Count];
			for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
			{
				printf("%s, ",Genome_Offsets[i].Genome);
				String_Hash[Genome_Offsets[i].Genome]=i;
				//unsigned char H=Hash(Genome_Offsets[i].Genome);
				Genome_List[i].Genome=(Split_Point+=Genome_Offsets[i].Offset);
				Genome_List[i].Index=i;
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
			
			//Requested by Jinyan
			long non_met_CG=0;
			long met_CG=0;
			long non_met_CHG=0;
			long met_CHG=0;
			long non_met_CHH=0;
			long met_CHH=0;

			//Added to removePCRdup from all temp files at once
			for (int tempFileNum=5; tempFileNum<argc; tempFileNum++) {
				printf("\nProcessing %d out of %d. File: %s\n\n", tempFileNum-4, argc-5, argv[tempFileNum]);
				INFILE=File_Open(argv[tempFileNum],"r");
			int Progress=0;unsigned Number_of_Tags=INITIAL_PROGRESS_READS;
			fseek(INFILE, 0L, SEEK_END);off64_t File_Size=ftello64(INFILE);rewind(INFILE);//
			Init_Progress(); Progress_Reads=0;

			int Read_Len=Get_String_Length(INFILE);
			fgets(Buf,BATBUF,INFILE);//read first header marker..
			while (!feof(INFILE)) 
			{
				Total_Reads++;
				Progress_Reads++;
				Progress++;
				if (Progress>Number_of_Tags) 
				{
					off64_t Current_Pos=ftello64(INFILE);
					off64_t Average_Length=Current_Pos/Progress_Reads+1;//+1 avoids divide by zero..
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
				//space-trouble in read description
				for(int i=0;s2t[i];i++) if(s2t[i]==' ') s2t[i]='_';
				sscanf(s2t,"%s%s%s",Dummy,forReadString,Dummy);
				for (int i=0;i<Read_Len;i++) {rcReadString[Read_Len-i-1]=Char2Comp[forReadString[i]];}

				char* Genome;

				char read_Methyl_Info[cntHit][Read_Len+1];

				for (int j=0; j<cntHit; j++) 
				{
					//We init the string which will contain the methylation status of a read wrt a putative genomic location
					read_Methyl_Info[j][Read_Len]='\0';
					//read_Methyl_Info[j][0]='\0';

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
							if (!Comp_String[Num]) 
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
							Comp_String[Num]=1;
							N[l++].Pos=Num;
							
						}
						C+=(Skip+1);
					}
					while(l)Comp_String[N[--l].Pos]=0;//restore Comp_String...

					if (hitType==3 || hitType==4) readString = rcReadString;
					else readString=forReadString;
					//unsigned char H=Hash(Chrom);
					int H=String_Hash[Chrom];
					Genome=Genome_List[H].Genome;//load current genome..
					int Hash_Index=Genome_List[H].Index;//load current genome..
					hitsToken[j][0]=hitType;hitsToken[j][1]=H;hitsToken[j][2]=Strand;hitsToken[j][3]=pos;
					for (int k=0; k<Read_Len; k++) 
					{
						read_Methyl_Info[j][k] = '=';
						if (pos+k >= Genome_Offsets[Hash_Index+1].Offset) {hits[j].clear();break;}

						char genome_Char = toupper(Genome[pos+k]);
						if (hitType==1 || hitType==3) {
						
							if (readString[k]=='C' && genome_Char=='C') read_Methyl_Info[j][k] = 'M';
							else if (readString[k]=='T' && genome_Char=='C') read_Methyl_Info[j][k] = 'U';
							else if (readString[k] != genome_Char) {
								read_Methyl_Info[j][k] = genome_Char;
								if (readString[k]=='C' && genome_Char=='T') true_matches++;
							}
						}
						else if (hitType==2 || hitType==4) {
							if (readString[k]=='G' && genome_Char=='G') read_Methyl_Info[j][k] = 'M';
							else if (readString[k]=='A' && genome_Char=='G') read_Methyl_Info[j][k] = 'U';
							else if (readString[k] != genome_Char) {
								read_Methyl_Info[j][k] = genome_Char;
								if (readString[k]=='G' && genome_Char=='A') true_matches++;
							}
						}
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
						if (tempNum==min_Mis) 
						{
							gdHit++;
						}
						else if (tempNum<min_Mis) 
						{
							min_Mis=tempNum;
							ind=i;
							gdHit=1;
						}
					}
				}
				if (gdHit==1 && ind!=-1 && calc_Entropy(readString, Read_Len)>= ENTROPY_CUTOFF) 
				{

					Tot_Unique_Org++;
					string s3t = hits[ind];
					unsigned *HL = hitsToken[ind];char Flag=0;
					unsigned char H=HL[1];Genome=Genome_List[H].Genome;//load current genome..
					unsigned G_Skip=Genome-Org_Genome;
					if(HL[0] == '1' || HL[0] == '3' ) Flag=2; else Flag=4;
					char Mark=Marked_Genome[HL[3]+G_Skip];
					if ((!Mark || !(Mark & Flag)) || !REMOVE_DUP)//Remove duplicate hits..
					{

						Genome=Genome_List[H].Genome;//load current genome..
						int Hash_Index=Genome_List[H].Index;//load current genome..
						unsigned pos=HL[3];
						
						int hitType=HL[0];
						if (hitType==3 || hitType==4) readString = rcReadString;
						else readString=forReadString;
						
						for (int k=0; k<Read_Len; k++) 
						{
							if (pos+k >= Genome_Offsets[Hash_Index+1].Offset || pos+k < 2) continue;

							char genome_Char = toupper(Genome[pos+k]);
							char string_char = toupper(readString[k]);
							char genome_CharFor1 = toupper(Genome[pos+k+1]);
							char genome_CharFor2 = toupper(Genome[pos+k+2]);
							char genome_CharBac1 = toupper(Genome[pos+k-1]);
							char genome_CharBac2 = toupper(Genome[pos+k-2]);
							
							if(hitType==1) {
								if(genome_Char=='C'&&string_char=='C') {
									if(genome_CharFor1=='G') met_CG++;
									else if(genome_CharFor1!='G' && genome_CharFor2!='G') met_CHH++;
									else if(genome_CharFor1!='G' && genome_CharFor2=='G') met_CHG++;
								}
								else if(genome_Char=='C'&&string_char=='T') {
									if(genome_CharFor1=='G') non_met_CG++;
									else if(genome_CharFor1!='G' && genome_CharFor2!='G') non_met_CHH++;
									else if(genome_CharFor1!='G' && genome_CharFor2=='G') non_met_CHG++;
								}
							}
							
							else if(hitType==2) {
								if(genome_Char=='G'&&string_char=='G') {
									if(genome_CharBac1=='C') met_CG++;
									else if(genome_CharBac1!='C' && genome_CharBac2!='C') met_CHH++;
									else if(genome_CharBac1!='C' && genome_CharBac2=='C') met_CHG++;
								}
								else if(genome_Char=='G'&&string_char=='A') {
									if(genome_CharBac1=='G') non_met_CG++;
									else if(genome_CharBac1!='C' && genome_CharBac2!='C') non_met_CHH++;
									else if(genome_CharBac1!='C' && genome_CharBac2=='C') non_met_CHG++;
								}
							}
							
							else if(hitType==3) {
								if(genome_Char=='C'&&string_char=='C') {
									if(genome_CharBac1=='G') met_CG++;
									else if(genome_CharBac1!='G' && genome_CharBac2!='G') met_CHH++;
									else if(genome_CharBac1!='G' && genome_CharBac2=='G') met_CHG++;
								}
								else if(genome_Char=='C'&&string_char=='T') {
									if(genome_CharBac1=='G') non_met_CG++;
									else if(genome_CharBac1!='G' && genome_CharBac2!='G') non_met_CHH++;
									else if(genome_CharBac1!='G' && genome_CharBac2=='G') non_met_CHG++;
								}
							}
							
							else if(hitType==4) {
								if(genome_Char=='G'&&string_char=='G') {
									if(genome_CharFor1=='C') met_CG++;
									else if(genome_CharFor1!='C' && genome_CharFor2!='C') met_CHH++;
									else if(genome_CharFor1!='C' && genome_CharFor2=='C') met_CHG++;
								}
								else if(genome_Char=='G'&&string_char=='A') {
									if(genome_CharFor1=='C') non_met_CG++;
									else if(genome_CharFor1!='C' && genome_CharFor2!='C') non_met_CHH++;
									else if(genome_CharFor1!='C' && genome_CharFor2=='C') non_met_CHG++;
								}
							}
						}

						fprintf(OUTFILE,"@\n%s%u\t%s\t%c\t%u\t%u\t%d\t%s\n",s2t,HL[0],Genome_Offsets[Genome_List[HL[1]].Index].Genome,HL[2],HL[3],min_Mis,Read_Len,read_Methyl_Info[ind]);
						Tot_Unique_Remdup++;
					}
					Marked_Genome[HL[3]+G_Skip] |= Flag;
				}
			}
			Done_Progress();
			//Print_Mismatch_Quality(OUTFILE_MM_QUALITY, Read_Len);
			}
			fclose(OUTFILE);
			
			printf("Raw count of Met_C in CG: %d\n",met_CG);
			printf("Raw count of Non_Met_C in CG: %d\n",non_met_CG);
			printf("Raw count of Met_C in CHG: %d\n",met_CHG);
			printf("Raw count of Non_Met_C in CHG: %d\n",non_met_CHG);
			printf("Raw count of Met_C in CHH: %d\n",met_CHH);
			printf("Raw count of Non_Met_C in CHH: %d\n",non_met_CHH);
		
		}
		catch(char* Err)
		{
			printf(Err);
			exit(-1);
		}
		//printf("@>\nMet_loc / Non-Met_loc : %u\t%u\n ", met_cnt,non_met_cnt);
		//printf("Methylation Rate : %.01f \n",(float(non_met_cnt)*100)/float((non_met_cnt+met_cnt)));
		printf("Total Reads/Total unique mapped/Total after PCR-dup-removal/OutputFile : %u\t%u\t%u\t%s \n",Total_Reads,Tot_Unique_Org,Tot_Unique_Remdup,Output_Name);

		time(&End_Time);printf("Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	}
	////////////////////////////////////////////////////SOLiD Starts here//////////////////////////////////////////////////////////
	else if (argc >1 && !strcmp(argv[4], "n")) //this is branching into the SOLiD part
	{
		try
		{
			time(&Start_Time);
			char* Input_Name=argv[5];
			int UPPER_MAX_MISMATCH=atoi(argv[3]);
			//bool Single_Output_File;if(*argv[4]=='y') Single_Output_File=true; else Single_Output_File=false;/*jq: set to default*/ 
			//string G=argv[2];G+=".bin";
			//string L=argv[2];L+=".ann.location";

			//FILE* BINFILE=File_Open(G.c_str(),"r");
			//FILE* Location_File=File_Open(L.c_str(),"r");
			FILE* INFILE=File_Open(Input_Name,"r");
			FILE* OUTFILE=File_Open(Output_Name,"w");
			//fseek(BINFILE, 0L, SEEK_END);off64_t Genome_Size=ftello64(BINFILE);rewind(BINFILE);//load original genome..
			//char* Org_Genome=new char[Genome_Size];if(!Org_Genome) throw("Insufficient memory to load genome..\n"); 
			//if(!fread(Org_Genome,Genome_Size,1,BINFILE)) throw ("Error reading file..\n");/*jq: format of genome*/
			//char* Marked_Genome=new char[Genome_Size+1];if(!Marked_Genome) throw("Insufficient memory to Mark genome..\n"); 

			printf("Splitting Genome..\n");
			struct Offset_Record
			{
				char Genome[40];
				unsigned Offset;
			} Genome_Offsets[Genome_CountX+1];
			int Genome_Count=0;

			/*while (fgets(Genome_Offsets[Genome_Count].Genome,39,Location_File)!=0 && Genome_Count<80)
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
			Genome_Count--;*/

			//char* Split_Point=Org_Genome;//split and write...
			/*
			struct Gene_Hash
			{
				char* Genome;
				int Index;
			};
			Gene_Hash Genome_List[UCHAR_MAX];for (int i=0;i<UCHAR_MAX;Genome_List[i++].Index=CHAR_MAX); //for collision detection
			for ( int i=0;i<Genome_Count;i++)//Stores the location in value corresponding to has..
			{
				printf("%s, ",Genome_Offsets[i].Genome);
				unsigned char H=Hash(Genome_Offsets[i].Genome);
				if(Genome_List[H].Index != CHAR_MAX) throw ("Collision in Hash: Please rename your chromosomes..\n");
				Genome_List[H].Genome=(Split_Point+=Genome_Offsets[i].Offset);
				Genome_List[H].Index=i;
			}
			printf("Loaded\n");
			*/

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
			char readString[100], genomeString[100];
			int hitType=0;
			int n_mismatch[MAX_HITS_ALLOWED];
			string hits[MAX_HITS_ALLOWED];
			unsigned hitsToken[MAX_HITS_ALLOWED][5];
			//start to read batman hit file
			int jqread=0,Read_No;
			char Buf[BATBUF],s2t[BATBUF],Dummy[BATBUF],forReadString[BATBUF],rcReadString[BATBUF],Chrom[CHROMSIZE],Strand;

			//Added to removePCRdup from all temp files at once
			for (int tempFileNum=5; tempFileNum<argc; tempFileNum++) {
				printf("\nProcessing %d out of %d. File: %s\n\n", tempFileNum-4, argc-5, argv[tempFileNum]);
				INFILE=File_Open(argv[tempFileNum],"r");
			int Progress=0;unsigned Number_of_Tags=INITIAL_PROGRESS_READS;
			fseek(INFILE, 0L, SEEK_END);off64_t File_Size=ftello64(INFILE);rewind(INFILE);//
			Init_Progress(); Progress_Reads=0;

			int Read_Len=Get_String_Length(INFILE);
			fgets(Buf,BATBUF,INFILE);//read first header marker..
			while (!feof(INFILE)) 
			{
				Total_Reads++;
				Progress_Reads++;
				Progress++;
				if (Progress>Number_of_Tags) 
				{
					off64_t Current_Pos=ftello64(INFILE);
					off64_t Average_Length=Current_Pos/Progress_Reads+1;//+1 avoids divide by zero..
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
				//space-trouble in read description
				//for(int i=0;s2t[i];i++) if(s2t[i]==' ') s2t[i]='_';
				//sscanf(s2t,"%s%s%s%s",Dummy,forReadString,Dummy, Dummy);
				//for (int i=1;i<Read_Len;i++) {rcReadString[i]=forReadString[Read_Len-i];}

				char* Genome;
				
				//start of stage-filtering
				bool stage_hit_mark[MAX_HITS_ALLOWED]; for (int t=0; t<MAX_HITS_ALLOWED; t++) stage_hit_mark[t]=false;
				for (int mm_stage=0; mm_stage<9; mm_stage++) {
					for (int j=0; j<cntHit; j++) 
					{
						//we check for the mismatches information here and decide if we should carry them over for future computation
			
						if(hits[j].empty()) continue;
						if (stage_hit_mark[j]) continue;

						sscanf(hits[j].c_str(),"%d%s%s%u%d%d%d%s%s",&hitType,Chrom,&Strand,&pos,&true_matches,&Read_Len,&Read_No, readString, genomeString);
						//skip a non-current stage hit
						if (true_matches!=mm_stage) continue;
						stage_hit_mark[j]=true; //mark this so we do not keep sscanf-ing it
						
						true_matches=0;
						
						char const *C=hits[j].c_str();char Mark,Mis_Letter;
						int Num,Skip,k=0,l=0;
						for (int Tab_Count=0;Tab_Count!=9;C++) if (*C=='\t') Tab_Count++; //seek tabs..

						bool mismatch_Position[Read_Len]; for (int i=0; i<Read_Len; i++) mismatch_Position[i]=false;
						
						while(sscanf(C,"%d%c%c%c",&Num,&Mark,&Mis_Letter,&Skip)==3)	
						{
							if ('>' == Mark && Num!='F')
							{
								printf("%d=",Num);
								mismatch_Position[Num]=true;
							}
							if (Skip!='\t') break;
							else C+=(Skip+1);
						}

						//unsigned char H=Hash(Chrom);
						int H=String_Hash[Chrom];
					//	Genome=Genome_List[H].Genome;//load current genome..
					//	int Hash_Index=Genome_List[H].Index;//load current genome..
						hitsToken[j][0]=hitType;hitsToken[j][1]=H;hitsToken[j][2]=Strand;hitsToken[j][3]=pos;
						for (int k=1; k<Read_Len; ) 
						{
							if(mismatch_Position[k]&&mismatch_Position[k-1]) 
							{
								if ((hitType==1 || hitType==3) && !(toupper(readString[k])=='T' && toupper(genomeString[k])=='C') && !(toupper(readString[k])==toupper(genomeString[k]))) true_matches++;
								else if ((hitType==2 || hitType==4) && !(toupper(readString[k])=='A' && toupper(genomeString[k])=='G') && !(toupper(readString[k])==toupper(genomeString[k]))) true_matches++;
								if(mismatch_Position[k+1]&&mismatch_Position[k+2]) k++;
							}
							k++;
						}
						hitsToken[j][4] = true_matches;
						if(true_matches>UPPER_MAX_MISMATCH) hits[j].clear();
					}
						
					for (int i=0; i<cntHit;i++) //dup hit fix
					{
						if (hits[i].empty() || !stage_hit_mark[i]) continue;
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
						if (!hits[i].empty() && stage_hit_mark[i]) 
						{
							int tempNum = hitsToken[i][4];
							if (tempNum==min_Mis) 
							{
								gdHit++;
							}
							else if (tempNum<min_Mis) 
							{
								min_Mis=tempNum;
								ind=i;
								gdHit=1;
							}
						}
					}
					if (gdHit==1 && ind!=-1) 
					{
						Tot_Unique_Org++;
						string s3t = hits[ind];
						//unsigned *HL = hitsToken[ind];char Flag=0;
						//unsigned char H=HL[1];Genome=Genome_List[H].Genome;//load current genome..
						//unsigned G_Skip=Genome-Org_Genome;
						//if(HL[0] == '1' || HL[0] == '3' ) Flag=2; else Flag=4;
						//char Mark=Marked_Genome[HL[3]+G_Skip];
						//if ((!Mark || !(Mark & Flag)) || true)//Remove duplicate hits..
						{

							//Genome=Genome_List[H].Genome;//load current genome..
							//int Hash_Index=Genome_List[H].Index;//load current genome..
							//unsigned pos=HL[3];

							//int hitType=HL[0];
							fprintf(OUTFILE,"@\n%s%s", s2t,s3t.c_str());//jqtest: output the ORIGINAL MISMATCH NUMBER
							//fprintf(OUTFILE,"@\n%s%u\t%s\t%c\t%u\t%u\t%d\t%d\n",s2t,HL[0],Genome_Offsets[Genome_List[HL[1]].Index].Genome,HL[2],HL[3],min_Mis,Read_Len,Read_No);
							Tot_Unique_Remdup++;
						}
						//Marked_Genome[HL[3]+G_Skip] |= Flag;
						break;
					}
				}
			}
			Done_Progress();
			fclose(INFILE);
			}
			
			fclose(OUTFILE);
		}
		catch(char* Err)
		{
			printf(Err);
			exit(-1);
		}
		//printf("@>\nMet_loc / Non-Met_loc : %u\t%u\n ", met_cnt,non_met_cnt);
		//printf("Methylation Rate : %.01f \n",(float(non_met_cnt)*100)/float((non_met_cnt+met_cnt)));
		printf("Total Reads/Total unique mapped/Total after PCR-dup-removal/OutputFile : %u\t%u\t%u\t%s \n",Total_Reads,Tot_Unique_Org,Tot_Unique_Remdup,Output_Name);

		time(&End_Time);printf("Time Taken  - %.0lf Seconds ..\n ",difftime(End_Time,Start_Time));
	}
	else {/*printf("Command Format :  split <Output_File> GENOME <number of mismatches> y <Final_result>\n");*/exit(-1);}

}

void Print_Mismatch_Quality(FILE* OUTFILE_MM, int L) {
	char bases[] = "ACGT";

	for (int i=0; i<L; i++) {
		for (int j=0; j<4; j++) {
			for (int k=0; k<4; k++) {
				fprintf(OUTFILE_MM,"%u\t", Mismatch_Qual[i][bases[j]][bases[k]]);
			}
		} fprintf(OUTFILE_MM,"\n");
	}

}

inline unsigned char Hash(char* S)
{
	assert(false);
	unsigned char C=0;
	for (int i=2;S[i];C+=i*S[i++]);
	return C;
}

inline float calc_Entropy (string readString, int L)
{ 
	short entropy_arr[255]={0};
	//int length=strlen(readString);
	for(int i=0; i<L; i++) entropy_arr[readString[i]]++;
	
	float entropy=0.0;
	for(int i=0; i<4; i++) {
		double p = 1.0*entropy_arr["ACGT"[i]]/L;
		if(p>0) entropy-=p*log(p);
	}
	return entropy;
}

int Get_String_Length(FILE* INFILE)
{
	char Buf[BATBUF],Dummy[BATBUF],Tag[BATBUF];int L;
	fgets(Buf,BATBUF,INFILE);
	fgets(Buf,BATBUF,INFILE);
	//space-trouble in read description
	for(int i=0;Buf[i];i++) if(Buf[i]==' ') Buf[i]='_';
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
