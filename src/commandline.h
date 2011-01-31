#ifndef CMDLN
#define CMDLN
void Parse_Command_line(int argc, char* argv[],BATPARAMETERS & P,LEN & L,FMFILES & FM);
void Parse_Command_line_Decode(int argc, char* argv[],char* & INPUTFILE,FMFILES & FM,char & SAM,char & PLUSSTRAND,int & Force_Maxhits,unsigned & Offset, char* & OUTPUTFILE,char &  VERIFY);
#endif
