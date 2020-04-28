//To process CORSIKA showers generated using "Thining" option
//To compile: g++ test_With_Thining.cxx
//To run: ./a.out

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define	pi		    3.14159265	//value of pi.

#define RH		211285 
#define EH		217433
#define EE		3397
#define RE		3301
#define Mu_AddInfo_Low		75000	//To reject muon additional info stored in the DATxxxxx file
#define Mu_AddInfo_High		77000	//	,,	,,	,,

#define	Bytes_Per_Word		4	//No. of bytes for 1 word
#define Words_Per_SBlock	312 	//No. of words in a sub-block
#define	Words_Per_Particle	8	//No of words per particle
#define	SBlock_Per_Block	21	//No. of sub-blocks in one block

#define infilepath "DAT000001"
#define outfilepath "test.dat"

class PROCESS_SHOWER_THIN
{
    private:
        struct SPACE	//To read the blank space provided between data blocks
		{
			float word ;
		} ;
		struct SPACE space ;
        
        struct DATA		//To read data sub-block
		{
			float word[Words_Per_SBlock] ;
		} ;
        struct DATA data;
		struct DATA run_header;		//Run Header
		struct DATA event_header;	//Event Header
		struct DATA event_end;		//Event End
		struct DATA run_end;		//Run End
		
        struct PRI_PARTICLE //For primary particle
        {
            int evt_no ;    //Event number
            int id ;        //Particle ID
            float energy ;  //Total energy in GeV
            float zenith ;  //Corsika Zenith angle in radian
            float azimuth ; //Corsika Azimuth angle in radian :measured anticlockwise from north (X-axis) towards west (Y-axis)
        } ;
        struct PRI_PARTICLE Pparticle ;
        
		struct PARTICLE	//For air shower particles in the CORSIKA coordinate system: (X-axis:North),(Y-axis:West)
		{
			int CORSIKA_id ;	//Particle description code: (particle idx1000+....)
			int id ;	//Particle ID
			float px ;	//Momentum in x-direction in GeV/c
			float py ;	//Momentum in y-direction in GeV/c
			float pz ;	//Momentum in z-direction in GeV/c
			float x ;	//x position coordinate in cm
			float y ;	//y position coordinate in cm
			int t ;		//time since first interaction in nsec
            float wt ;	//Weight parameter for each particle
		} ;
		struct PARTICLE particle ;

    public:
		void Read_Data(const char*,const char*) ;	//Reading data file
} ;

void PROCESS_SHOWER_THIN::Read_Data(const char*infile, const char*outfile)
{
    int subblock ;	//No. of subblocks read
	subblock=0 ;
    
	FILE *in=fopen(infile,"rb") ;		//CORSIKA generated DATxxxxxx raw file
	FILE *out=fopen(outfile,"w") ;	//Output ascii file
	if(fread(&space,sizeof(struct SPACE),1,in)<0)	//Reading space at the start
    printf("Error reading %s: First space\n",infile) ;
	while(!feof(in))
	{
	    if(fread(&data,sizeof(struct DATA),1,in)<0) //Reading 1 data sub-block
        printf("Error reading %s: Sub-block\n",infile) ;
        switch((long int)data.word[0])
		{
			case RH:
				memcpy(&run_header,&data,sizeof(struct DATA)) ;
				break ;

			case EH:
				memcpy(&event_header,&data,sizeof(struct DATA)) ;
                Pparticle.evt_no=int(event_header.word[1]) ;
                Pparticle.id=int(event_header.word[2]) ;
                Pparticle.energy=event_header.word[3] ;
                Pparticle.zenith=event_header.word[10] ;
                Pparticle.azimuth=event_header.word[11] ;
				printf("Event No.=%d\n",Pparticle.evt_no) ;
				printf("Primary particle ID=%d\n",Pparticle.id) ;
				printf("Energy=%e eV\n",Pparticle.energy*1.e9) ;
				printf("Zenith=%f radian = %f deg\n",Pparticle.zenith,Pparticle.zenith*180.0/pi) ;
				printf("Azimuth=%f radian = %f deg\n",Pparticle.azimuth,Pparticle.azimuth*180.0/pi) ;
            fprintf(out,"%5d\t%4d\t%e\t%e\t%e\n",Pparticle.evt_no,Pparticle.id,Pparticle.energy,Pparticle.zenith,Pparticle.azimuth) ;
      
                break ;

            case EE:
                fclose(out) ;
                memcpy(&event_end,&data,sizeof(struct DATA)) ;
                break ;

            case RE:
				memcpy(&run_end,&data,sizeof(struct DATA)) ;
				break ;

			default:
				for(int j=0;j<Words_Per_SBlock;j=j+Words_Per_Particle)
				{
					if(round(data.word[j+0])==0 || (round(data.word[j+0])>=Mu_AddInfo_Low && round(data.word[j+0])<=Mu_AddInfo_High))	//Retaining only particle data
						continue ;
					particle.CORSIKA_id=int(data.word[j+0]) ;
					particle.id=int(data.word[j+0]/1000) ;
					particle.px=data.word[j+1] ;
					particle.py=data.word[j+2] ;
					particle.pz=data.word[j+3] ;
					particle.x=data.word[j+4] ;
					particle.y=data.word[j+5] ;
					particle.t=int(data.word[j+6]) ;
					particle.wt=data.word[j+7] ;
					fprintf(out,"%5d\t%4d\t%e\t%e\t%e\t%e\t%e\t%5d\t%e\n",particle.CORSIKA_id,particle.id,particle.px/1e3,particle.py/1e3,particle.pz/1e3,particle.x,particle.y,particle.t,particle.wt) ;
                }
                break ; 
        }
        subblock++ ;
        if(subblock==round(SBlock_Per_Block))	//End of a block
		{
			fread(&space,sizeof(struct SPACE),1,in) ;	//[Two spaces between blocks
			fread(&space,sizeof(struct SPACE),1,in) ;	//  --------------]
			subblock=0 ;
		}
	}
        printf("Reached end of file.");
        printf("Reached end of file.");	
        fclose(in) ;
}       

int main(int argc,char** argv)
{
   if (argc==3)   // batch mode  
   {
    PROCESS_SHOWER_THIN *process_shower_thin ;
	process_shower_thin=new PROCESS_SHOWER_THIN() ;	//To processing showers data in DATxxxxxx files
    process_shower_thin->Read_Data(argv[1],argv[2]) ; 
   }
}
