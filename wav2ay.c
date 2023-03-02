/***

latest official release at: https://github.com/EdouardBERGE/wav2ay

-----------------------------------------------------------------------------------------------------
This software is using MIT "expat" license

« Copyright © BERGE Edouard (roudoudou)

Permission  is  hereby  granted,  free  of charge,to any person obtaining a copy  of  this  software
and  associated  documentation/source   files   of RASM, to deal in the Software without restriction,
including without limitation the  rights  to  use, copy,   modify,   merge,   publish,    distribute,
sublicense,  and/or  sell  copies of the Software, and  to  permit  persons  to  whom the Software is
furnished  to  do  so,  subject  to  the following conditions:

The above copyright  notice  and  this  permission notice   shall   be  included  in  all  copies  or
substantial portions of the Software.
The   Software   is   provided  "as is",   without warranty   of   any   kind,  express  or  implied,
including  but  not  limited  to the warranties of merchantability,   fitness   for   a    particular
purpose  and  noninfringement.  In  no event shall the  authors  or  copyright  holders be liable for
any  claim, damages  or other  liability,  whether in  an  action  of  contract, tort  or  otherwise,
arising from,  out of  or in connection  with  the software  or  the  use  or  other  dealings in the
Software. »
-----------------------------------------------------------------------------------------------------
***/
#ifdef _WIN32
#define OS_WIN 1
#endif

#ifdef _WIN64
#define OS_WIN 1
#endif

#ifdef OS_WIN
#define _USE_MATH_DEFINES
#include<io.h>
#include<fcntl.h>
#endif

#include<stdlib.h>
#include<string.h>
#include<stdio.h>
#include<errno.h>
#include<math.h>

#define MAXCHANNEL 3

/*************************************************************************
  CSV Export Func
*************************************************************************/
void export_csv(char *filename,double *data,int n)
{
	char temp[1024];
	FILE *f;
	int i;

	f=fopen(filename,"wb");

	for (i=0;i<n;i++) {
		sprintf(temp,"%.1lf\n",data[i]);
		fputs(temp,f);
	}
	fclose(f);
}

/*************************************************************************
  WAV Loading Func
*************************************************************************/
struct s_wav_header {
char ChunkID[4];
unsigned char ChunkSize[4];
char Format[4];
char SubChunk1ID[4];
unsigned char SubChunk1Size[4];
unsigned char AudioFormat[2];
unsigned char NumChannels[2];
unsigned char SampleRate[4];
unsigned char ByteRate[4];
unsigned char BlockAlign[2];
unsigned char BitsPerSample[2];
unsigned char SubChunk2ID[4];
unsigned char SubChunk2Size[4];
};

/*
 * strange WAV header
 *
00000000  52 49 46 46 b2 83 c8 00  57 41 56 45 66 6d 74 20  |RIFF....WAVEfmt |
00000010  12 00 00 00 01 00 01 00  22 56 00 00 44 ac 00 00  |........"V..D...|
00000020  02 00 10 00 00 00 66 61  63 74 04 00 00 00 c0 41  |......fact.....A|
00000030  64 00 64 61 74 61 80 83  c8 00 00 00 00 00 00 00  |d.data..........|
00000040  00 00 00 00 00 00 00 00  00 00 00 00 00 00 00 00  |................|
*/

void push_wav(char *outputfilename, short int *data, int n) {
	struct s_wav_header wav_header={0};
	unsigned int datasize,riffsize;
	unsigned char chunksize[4];
	char RIFF[5];
	FILE *f;
#ifdef OS_WIN
	int sr;
#endif

	datasize=n*2;

	memcpy(wav_header.ChunkID,"RIFF",4);
	riffsize=datasize+16+8+datasize+8;
	wav_header.ChunkSize[0]=riffsize&0xFF;
	wav_header.ChunkSize[1]=(riffsize>>8)&0xFF;
	wav_header.ChunkSize[2]=(riffsize>>16)&0xFF;
	wav_header.ChunkSize[3]=(riffsize>>24)&0xFF;
	memcpy(wav_header.Format,"WAVE",4);
	memcpy(wav_header.SubChunk1ID,"fmt ",4);
	wav_header.SubChunk1Size[0]=16;
	wav_header.AudioFormat[0]=1; // integer
	wav_header.NumChannels[0]=1;
	wav_header.SampleRate[0]=0x44;
	wav_header.SampleRate[1]=0xAC; // 44.1KHz
	wav_header.ByteRate[0]=0x88;
	wav_header.ByteRate[1]=0x58;
	wav_header.ByteRate[2]=0x01; // 88.2Kbio
	wav_header.BlockAlign[0]=2;
	wav_header.BitsPerSample[0]=16;
	memcpy(wav_header.SubChunk2ID,"data",4);
	wav_header.SubChunk2Size[0]=datasize&0xFF;
	wav_header.SubChunk2Size[1]=(datasize>>8)&0xFF;
	wav_header.SubChunk2Size[2]=(datasize>>16)&0xFF;
	wav_header.SubChunk2Size[3]=(datasize>>24)&0xFF;

	// dataaaaa output
	f=fopen(outputfilename,"wb");
#ifdef OS_WIN
	sr=_setmode(_fileno(f), _O_BINARY );
	if (sr==-1) {
		fprintf(stderr,"wavoutput windows binary mode got problem!\n");
		exit(1);
	}
#endif
	if (!f) {
		fprintf(stderr,"wavoutput got problem!\n");
		exit(1);
	}
	fwrite(&wav_header,1,sizeof(struct s_wav_header),f);
	fwrite(data,1,datasize,f);
	fclose(f);
}

double (*_internal_getsample)(unsigned char *data, int *idx);

double __internal_getsample8(unsigned char *data, int *idx) {
        double v;
        v=data[*idx]-128;*idx=*idx+1;return v;
}
double __internal_getsample16(unsigned char *data, int *idx) {
	char *most;
        double v;
	most=(char*)data;
        v=most[*idx+1];v+=data[*idx]/256.0;*idx=*idx+2;return v;
}

double *load_wav(char *filename, int *n, double *acqui) {
	struct s_wav_header *wav_header;
	int controlsize,nbchannel,wFormat,frequency,bitspersample,nbsample;
        unsigned char *subchunk;
	unsigned char *data;
        int subchunksize;
	int filesize,idx;
	int i,num;
	FILE *f;
	double *wav;

	f=fopen(filename,"rb");
	if (!f) {
		fprintf(stderr,"file [%s] not found\n",filename);
		return NULL;
	}
	fseek(f,0,SEEK_END);
	filesize=ftell(f);
	fseek(f,0,SEEK_SET);

        if (filesize<sizeof(struct s_wav_header)) {
                fprintf(stderr,"WAV import - this file is too small to be a valid WAV!\n");
                return NULL;
        }

	data=(unsigned char *)malloc(filesize);
	fread(data,1,filesize,f);
	fclose(f);

        wav_header=(struct s_wav_header *)data;
        if (strncmp(wav_header->Format,"WAVE",4)) {
                fprintf(stderr,"WAV import - unsupported audio sample type (format must be 'WAVE')\n");
		free(data);
                return NULL;
        }
        controlsize=wav_header->SubChunk1Size[0]+wav_header->SubChunk1Size[1]*256+wav_header->SubChunk1Size[2]*65536+wav_header->SubChunk1Size[3]*256*65536;
        if (controlsize!=16) {
                fprintf(stderr,"WAV import - invalid wav chunk size (subchunk1 control) => %d!=16\n",controlsize);
        }
        if (strncmp(wav_header->SubChunk1ID,"fmt",3)) {
                fprintf(stderr,"WAV import - unsupported audio sample type (subchunk1id must be 'fmt')\n");
		free(data);
                return NULL;
        }

	if (controlsize==16) {
		subchunk=(unsigned char *)&wav_header->SubChunk2ID;
	} else {
		subchunk=(unsigned char *)&wav_header->SubChunk2ID-16+controlsize;
	}
	// tant qu on n est pas sur le chunk data, on avance dans le "TIFF"
	while (strncmp((char *)subchunk,"data",4)) {
		subchunksize=8+subchunk[4]+subchunk[5]*256+subchunk[6]*65536+subchunk[7]*256*65536;
		if (subchunksize>=filesize) {
			fprintf(stderr,"WAV import - data subchunk not found\n");
			free(data);
			return NULL;
		}
		subchunk+=subchunksize;
	}
	subchunksize=subchunk[4]+subchunk[5]*256+subchunk[6]*65536+subchunk[7]*256*65536;
	controlsize=subchunksize; // taille des samples

        nbchannel=wav_header->NumChannels[0]+wav_header->NumChannels[1]*256;
        if (nbchannel<1) {
                fprintf(stderr,"WAV import - invalid number of audio channel\n");
		free(data);
                return NULL;
        }

        wFormat=wav_header->AudioFormat[0]+wav_header->AudioFormat[1]*256;
        if (wFormat!=1) {
                fprintf(stderr,"WAV import - invalid or unsupported wFormatTag (%04X)\n",wFormat);
		free(data);
                return NULL;
        }

        frequency=wav_header->SampleRate[0]+wav_header->SampleRate[1]*256+wav_header->SampleRate[2]*65536+wav_header->SampleRate[3]*256*65536;
        bitspersample=wav_header->BitsPerSample[0]+wav_header->BitsPerSample[1]*256;

	switch (bitspersample) {
		case 8:_internal_getsample=__internal_getsample8;break;
		case 16:_internal_getsample=__internal_getsample16;break;
		default:fprintf(stderr,"unsupported bits per sample size %d\n",bitspersample); free(data);return NULL;
	}

        nbsample=controlsize/nbchannel/(bitspersample/8);
        if (controlsize+sizeof(struct s_wav_header)>filesize) {
                fprintf(stderr,"WAV import - cannot read %d byte%s of audio whereas the file is %d bytes big!\n",controlsize,controlsize>1?"s":"",filesize);
		free(data);
                return NULL;
        }

	*n=nbsample;

	fprintf(stderr,"nbchannel=%d | bitpersample=%d | Format=%s | frequency=%d\n",nbchannel,bitspersample,wFormat==1?"PCM":"IEEE Float",frequency);

	*acqui=frequency;
	wav=(double *)malloc(nbsample*sizeof(double));
	idx=subchunk+8-data;

       for (i=0;i<nbsample;i++) {
		/* downmixing */
		double accumulator=0.0;
		for (num=0;num<nbchannel;num++) {
			accumulator+=_internal_getsample(data,&idx);
		}
		wav[i]=accumulator/nbchannel;
	}

	return wav;
}

/*****************************************************
 *
 *          Fourier stuff (basic transform)
 *
*****************************************************/

/* slow but understandable Fourier transform */
void calcule_fourier(double *data, int n, double *dataout, int clean) {
	double rel,img,inverse;
	int i,j;

	inverse=1.0/n;

	// skip everything too low according to lowcut filter
	for (i=0;i<clean;i++) dataout[i]=0.0;

	for (;i<=n/2;i++) {
		rel=img=0.0;
		// cumul
		for (j=0;j<n;j++) {
			rel+=data[j]*cos((2*M_PI*j*i)/n);
			img+=data[j]*sin((2*M_PI*j*i)/n);
		}
		// norme
		rel*=inverse;
		img*=inverse;
		dataout[i]=2.0*sqrt(rel*rel+img*img);
	}
}

/* partial Fourier trick to compute only peak neighborhood */
double calcule_fourier_precision(double *data, int n, int peak, double cutlow, double peakHz) {
        double v[16384],rel,img,vmax,inverse;
        int i,j,imax;
	int lp,mp,ws,bornemin;

	bornemin=(double)peak*cutlow/peakHz*16384/n;
	if (bornemin<1) bornemin=1;

	// [lp:mp] let us zoom on previous (small) fourier buffer
	ws=16384/n+1;
	lp=peak*16384/n-ws/2;
	if (lp<bornemin) lp=bornemin; // may happend with strong bass

	mp=peak*16384/n+ws/2;
	if (mp>=8192) mp=8191;

	inverse=1.0/n;

	// do Fourier on that segment
        for (i=lp;i<=mp;i++) {
                rel=img=0.0;
                for (j=0;j<n;j++) {
                        rel+=data[j]*cos((2*M_PI*j*i)/16384.0);
                        img+=data[j]*sin((2*M_PI*j*i)/16384.0);
                }
		// energy
		rel*=inverse;
		img*=inverse;
		v[i]=2.0*sqrt(rel*rel+img*img);
        }

	// max is max
	vmax=v[lp];imax=lp;
	for (i=lp+1;i<=mp;i++) {
		if (v[i]>vmax) {
			vmax=v[i];
			imax=i;
		}
	}

	// precise peak
	return imax*n/16384.0;
}


void passe_bande(double acqui,double basse,double haute,double *datain,double *dataout,int nb)
{
        double *datatmp=NULL;
        double A1,A2,B0,B1,B2;
        double a,b;
        int i;
        /***/
        double e,old_e,old_old_e;
        double s,old_s,old_old_s;

        if (nb<=0) return;
        if (basse<0.0001) basse=0.0001;
        if (haute<0.0001) haute=0.0001;
        if (basse>acqui) basse=acqui;
        if (haute>acqui) haute=acqui;

        datatmp=(double *)malloc(sizeof(double)*nb);
        if (!datatmp) return;

        a=M_PI*basse/acqui;
        b=M_PI*haute/acqui;
        a=sin(a)/cos(a);
        b=sin(b)/cos(b);
        B0= -b / ((1 + a) * (1 + b));
        B1= 0;
        B2= b / ((1 + a) * (1 + b));
        A1= ((1 + a) * (1 - b) + (1 - a) * (1 + b)) / ((1 + a) * (1 + b));
        A2= -(1 - a) * (1 - b) / ((1 + a) * (1 + b));

        /* init */
        old_e=old_old_e=datain[0];
        old_s=old_old_s=0;

	/* a filter phases out signal */
        for (i=0;i<nb;i++) {
                e=datain[i];
                s=e*B0+old_e*B1+old_old_e*B2+old_s*A1+old_old_s*A2;
                datatmp[i]=s;
                old_old_s=old_s;
                old_s=s;
                old_old_e=old_e;
                old_e=e;
        }
        old_e=old_old_e=e;
        old_s=old_old_s=0;

	/* so we filter again reverse */
        while (i>0) {
                i--;
                e=datatmp[i];
                s=e*B0+old_e*B1+old_old_e*B2+old_s*A1+old_old_s*A2;
                dataout[i]=s;
                old_old_s=old_s;
                old_s=s;
                old_old_e=old_e;
                old_e=e;
        }
        free(datatmp);
}


// AY period related to a frequency
struct s_ay_period {
	double *fourier;
};

double * compute_AY(int ws, int replay, int clean, double cuthigh, double cutlow,int periode, double workingfreq) {
	double *retfour;
	double *minisignal;
	double freq,target;
	int i,j,k;
	double vtic,tic,mtic;

	freq=replay*ws; // echantillons/seconde

	minisignal=(double *)malloc(ws*sizeof(double));

	i=periode;

	// generate AY signal for this period, windows size and replay freq
	if (!i) target=workingfreq; else target=workingfreq/i;

	// keep shanon compliant
	if (target*2>freq) {
		retfour=(double *)malloc(ws*sizeof(double));
		memset(retfour,0,ws*sizeof(double));
		free(minisignal);
		return retfour;
	}

	tic=mtic=0.5*freq/target;
	vtic=80.0;

	for (j=0;j<ws;j++) {
		minisignal[j]=vtic;
		tic-=1.0;
		if (tic<=0.0) {
			vtic=-vtic;
			tic+=mtic;
		}
	}
	retfour=(double *)malloc(ws*sizeof(double));
	calcule_fourier(minisignal,ws,retfour,clean);

	free(minisignal);

	return retfour;
}

int getpeak(double *fourier, int ws) {
	double vmax;
	int imax;
	int j;

	imax=0;
	vmax=fourier[0];
	for (j=1;j<=ws/2;j++) {
		if (fourier[j]>vmax) {
			vmax=fourier[j];
			imax=j;
		}
	}
	return imax;
}

int getvolume(double level) {
	int volume;

	if (level>=64) volume=15; else
	if (level>=45) volume=14; else
	if (level>=32) volume=13; else
	if (level>=22) volume=12; else
	if (level>=16) volume=11; else
	if (level>=11) volume=10; else
	if (level>=8) volume=9; else
	if (level>=6) volume=8; else
	if (level>=4) volume=7; else
	if (level>=3.3) volume=6; else
	if (level>=2) volume=5; else
	if (level>=1.7) volume=4; else
	if (level>=1) volume=3; else
	if (level>=0.7) volume=2; else
		volume=1;
	return volume;
}

void do_sample(double *data,int n, double pw, double cutlow, double cuthigh, double acqui, double replay, double preamp, int info, double workingfreq, int nbchannel, double treshold, int dmalist, char *wavout_filename, int *channel_list, int cpclist) {
	double *fourier,*oldfourier;
	double *newdata,subcoef;
	double vmax,resolution,picfreq;
	int nbwin,i,j,k,imax,ws,clean;
	struct s_ay_period *ay=NULL;
	short int *wavout;
	int wavout_n=0;

	// Amstrad registers
	int psgperiod,psgvolume;
	int nbpausedma;
	int channel;

	int AYprevperiod[MAXCHANNEL];
	int AYprevvolume[MAXCHANNEL];
	int AYperiod[MAXCHANNEL];
	int AYvolume[MAXCHANNEL];
	int nbchanges;
	int makenoise=0;
	int hadnoise=0;

	if (pw<0.0 || pw>0.75) {
		fprintf(stderr,"previous weight not in [0.0:0.75] interval. Default value is 0.25\n");
		pw=0.25;
	}
	if (acqui<4000 || acqui>44100) {
		fprintf(stderr,"wrong acquisition frequency. Default value is 15.6KHz\n");
		acqui=15600.0;
	}
	if (cuthigh<cutlow || cutlow<0.0 || cuthigh>acqui || cutlow<replay) {
		fprintf(stderr,"wrong band-pass filter frequencies. Default value are replay-2500Hz\n");
		cuthigh=2500.0;
		cutlow=replay;
	}
	if (cutlow<workingfreq/4095.0*1.3) cutlow=workingfreq/4095.0*1.3;

	ws=acqui/replay; // 312 samples pour 15000Hz avec replay a 50Hz
	if (ws<200) {
		fprintf(stderr,"window size is small, results may be very innacurate\n");
	}
	if (ws<100) {
		fprintf(stderr,"wont compute with such a small window, raise input frequency\n");
		return;
	}

	// compléter le nombre de sample pour arriver pile sur un multiple de ws (window size) | default:312
	nbwin=n/ws;
	while (nbwin*ws<n) {
		nbwin++;
	}
	if (nbwin*ws!=n) {
		data=(double *)realloc(data,sizeof(double)*nbwin*ws);
		for (i=n;i<nbwin*ws;i++) {
			data[i]=data[i-1]*0.95; // histoire de terminer salement mais en evitant le poc => normalement le sample est ok de base!
		}
		n=nbwin*ws;
	}

	newdata=(double *)malloc(sizeof(double)*n);

	// calculer la valeur du nettoyage a faire sur la transformee pour eviter les basses frequences mal filtrees
	resolution=acqui/(double)ws;
	clean=1+floor(replay/resolution+0.5);

	// precalc des fouriers de l'AY sur la fenetre utilisee
	ay=(struct s_ay_period *)malloc(4096*sizeof(struct s_ay_period));
	memset(ay,0,4096*sizeof(struct s_ay_period));

	// par defaut il faut filtrer en dessous de la frequence de replay et au dessus de 2500Hz a cause de la precision de restitution
	fprintf(stderr,"band-pass %.1lf|%.1lf window size=%d preamp=%.1lf AY frequency=%.1lfMHz\n",cutlow,cuthigh,ws,preamp,workingfreq*16.0);
	passe_bande(acqui,cutlow,cuthigh,data,newdata,n);
	memcpy(data,newdata,n*sizeof(double));

	// init temporal buffer	
	fourier=(double *)malloc(sizeof(double)*ws);
	oldfourier=(double *)malloc(sizeof(double)*ws);
	calcule_fourier(data,ws,oldfourier,clean);

	fprintf(stderr,"%d windows for a %.1lfs sample\n",nbwin,n/acqui);

	for (channel=0;channel<nbchannel;channel++) {
		AYprevperiod[channel]=0xFFFF;
		AYprevvolume[channel]=0;
	}

	if (!dmalist && !cpclist) {
		printf("idx");
		for (channel=0;channel<nbchannel;channel++) {
			printf(";reg;value");
			printf(";reg;value");
			printf(";reg;value");
		}
		printf("\n");
	}

	if (wavout_filename) {
		wavout=(short int *)malloc(nbwin*(sizeof(short int)*workingfreq*2.0/replay+1));
		wavout_n=0;
	}

	if (preamp!=1.0) for (j=0;j<n;j++) data[j]*=preamp;

	for (i=0;i<nbwin;i++) {

		// calcul de fourier sur le segment
		calcule_fourier(&data[ws*i],ws,fourier,clean);

		for (channel=0;channel<nbchannel;channel++) {
			/*
			// mix avec le segment precedent si besoin
			if (pw>0.0) {
				for (j=0;j<ws/2;j++) {
					fourier[j]=oldfourier[j]*pw+fourier[j]*(1.0-pw);
				}
			}
			*/

			// trouver le pic
			imax=getpeak(fourier,ws);
			// energie suffisante pour s'en soucier?
			if (fourier[imax]>treshold) {
				int volume;
				int plow,phigh;

				if (info) fprintf(stderr,"segment %3d pic en %.0lfHz (resolution=%.1lf) norme=%.1lf\n",i,imax*resolution,resolution,fourier[imax]);

				AYvolume[channel]=getvolume(fourier[imax]);
				k=imax;

//export_csv("fourier_ref.csv",fourier,ws/2+1);

				// chercher la precision sur le pic
				picfreq=calcule_fourier_precision(&data[ws*i],ws,imax,cutlow,imax*resolution)*resolution;
				if (info) fprintf(stderr," => recherche de precision donne %.1lfHz\n",picfreq);

				if (picfreq>2500) makenoise++;

				// appliquer la soustraction AY sur la transformee de Fourier pour la suite
				imax=workingfreq/picfreq; // periode AY
				if (imax>4095) imax=4095;

				if (!ay[imax].fourier) ay[imax].fourier=compute_AY(ws,replay,clean,cuthigh,cutlow,imax,workingfreq);

				subcoef=fourier[k]/ay[imax].fourier[k];
				for (j=0;j<=ws/2;j++) fourier[j]-=ay[imax].fourier[j]*subcoef;

//export_csv("fourier_ay.csv",ay[imax].fourier,ws/2+1);
//export_csv("fourier_sortie.csv",fourier,ws/2+1);
//exit(1);

				// produire les registres AY (imax=periode)
				AYperiod[channel]=imax;
			} else {
				AYvolume[channel]=0;
				/* volume a zero on ne touchera pas la periode! */
				AYperiod[channel]=AYprevperiod[channel];
			}
		}
		/*********************************
		       Channel optimisation   
		*********************************/

		// d'abord un tri pour conserver une harmonie éventuelle sur le volume!
		if (nbchannel>1) {
			if (AYperiod[0]<AYperiod[1]) {
				imax=AYperiod[1];
				AYperiod[1]=AYperiod[0];
				AYperiod[0]=imax;
				imax=AYvolume[1];
				AYvolume[1]=AYvolume[0];
				AYvolume[0]=imax;
			}
			if (nbchannel>2) {
				if (AYperiod[0]<AYperiod[2]) {
					imax=AYperiod[2];
					AYperiod[2]=AYperiod[0];
					AYperiod[0]=imax;
					imax=AYvolume[2];
					AYvolume[2]=AYvolume[0];
					AYvolume[0]=imax;
				}
				if (AYperiod[1]<AYperiod[2]) {
					imax=AYperiod[2];
					AYperiod[2]=AYperiod[1];
					AYperiod[1]=imax;
					imax=AYvolume[2];
					AYvolume[2]=AYvolume[1];
					AYvolume[1]=imax;
				}
			}
		}
		/* switch channel if period are exactly the same */
		if (nbchannel>1) {
			if (AYperiod[0]==AYprevperiod[1]) {
				imax=AYperiod[1];
				AYperiod[1]=AYperiod[0];
				AYperiod[0]=imax;
				imax=AYvolume[1];
				AYvolume[1]=AYvolume[0];
				AYvolume[0]=imax;
			}
			if (nbchannel>2) {
				if (AYperiod[0]==AYprevperiod[2]) {
					imax=AYperiod[2];
					AYperiod[2]=AYperiod[0];
					AYperiod[0]=imax;
					imax=AYvolume[2];
					AYvolume[2]=AYvolume[0];
					AYvolume[0]=imax;
				} else if (AYperiod[1]==AYprevperiod[2]) {
					imax=AYperiod[2];
					AYperiod[2]=AYperiod[1];
					AYperiod[1]=imax;
					imax=AYvolume[2];
					AYvolume[2]=AYvolume[1];
					AYvolume[1]=imax;
				}
			}
		}

		if (wavout_filename) {
			int ref[16]={0,231,695,1158,2084,2779,4168,6716,8105,13200,18294,24315,32189,40757,52799,65535};
			int pulse,maxp;
			int acc;
			int tic[MAXCHANNEL]={0};
			int vchan[MAXCHANNEL];

			maxp=workingfreq*2.0/replay;
			for (channel=0;channel<nbchannel;channel++) {
				vchan[channel]=ref[AYvolume[channel]]/2;
				tic[channel]=AYperiod[channel];
			}

			for (pulse=0;pulse<maxp;pulse++) {
				acc=0;
				for (channel=0;channel<nbchannel;channel++) {
					acc+=vchan[channel];
					tic[channel]--;
					if (!tic[channel]) {
						tic[channel]=AYperiod[channel];
						vchan[channel]=-vchan[channel];
					}
				}
				wavout[wavout_n++]=acc/nbchannel;
			}
		}

		/*********************************
		      Amstrad CPC list output
		*********************************/
		if (cpclist) {
			int packed_exec=0;
			int minilist[10];
			int ilist=0;
			int zebit=128;

			if (!i) {
				printf("defb %d ; nombre d'iterations\n",nbwin);
			}

			for (channel=0;channel<3;channel++) {
				if (AYvolume[channel]!=AYprevvolume[channel]) {
					packed_exec|=zebit;
					minilist[ilist++]=AYvolume[channel];
				}
				zebit>>=1;
				if ((AYperiod[channel]>>8)!=(AYprevperiod[channel]>>8)) {
					packed_exec|=zebit;
					minilist[ilist++]=AYperiod[channel]>>8;
				}
				zebit>>=1;
				minilist[ilist++]=AYperiod[channel]&0xFF; // au minimum on enverra la frequence "basse"
			}

			if (hadnoise && !makenoise) {
				packed_exec|=2;
				minilist[ilist++]=56;
			}
			if (!hadnoise && makenoise) {
				packed_exec|=2;
				minilist[ilist++]=0; // bruit partout sinon on n'entend rien...
			}
			printf("defb #%02X",packed_exec);
			for (channel=0;channel<ilist;channel++) printf(",#%02X",minilist[channel]);
			printf(" ; %5d-%02d %5d-%02d %5d-%02d\n",AYperiod[0],AYvolume[0],AYperiod[1],AYvolume[1],AYperiod[2],AYvolume[2]);

			hadnoise=makenoise;
			makenoise=0;

			for (channel=0;channel<nbchannel;channel++) {
				AYprevperiod[channel]=AYperiod[channel];
				AYprevvolume[channel]=AYvolume[channel];
			}
		} else
		/*********************************
		      Amstrad Plus DMA output
		*********************************/
		if (dmalist) {
			psgperiod=0;
			psgvolume=8;
			nbchanges=1;

			for (channel=0;channel<nbchannel;channel++) {
				if (AYvolume[channel]) {
					if (AYperiod[channel]!=AYprevperiod[channel]) {
						if ((AYperiod[channel]>>8)!=(AYprevperiod[channel]>>8)) {
							printf("defb %d,%d : ",AYperiod[channel]>>8,psgperiod+1);
							nbchanges++;
						}
						if ((AYperiod[channel]&0xFF)!=(AYprevperiod[channel]&0xFF)) {
							printf("defb %d,%d : ",AYperiod[channel]&0xFF,psgperiod);
							nbchanges++;
						}
					}
				}
				if (AYvolume[channel]!=AYprevvolume[channel]) {
					printf("defb %d,%d : ",AYvolume[channel],psgvolume);
					nbchanges++;
				}
				psgperiod+=2;
				psgvolume+=1;
			}

			nbpausedma=312.0/replay*50.0-nbchanges;
			// no control on this value
			printf("defb %d,#30+%d\n",nbpausedma&0xFF,(nbpausedma>>8)&0xF); // pause 50Hz

			for (channel=0;channel<nbchannel;channel++) {
				AYprevperiod[channel]=AYperiod[channel];
				AYprevvolume[channel]=AYvolume[channel];
			}
		} else {
		/*********************************
		            CSV Export
		*********************************/
			if (!i) {
			}
			printf(";%d",i);
			for (channel=0;channel<nbchannel;channel++) {
				printf(";%d;%d",channel_list[channel]  ,AYvolume[channel]);
				printf(";%d;%d",0+channel*2,AYperiod[channel]&0xFF);
				printf(";%d;%d",1+channel*2,AYperiod[channel]>>8);
			}
			if (hadnoise && !makenoise) {
				printf(";7;0");
			}
			if (!hadnoise && makenoise) {
				printf(";7;40");
			}
			if (makenoise) printf(" // noise");
			hadnoise=makenoise;
			makenoise=0;
			printf("\n");
		}
	}

	if (wavout_filename) {
		double acc,pos,step,aron,frac,integer_part;
		int iidx,oidx=0;
		step=2.0*workingfreq/44100.0;

		for (pos=0;pos<wavout_n-step;pos+=step) {
			frac=modf(pos,&integer_part);
			iidx=integer_part;
			acc=(1.0-frac)*wavout[iidx];
			while (integer_part+1.0<pos+step) {
				acc+=wavout[++iidx];
				integer_part+=1.0;
			}
			frac=modf(pos,&integer_part);
			acc+=frac*wavout[++iidx];
			acc/=step;
			wavout[oidx++]=acc;
		}
		push_wav(wavout_filename,wavout,oidx);
		free(wavout);
	}
}





void usage() {
	printf("=========== conversion from WAV to AY registers ===========\n");
	printf("usage: wav2ay.exe <wavfile> <options>\n");
	printf("options:\n");
	printf("-preamp <value>  amplification | default 1.0\n");
	printf("-tresh  <value>  minimal energy| default 0.25 is minimum (max usable approx 15)\n");
	printf("-replay <value>  frequency play| default 10\n");
	printf("-high   <value>  highcut filter| default 4000\n");
	printf("-low    <value>  lowcut filter | default is replay (min: 20Hz)\n");
//	printf("-pw     <value>  sound inertia | default 0.0 max 0.75\n");
	printf("-wfreq  <value>  AY frequency  | default 1000000 (1MHz)\n");
	printf("-nbchan <value>  nb channel    | default 3\n");
	printf("-chans  <value>  channel used  | default 'ABC'\n");
	printf("-dmalist         output optimised DMA list\n");
	printf("-cpclist         output optimised list for CPC replay\n");
	printf("-wavout <file>   output WAV preview\n");
	printf("-verbose\n");
	printf("\n");
	exit(1);
}

void main(int argc, char **argv) {
	int info,n,i,ifilename=-1,idx;
	double *data;
	// conversion param
	double preamp,replay,acqui,cuthigh,cutlow,pw,workingfreq,treshold;
	int nbchannel,dmalist,cpclist;
	int channel_list[3],ichan=0;
	char *wavoutfilename=NULL;

	dmalist=cpclist=0;
	preamp=1.0;
	replay=10;
	cuthigh=4000;
	cutlow=replay;
	pw=0.0;
	info=0;
	workingfreq=62500.0;
	nbchannel=3;
	treshold=0.5;

	// default channels
	channel_list[0]=8;
	channel_list[1]=9;
	channel_list[2]=10;

	for (i=1;i<argc;i++) {
		if (strcmp(argv[i],"-dmalist")==0) {
			dmalist=1;
		} else if (strcmp(argv[i],"-cpclist")==0) {
			cpclist=1;
		} else if (strcmp(argv[i],"-verbose")==0) {
			info=1;
		} else if (strcmp(argv[i],"-tresh")==0 && i+1<argc) {
			treshold=atof(argv[++i]);
			if (treshold<0.5) treshold=0.5;
		} else if (strcmp(argv[i],"-wavout")==0 && i+1<argc) {
			wavoutfilename=argv[++i];
		} else if (strcmp(argv[i],"-wfreq")==0 && i+1<argc) {
			workingfreq=atof(argv[++i])/16.0;
		} else if (strcmp(argv[i],"-preamp")==0 && i+1<argc) {
			preamp=atof(argv[++i]);
		} else if (strcmp(argv[i],"-chans")==0 && i+1<argc) {
			for (idx=0;argv[i+1][idx];idx++) {
				switch (argv[i+1][idx]) {
					case 'a':case 'A':case '0':channel_list[ichan++]=8;break;
					case 'b':case 'B':case '1':channel_list[ichan++]=9;break;
					case 'c':case 'C':case '2':channel_list[ichan++]=10;break;
					default:
					   fprintf(stderr,"Error defining channels, must be like -chans ABC or -chans 012 or -chans b\n");
					   exit(0);
				}
				if (ichan==3) break; // skip
			}
			i++;
		} else if (strcmp(argv[i],"-nbchan")==0 && i+1<argc) {
			nbchannel=atoi(argv[++i]);
			if (nbchannel<1 || nbchannel>MAXCHANNEL) usage();
		} else if (strcmp(argv[i],"-replay")==0 && i+1<argc) {
			replay=atof(argv[++i]);
		} else if (strcmp(argv[i],"-high")==0 && i+1<argc) {
			cuthigh=atof(argv[++i]);
		} else if (strcmp(argv[i],"-low")==0 && i+1<argc) {
			cutlow=atof(argv[++i]);
		} else if (strcmp(argv[i],"-pw")==0 && i+1<argc) {
			pw=atof(argv[++i]);
		} else if (ifilename==-1) {
			ifilename=i;
		} else usage();
	}

	if (ifilename==-1) usage();

	if ((data=load_wav(argv[ifilename],&n,&acqui))!=NULL) {
		do_sample(data,n,pw,cutlow,cuthigh,acqui,replay,preamp,info,workingfreq,nbchannel,treshold,dmalist,wavoutfilename,channel_list,cpclist);
	}
}

