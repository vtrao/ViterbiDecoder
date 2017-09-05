/*-------------------------------------------------------------------------------------------------------------*/

/* Convolutional Encoder and Viterbi Decoder with Tail Bitting in AWGN Channel(Additive Gaussian Noise Channel)*/

/*-------------------------------------------------------------------------------------------------------------*/

#include<stdlib.h>
#include<time.h>
#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#define PI 3.14159
#define N 20 			/* Maximum output bits to represent k input symbols */
#define L 20 			/* Maximum length of polynomial i.e.., constraint length */

/*-----------------------------Function Prototypes of various modules of Project-------------------------------*/

int cti(char c);
void genranddata(int framesize);
void shiftnencode(int k,int K,int n);
void encode();
float gauss(float mean, float sigma);
void addnoise();
int btd(int *b, int size);
void dtb(int d,int x,int y,int *b);
void extract(int x,int y,int *source,int *destination);
void constructstatetable();
void decode();
float euclidean(int op,float *l);
int statetrans(int i,int j,int *bt);
void calculatebiterrorrate();

/*------------------------------------------------------------------------------------------------------------*/

FILE *of;
int **sim;	                  /*State Input Matrix*/
int **som;	                  /*State Output Matrix*/
int **nsm;	                  /*Next State Matrix*/
int n,k,K,pow1,pow2,limit,dbs,ofs;
int arr[L-1];
int a[N][L];
int shiftreg[L];
float *gaussop,snr=00.0f;
float str[N];
struct staten           /* This is the Data Structure used for decoding representing a node at time instant t*/  
	{
	int prev;
	int value;
	int input;		/*new entry #1*/
	float accum;
	};

/*------------------------------------------------------------------------------------------------------------*/

		
int cti(char c)                /* This Function converts character to integer in various parts of the Project */
	{
	if(c=='1')
		return(1);
	else
		return(0);
	}

/*------------------------------------------Input Generator Starts here--------------------------------------*/

void genranddata(int framesize) /* This is for the generation of random input */
	{
	int i,x;
	static int j1=2;
	FILE *fp;
	fp=fopen("tinput.txt","w");
        srand(j1++);            /* This is to set seed value to different ones for diff input frames so j static*/
       	printf("\n");
        for(i=0;i<framesize;i++)
       	        {
               	x=rand()/999999999;
                if(x==2) x=1;
		fprintf(fp,"%d",x);
		printf("%d",x);
               	}
	fputc(' ',fp);
	fclose(fp);
	return;
	}

/*---------------------------------------------Encoder Starts Here-------------------------------------------*/

void shiftnencode(int k,int K,int n) /* This function is called from encode() for each k bit blocks of input */
	{
	int i,j,m,y;
	int and[K-1];
	printf("\ninput\t");
	for(i=k;i>=0;i--)
		printf("%d",arr[i]);
	printf("\nbefore shifti\t");
	for(i=0;i<=K;i++)
		printf("%d",shiftreg[i]);
	printf("\nafter shift\t");   /* here shifting starts */
	for(i=k;i>=0;i--)
		{
		for(j=K;j>0;j--)
		    shiftreg[j]=shiftreg[j-1];
		shiftreg[0]=arr[i];
		}
	for(i=0;i<=K;i++)
                printf("%d",shiftreg[i]);
	printf("\nThe and Operation yeilds y as \t");
	for(i=0;i<=n;i++)
		{
		for(j=0;j<=K;j++)
			{
			and[j]=shiftreg[j] & a[i][j]; /* And Operation is done here */
			}
		for(m=0;m<=K;m++)
			{
			if(m==0)
				y=and[m];
			else
				y^=and[m];  /* The XOR operation takes place here */
			}
		printf(" %d",y);
		fprintf(of,"%d",y);
		}	
	}

void encode()                   /* This function actually Encodes the input in the input.txt file */
        {
        int arr1[L-1];
        FILE *fp;
        char c;
        int i,j,x,y,m,p;

        fp=fopen("tinput.txt","r");

        fseek(fp,-(K+1),SEEK_END);
        y=K;
        while((c=fgetc(fp))!=EOF)
                {
                shiftreg[y-1]=cti(c);  /* setting the contents of Shift Register */
                --y;
                }
        rewind(fp);
        for(i=0;i<dbs;i+=k)
                {
                for(y=0;y<k;y++)
                        {
                        c=fgetc(fp);
                        arr1[y]=cti(c);
                        }
                for(p=0,m=k-1;p<k;p++,m--)
                        arr[m]=arr1[p];
                shiftnencode(k-1,K-1,n-1); /* Shiftandencode funtion is called for a particular k bit block */
                }
        fclose(fp);
        //constructstatetable();
        }

/*--------------------------Noise Addition Starts here(Transmit over the channel)--------------------------*/

float gauss(float mean, float sigma)    /* This function is called from addnoise to generates gauss output */
	{
        double u,r,g;
        u=(double)rand()/RAND_MAX;
        if(u==1.0)
        u=0.99999;
        r=sigma*sqrt(2*log(1/(1-u)));
        u=(double)rand()/RAND_MAX;
        if(u==1.0)
        u=0.999999;
        g=(float)(mean + r*sin(2*PI*u));
        return(g);
	}

void addnoise()                 /* This function actually adds noise to the encoded.txt file, generates gaussop.txt*/
        {
        FILE *fp,*fof;
        float b,temp;
        int ci;
        char v;
        float sigma,mean=0,sn_r,es=1;
        sn_r=(float)pow(10,(snr/10));
        sigma=(float)sqrt((es/(2*sn_r)));
        if((fp=fopen("encoded.txt","r"))==NULL)
                {
                printf("\n file does not exist\n");
                exit(0);
                }
        fof=fopen("gaussop.txt","w");
        while((v=fgetc(fp))!=EOF)
                {
                ci=cti(v);
                b=2*ci-1;
                fprintf(fof,"%f ",b+gauss(mean,sigma));
		printf("\t%f",b+gauss(mean,sigma));
                }
        fclose(fp);
        fclose(fof);
        }

int btd(int *b, int size)       /* This is for binary to Decimal calculation */
	{
    	int i, d;
        d = 0;
        for (i = 0; i < size; i++)
    	   	d += b[i] << (size - i - 1);
    	return(d);
	}

void dtb(int d,int x,int y,int *b)  /* This is for decimal to binary conversion */
	{
	int i;
    	for(i = x; i < y; i++)
        	b[i] = 0;
    	b[y - 1] = d & 0x01;
    	for (i = y - 2; i >= x; i--) 
		{
        	d = d >> 1;
        	b[i] = d & 0x01;
    		}
	}

/*---------------------------Data Structure Construction for Decoder Starts Here-----------------------------*/

void extract(int x,int y,int *source,int *destination) /* This for extraction of particular bits from source */
	{
	int v7,v8;
	for(v7=x,v8=0;v7<=y;v7++,v8++)
		{
		destination[v8]=source[v7];
		}	
	}

void constructstatetable()       /* This is to construct various tables required to decode using k,n & K */ 
        {
	int x,y,i,j,l,w,v1,v2,v3,v4,v6,v5;
	int sr[K],and1[K],tempi[n],state[K-k];
        x=pow(2,K-k); 
        y=pow(2,k);  
        sim=(int **)malloc(x*(sizeof(int**)));
        som=(int **)malloc(x*(sizeof(int**)));
        nsm=(int **)malloc(x*(sizeof(int**)));
        for(i=0;i<x;i++)
                {
                sim[i]=(int*)malloc(x*sizeof(int*));
                som[i]=(int*)malloc(y*sizeof(int*));
                nsm[i]=(int*)malloc(y*sizeof(int*));
		for(l=0;l<x;l++)		
			sim[i][l]=-1;
		}
	for(i=0;i<x;i++)
		{
		for(j=0;j<y;j++)
			{
			dtb(i,k,K,sr);
			dtb(j,0,k,sr);
			for(v1=0;v1<n;v1++)
				{
				for(v2=0;v2<K;v2++)
					and1[v2] = sr[v2] & a[v1][v2];
				for(v3=0;v3<K;v3++)
					if(v3==0) v4=and1[v3];
					else v4^=and1[v3];
				tempi[v1]=v4;
				}
			som[i][j]=btd(tempi,n);       /* state Output matrix is calculated here */
			for(v5=K-k,v6=1;v5>=0;v5--,v6++)
                                sr[K-v6]=sr[v5];

                        extract(k-1,K-1,sr,state);
                        nsm[i][j]=btd(state,K-k);     /* Next State matrix is calculated here */
			}
		}
	for(i=0;i<x;i++)
		{
		for(j=0;j<y;j++)
			{
			for(v1=0;v1<x;v1++)
				{
				if(nsm[i][j]==v1)
					sim[i][v1]=j; /* State Input matrix is calculated here */ 
				}
			}
		}
	printf("\nThe State Output Matrix\n"); 
        for(i=0;i<x;i++)
                {
		printf("\n");
		for(j=0;j<y;j++)
			{
			printf(" %6d",som[i][j]);
			}
                }
	printf("\nThe Next State Matrix\n");
	for(i=0;i<x;i++)
                {
                printf("\n");
                for(j=0;j<y;j++)
                        {
                        printf(" %6d",nsm[i][j]);
                        }
                }
	printf("\nThe State Input Matrix\n");
	for(i=0;i<x;i++)
                {
                printf("\n");
                for(j=0;j<x;j++)
                        {
                        printf(" %4d",sim[i][j]);
                        }
                }
        }

/*------------------------------------------THE DECODER STARTS HERE-------------------------------------------*/

float euclidean(int op,float *l)  /* This is to calculate the free distance between the given outputs */
        {
        float h=0.00f;
        int a[n];
        int i=0;
        dtb(op,0,n,a);
        for(i=0;i<n;i++)
                {
                a[i]=2*a[i]-1;//printf("\t%d",a[i]);
                h+=pow((l[i]-a[i]),2);
                }
        h=sqrt(h);
        return(h);
        }

int statetrans(int i,int j,int *bt) /* This is to calculate minimum distance input between the given twostates*/
        {
        int m,x1,y1,flag3=0;
	float min,euc;
        m=sim[i][j];
        if(m==-1)
                return(m);
        else
                {
		for(x1=0;x1<pow2;x1++)
			{
			if(nsm[i][x1]==j)
				{
				y1=som[i][x1];
				euc=euclidean(y1,str);
				printf("\niBest Input %f %d",euc,x1);
				if(flag3==0)
					{
					min=euc;
					*bt=x1;
					flag3=1;
					}
				else
					{
					if(min>=euc)
						{
						*bt=x1;
						min=euc;
						}
					}
				}
			}
                return(som[i][*bt]);
                }
        }

void decode()                      /* This is the actually Decoder of the Project */
        {
	struct vv
		{
		float accum;
		int ssa[limit+1];
		}states[pow1];
        FILE *fp,*dp;
        struct staten stprev[pow1][limit+1];
        int cur,t,x1,x2,z=0,i=0,j=0;
        int bestinputb,bestinputa,bestinput;                            /* new entry #1*/
        int op,flag=0,flag1=0,minst,*vvv,q1=0,q2=0;
	vvv=(int *)malloc(k*sizeof(int));
        float wt,mineu,fl;
        for (cur=0;cur<pow1;cur++)
		{
        	printf("\ncur %d here i am",cur);
        	fp=fopen("gaussop.txt","r");
	        for(i=0;i<n;i++)
        	        {
                	fscanf(fp,"%f",&str[i]);
        	        }
	        for(t=1;t<limit+1;t++)
        	        {
                	printf("\n time instant\n");
	                for(x1=0;x1<pow1;x1++)
        	                {
                	        if(z==0)
                        	        {
                                	op=statetrans(cur,x1,&bestinput);
	                                if(op==-1)
        	                                {
                	                        stprev[x1][t].value=-1;
                        	                continue;
                                	        }
	                                bestinputa=bestinput;   /*new entry #1*/
        	                        wt=euclidean(op,str);
                	                stprev[x1][t].accum=wt;printf("\n%d %d %f %d",t,x1,wt,bestinputa);
                        	        stprev[x1][t].value=1;
                                	stprev[x1][t].input=bestinputa;   /*new entry #1*/
	                                stprev[x1][t].prev=cur;
        	                        }
                	        else
                        	        {
                                	flag=0;
	                                for(x2=0;x2<pow1;x2++)
        	                                {
                	                        if(stprev[x2][t-1].value == -1)
                        	                        continue;
                                	        else
                                        	        {
                                                	op=statetrans(x2,x1,&bestinput);
                                                	if(op==-1)
                                                        	continue;
		                                        flag1=1;          /*to determine path exixtence between i and j*/
                	                                bestinputa=bestinput;   /*new entry #1*/
	                                                wt=euclidean(op,str)+stprev[x2][t-1].accum;
        	                                        printf("\n%d %d %d %f %d",t,x1,x2,wt,bestinputa);
                	                                if(flag==0)
                        	                                {
                                	                        mineu=wt;
                                        	                bestinputb=bestinputa;  /*new entry #1*/
                                                	        minst=x2;
                                                        	flag=1;
	                                                        }
        	                                        else
                	                                        {
                        	                               if(mineu>=wt)
                                	                                {
                                        	                        bestinputb=bestinputa; /*new entry #1*/
                                                	                mineu=wt;
                                                        	        minst=x2;
                                                                	}
                                                        	}

                                                	}
                                        	}               /*x2 loop ends here*/
                                	if(flag1==1)
                                        	{
                                        	stprev[x1][t].prev= minst;
                                        	stprev[x1][t].value= 1;
                                        	stprev[x1][t].input=bestinputb; /*new entry #1*/
	                                        stprev[x1][t].accum= mineu;   //stprev[minst][t-1].accum + mineu;
        		                        printf("\n%d %d %f %f %d %d",t,x1,stprev[x1][t].accum,mineu,stprev[x1][t].prev,bestinputb);
                        	                }
                                	else
                                        	stprev[x1][t].value=-1;
                                	flag1=0;
                                	}
                        	printf("\n");
                        	} /*x1 loop ends here */
                        	z=1;
                		for(i=0;i<n;i++)
                        		{
                        		fscanf(fp,"%f",&str[i]);
                        		}
                	} /* t loop ends here */
        	t=t-1;
		states[cur].accum==stprev[cur][t].accum;;
		states[cur].ssa[limit]=cur;
		states[cur].ssa[0]=cur;
		for(i=limit-1;i>=0;i--)
	                {
	                states[cur].ssa[i]=stprev[states[cur].ssa[i+1]][i+1].prev;
        	        }
		fclose(fp);
/*		for(q1=0;q1<pow1;q1++)
			for(q2=0;q2<limit+1;q2++)
				{
				stprev[q1][q2].value=0;
				stprev[q1][q2].prev=0;
				stprev[q1][q2].input=0;
				stprev[q1][q2].accum=0;
				}*/
		}
	flag=0;
        for(i=0;i<pow1;i++)
                {
                fl=states[i].accum; 
                if(flag==0)
                	{
                        mineu=states[i].accum;
                        minst=i;
                        flag=1;
                        }
                else
                	if(mineu>=fl)
                        	{
                                mineu=states[i].accum;
                                minst=i;
                                }
                        
                }

        dp=fopen("decoded.txt","a");
	//fputc('\n',fp);
        for(i=0;i<=limit;i++)
                {
                printf("\t%d",states[minst].ssa[i]);
                }
        printf("\n");
        for(i=1;i<=limit;i++)
                {
                x1=stprev[ states[minst].ssa[i] ] [i].input;
                dtb(x1,0,k,vvv);
                for(z=k-1;z>=0;z--)
                        {
                        fprintf(dp,"%d",vvv[z]);
                        printf("%d",vvv[z]);
                        }
                }
	fclose(dp);
        }

/*------------------------------------------Bit Error Rate Calculator-----------------------------------*/

void calculatebiterrorrate()
	{
	FILE *fp,*fp1,*lf;
	int e=0,i=0;
	char c1,c2;
	lf=fopen("logfile.txt","a");
	fp=fopen("decoded.txt","r");
	fp1=fopen("input.txt","r");
	while(feof(fp)==0)
		{
		i++;
		c1=fgetc(fp);
		c2=fgetc(fp1);
		if(c1==c2)
			{
			}
		else
			e++;
		}
	printf("\nThe error rate is %f%%\n",((float)e/i)*100);
	fprintf(lf,"\nThe error rate is %f%%\n",((float)e/i)*100);
	fclose(fp);
	fclose(fp1);
	fclose(lf);
	}

/*---------------------------------------THE MAIN FUNCTION STARTS HERE-------------------------------------*/

int main()
	{
	FILE *f1,*fp,*lf;
	char c;
        int i,j,nof,y;   /*ofs-Output Frame Size, ifs-Input Frame Size and nof-Number Of Frames to consider*/
	remove("decoded.txt");
	remove("input.txt");
	lf=fopen("logfile.txt","a");
	fprintf(lf,"\n\t---------------Latest Transmission if this is the last Entry------------------\t\n");
        printf("\nEnter the Rate of Encoder\n");
        inputvalidation1:

                printf("\n Enter the k and n values of (k,n) code\n");
                scanf("%d%d",&k,&n);
                if(k>=n)
                        {
                        printf("Invalid input k cant be Greater than or Equal to  n\n Re Enter");
                        goto inputvalidation1;
                        }
		fprintf(lf,"\nThe Rate of Encoder  \n\tk is %d\n\tn is %d\n",k,n);  

        inputvalidation2:

                printf("\nEnter the constraint length K\n");
                scanf("%d",&K);
                if(K<=k)
                        {
                        printf("Invalid input constraint length cant be less than k of (k,n)\nRe Enter");
                        goto inputvalidation2;
                        }
		fprintf(lf,"\nThe Constraint Length K is %d\n",K);

	for(i=0;i<n;i++)
                {
                printf("Enter the %d polynomial of %d polynomials \n",i+1,n);
		fprintf(lf,"\nThe %d polynomial of %d polynomials\n",i+1,n);
                for(j=0;j<K;j++)
                        {
                        scanf("%d",&a[i][j]);
                        if(a[i][j]==0 || a[i][j]==1)
				{
				fprintf(lf,"\t%d ",a[i][j]);
                                continue;
				}
                        else
                                {
                                j=j-1;
                                printf("input symbol cant be anything other than 0 or 1\n");
                                }
                        }
                }

	Inputvalidation3:
	        
		printf("\nEnter the Packet Size or OutPut Frame Size\n");
	        scanf("%d",&ofs);  /* ofs is x in our analysis*/
		if(ofs%n!=0)
			{
			printf("\nEnter the Packet Size multiple of n\n");
			goto Inputvalidation3;
			} 
	fprintf(lf,"\nThe Packet Size of the Transmission Channel is %d \n",ofs);
        printf("\nEnter the Number of (Packets)Frames to Transmit\n");
       	scanf("%d",&nof);  /* ifs is y in our analysis*/
	fprintf(lf,"\nNumber of Frames Transmitted are %d \n",nof);
	fprintf(lf,"\nThe Signal to Noise Ratio (SNR) is %fdb\n",snr);
	fclose(lf);
	y=(ofs*k)%n;
	dbs=y+ofs*k/n;
	y=dbs;

	limit= dbs/k;	
        pow1=pow(2,K-k);
        pow2=pow(2,k);
//	genranddata(dbs);       /*  This is for Analysis sake */
	constructstatetable();
	printf("\n%d %d %d %d",dbs,limit,pow1,pow2);
	for(i=0;i<nof;i++)
		{
		of=fopen("encoded.txt","w");
		genranddata(dbs);
		printf("\n");
		encode();
		fclose(of);
		addnoise();
		decode();
		fp=fopen("tinput.txt","r");
		f1=fopen("input.txt","a");
		//fputc('\n',f1);
		while(y!=0)
			{
			c=fgetc(fp);
			fprintf(f1,"%c",c);
			y--;	
			}
		y=dbs;
		fclose(fp);
		fclose(f1);
		//free(nsm);
		//free(sim);
		//free(som);
		}

	printf("\ny is %d and dbs is %d",y,dbs);
	calculatebiterrorrate();
	}

/*------------------------------------------------THE END----------------------------------------------------*/
