# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
/****************************************************************************/
/* compile with: gcc -Wall -O3 fibwl.c ran3n.c -lm                          */ 
/*                                                                          */
/* prints to                                                                */
/*  file tmp_rt: MC sweep number, f level, energy, size of largest cluster, */
/*               length of largest cluster, width of largest cluster        */
/*  file lng_k: estimate of ln g(E) at level k (k=INIF,...,FINF-1)          */  
/*  file tmp_cluster (if k=FINF-1): MC sweep number, energy,                */
/*               no. of clusters, cluster size, (cluster shape, 4 numbers)  */
/*                                                                          */
/* reads lng_k with k=INIF-1 if INIF>0                                      */  
/****************************************************************************/
# define N 128                     /* no. of sticks                         */
# define L 51                      /* box length                            */
# define V (L*L*L) 
# define LSQ (L*L)
# define INIF 0                    /* first ln(f)=1/2^INIF                  */
# define FINF 25                   /* last ln(f)=1/2^(FINF-1)               */
# define TSW 0.65                  /* temperature in cluster move           */
# define NOBS 5                    /* no. of observables                    */
# define CUT 700                   /* -energy cutoff                        */
# define NRT 100000                /* output to rt file                     */
/************* geometry *****************************************************/
int site[N];                       /* lattice site                          */
int d1[N],d2[N];                   /* directions                            */
int occ[V];                        /* occupancy of sites                    */
int dir[6][4];                     /* neighbor list                         */
double lng[CUT];                   /* ln(DOS)                               */
/************* interactions *************************************************/
int E;                             /* energy                                */
const int e_p=5;                   /* parallel beta                         */
const int e_ap=3;                  /* antiparallel beta                     */
const int e_hh=1;                  /* HH sidechain                          */
const int e_oo=0;                  /* OO sidechain                          */
const int e_ho=0;                  /* HO sidechain                          */
/************* misc *********************************************************/
long seed=-13;                     /* random number seed                    */
/*--------------------------------------------------------------------------*/
double ran3n(long *seed);
void init(void); 
int nb(int n,int dir);   
int chk(int h[]);
int energy(int iflag);
int flip_peptide(void);
int single_peptide(void);
int cluster(int iflag);
void measurements(long mc,int nf,double o[]); 
void printheading(void); 
void runtime(long it,int nf,double o[]);

int main ()
{
  int i,nf;
  long nflp=0,nsng=0,ncls=0,mc=0;
  int h[CUT];
  double accflp=0,accsng=0,acccls=0;
  double lnf=1,o[NOBS],tmp;
  double frflp=0.333,frsng=0.333;
  char filename[12];
  FILE *fp;

  printheading();
  init();
  for (i=0;i<INIF;i++) lnf/=2;
  measurements(-1,0,o);

  for (nf=INIF;nf<FINF;nf++) { 
    printf("beginning on f level %i, ln f=%10.8f\n",nf,lnf); fflush(0);
    for (i=0;i<CUT;i++) h[i]=0;
    while (chk(h)==0) { 
      for (i=0;i<N;i++) {
        if ((tmp=ran3n(&seed))<frflp) {
	  nflp++; accflp+=flip_peptide(); 
	}
	else if (tmp<frflp+frsng) {
	  nsng++; accsng+=single_peptide();
	}
	else {
 	  ncls++; acccls+=cluster(0);
 	}    
	h[-E]++;
        lng[-E]+=lnf;
      }   
      mc++;
      if (mc%NRT==0) {
	measurements(mc,nf,o);
	runtime(mc,nf,o);
      }
    }
    lnf/=2;
    tmp=lng[0]; for (i=0;i<CUT;i++) lng[i]-=tmp;
    sprintf(filename,"lng_%i",nf); 
    fp=fopen(filename,"w");
    for (i=0;i<CUT;i++) fprintf(fp,"%i %f\n",-i,lng[i]);
    fclose(fp);
  }

  printf("--- Acceptance ---\n"); 
  if (nflp>0) printf("nflp %10ld  accflp %5.3f\n",nflp,accflp/nflp);
  if (nsng>0) printf("nsng %10ld  accsng %5.3f\n",nsng,accsng/nsng);
  if (ncls>0) printf("ncls %10ld  acccls %5.3f\n",ncls,acccls/ncls);
  
  return 0;
}
/****************************************************************************/
int chk(int h[])
/* returns 1 if f is to be updated, 0 otherwise                             */
{
  double have=0,tmp;
  const int delta=3;
  int i;
  for (i=0;i<CUT;i++) have+=h[i];
  if (have<0.5) have=1; 
  have*=.7/CUT;
  tmp=0; for (i=0;i<2*delta;i++) tmp+=h[i];
  if (tmp/(2*delta)<have) return 0;
  for (i=delta;i<CUT-delta;i++) {
    tmp+=h[i+delta]-h[i-delta];
    if (tmp/(2*delta)<have) return 0;
  }
  return 1;
}   
/****************************************************************************/
/***** ENERGY ***************************************************************/
/****************************************************************************/
int energy(int iflag)
{
  extern int mepair(int i,int j,int rij);
  int e=0;
  int i,j,k,rij;

  if (iflag<0 || iflag==N) {
    if (iflag<0) {mepair(-1,0,0); printf("initialized energy\n");}
    for (i=0;i<N;i++) {
      for (k=0;k<4;k++) {
        rij=dir[d1[i]][k]; 
	if ((j=occ[nb(i,rij)])<0) continue;
        if (i>j) continue;
	e+=mepair(i,j,rij);
      }
    }  
  }
  else if (iflag>=0 && iflag<N) {
    for (k=0;k<4;k++) {
      rij=dir[d1[iflag]][k];
      if ((j=occ[nb(iflag,rij)])<0) continue;
      e+=mepair(iflag,j,rij);
    }    
  }
  else 
    printf("energy: error\n");

  return -e;
}    
/****************************************************************************/
int mepair(int i,int j,int rij)
{
  static int xp[6][6];
  static int r_p=1+e_p,r_ap=1+e_ap,r_hh=1+e_hh,r_ho=1+e_ho,r_oo=1+e_oo;  

  if (i<0) {
    xp[0][1]=2; xp[1][0]=5;
    xp[0][4]=5; xp[4][0]=2;
    xp[3][1]=5; xp[1][3]=2;
    xp[3][4]=2; xp[4][3]=5;
    xp[0][2]=4; xp[2][0]=1;
    xp[0][5]=1; xp[5][0]=4;
    xp[3][2]=1; xp[2][3]=4;
    xp[3][5]=4; xp[5][3]=1;
    xp[1][2]=0; xp[2][1]=3;
    xp[1][5]=3; xp[5][1]=0;
    xp[4][2]=3; xp[2][4]=0;
    xp[4][5]=0; xp[5][4]=3;
    printf("initialized mepair\n");
    return -1;
  }

  if (((d1[i]-d1[j]+3)%3)!=0) return 0;      
  if (((d2[i]-d2[j]+3)%3)!=0) return 2;
  if (((rij-d2[i]+3)%3)==0) {  
    if ((d1[i]==d1[j]) && (d2[i]==d2[j]))                  /* parallel beta */
      return r_p;                          
    else if ((d1[i]!=d1[j]) && (d2[i]!=d2[j]))         /* antiparallel beta */
      return r_ap;
  }
  else {
    if ((xp[d1[i]][d2[i]]==rij) && (xp[d1[j]][d2[j]]!=rij))           /* HH */
      return r_hh;
    else if ((xp[d1[i]][d2[i]]!=rij) && (xp[d1[j]][d2[j]]==rij))      /* OO */
      return r_oo;
    else                                                              /* HO */
      return r_ho; 
  }
  return 1;
}
/****************************************************************************/
/***** UPDATES **************************************************************/
/****************************************************************************/
int single_peptide(void)
{
  int i,k,d1o,d2o,de;
  int siteo,siten;

  i=N*ran3n(&seed);
  siteo=site[i];
  de=-energy(i); 
  k=6*ran3n(&seed);  
  if (occ[(siten=nb(i,k))]>=0) return 0;
  site[i]=siten;
  occ[siteo]=-1;
  occ[site[i]]=i;
  d1o=d1[i]; d2o=d2[i];
  d1[i]=6*ran3n(&seed); 
  while (((d2[i]=6*ran3n(&seed))-d1[i]+3)%3==0) {};
  de+=energy(i);
  if (E+de>-CUT && ran3n(&seed)<exp(lng[-E]-lng[-(E+de)])) {
    E+=de;
    return 1;
  }
  else {
    occ[siteo]=i;
    occ[site[i]]=-1;
    site[i]=siteo;
    d1[i]=d1o; d2[i]=d2o;
    return 0;
  }
}
/****************************************************************************/
int flip_peptide(void)
{
  int i,d1o,d2o,de;

  i=N*ran3n(&seed);
  de=-energy(i); 
  d1o=d1[i]; d2o=d2[i];
  d1[i]=6*ran3n(&seed); 
  while (((d2[i]=6*ran3n(&seed))-d1[i]+3)%3==0) {};
  de+=energy(i);  
  if (E+de>-CUT && ran3n(&seed)<exp(lng[-E]-lng[-(E+de)])) {
    E+=de;
    return 1;
  }
  else {
    d1[i]=d1o; d2[i]=d2o;
    return 0;
  }
}
/****************************************************************************/
void add_ngb(int ppt,int list[])
{
  extern int mepair(int i,int j,int rij); 
  int k,dn,pptn;
  static double prb[13];

  if (ppt<0) {
    for (k=0;k<13;k++) prb[k]=1-exp((-k)/TSW);
    return;
  } 

  for (k=0;k<4;k++) {
    dn=dir[d1[ppt]][k];
    if (((pptn=occ[nb(ppt,dn)])>=0) && (list[pptn]==1) && 
	(ran3n(&seed)<prb[mepair(ppt,pptn,dn)])) {
      list[pptn]=-1;
      add_ngb(pptn,list);
    } 
  }
}
/****************************************************************************/
int cluster(int iflag)
{
  extern void add_ngb(int n,int list[]);
  extern int rotd(int xyz,int d,int sgn);
  int Enew;
  int j,n,rij,tmp,ra;
  int nc=0;
  int list[N],move[N],siten[N],siteo[N];
  int d1o[N],d2o[N];
  int x0,y0,z0,xt,yt,zt,sgn;
  double type;
  static double bsw;

  if (iflag<0) {
      add_ngb(-1,list);
      bsw=1/TSW;
      printf("initialized cluster\n");
      return -1;
  }

  for (j=0;j<N;j++) list[j]=1;
  n=N*ran3n(&seed);
  list[n]=-1;
  add_ngb(n,list);
  for (j=0;j<N;j++) {
    if (list[j]==-1) {move[nc++]=j;}
  }

  if (nc==1) return 0;

  type=ran3n(&seed);
  if (type<.5) {                                     /***** translation *****/
    rij=6*ran3n(&seed);
    for (j=0;j<nc;j++) {
      if (((tmp=occ[(siten[j]=nb(move[j],rij))])>=0) && (list[tmp]==1)) 
	return 0;
    }
    for (j=0;j<nc;j++) {
      siteo[j]=site[move[j]]; 
      occ[siteo[j]]=-1;
    }
    for (j=0;j<nc;j++) {
      occ[siten[j]]=move[j];
      site[move[j]]=siten[j];
    } 
    Enew=energy(N);
    if (Enew>-CUT && ran3n(&seed)<exp(lng[-E]-lng[-Enew]+bsw*(Enew-E))) {
      E=Enew;
      return 1;
    }
    else {
      for (j=0;j<nc;j++) {
        site[move[j]]=siteo[j]; 
        occ[siten[j]]=-1;
      } 
      for (j=0;j<nc;j++) occ[siteo[j]]=move[j]; 
      return 0;
    }
  }   
  else {                                               /***** rotations *****/ 
    for (j=0;j<nc;j++) {d1o[j]=d1[move[j]]; d2o[j]=d2[move[j]];} 
    x0=site[n]%L; y0=(site[n]/L)%L; z0=site[n]/LSQ;
    sgn=2*ran3n(&seed); sgn=2*sgn-1;
    if (type<4./6.) {          
      ra=0;
      for (j=0;j<nc;j++) {
        yt=y0+sgn*((site[move[j]]/LSQ)-z0); yt=(yt+L)%L;
        zt=z0-sgn*(((site[move[j]]/L)%L)-y0); zt=(zt+L)%L;
        siten[j]=(site[move[j]]%L)+yt*L+zt*LSQ;
      }
    }
    else if (type<5./6.) {
      ra=1;
      for (j=0;j<nc;j++) {
        xt=x0+sgn*((site[move[j]]/LSQ)-z0); xt=(xt+L)%L; 
        zt=z0-sgn*((site[move[j]]%L)-x0); zt=(zt+L)%L;
        siten[j]=xt+((site[move[j]]/L)%L)*L+zt*LSQ;
      }
    }
    else {
      ra=2;
      for (j=0;j<nc;j++) {
        xt=x0+sgn*(((site[move[j]]/L)%L)-y0); xt=(xt+L)%L;
        yt=y0-sgn*((site[move[j]]%L)-x0); yt=(yt+L)%L;
        siten[j]=xt+yt*L+(site[move[j]]/LSQ)*LSQ;
      }
    }
    for (j=0;j<nc;j++) {
      if ((tmp=occ[siten[j]])>=0 && list[tmp]==1) return 0;
    }
    for (j=0;j<nc;j++) {
      siteo[j]=site[move[j]];
      occ[siteo[j]]=-1;
    } 
    for (j=0;j<nc;j++) {
      occ[siten[j]]=move[j];
      site[move[j]]=siten[j];
      d1[move[j]]=rotd(ra,d1[move[j]],sgn); 
      d2[move[j]]=rotd(ra,d2[move[j]],sgn); 
    }
    Enew=energy(N);
    if (Enew>-CUT && ran3n(&seed)<exp(lng[-E]-lng[-Enew]+bsw*(Enew-E))) {
      E=Enew;
      return 1;
    }
    else {
      for (j=0;j<nc;j++) occ[siten[j]]=-1; 
      for (j=0;j<nc;j++) {
        occ[siteo[j]]=move[j]; 
        site[move[j]]=siteo[j]; 
        d1[move[j]]=d1o[j]; 
	d2[move[j]]=d2o[j]; 
      }
      return 0;
    }
  }
}
/****************************************************************************/
int rotd(int xyz,int d,int sgn)
{
  int rd;
  if (xyz==0) {
    switch (d) {
    case 0: rd=0; break;
    case 1: rd=(sgn==1)?5:2;  break;
    case 2: rd=(sgn==1)?1:4;  break;
    case 3: rd=3; break;
    case 4: rd=(sgn==1)?2:5;  break; 
    case 5: rd=(sgn==1)?4:1;  break; 
    default: printf("error x rotd\n"); exit(-17); break;        
    }
  }
  else if (xyz==1) {
    switch (d) {
    case 0: rd=(sgn==1)?5:2;  break;
    case 1: rd=1; break;
    case 2: rd=(sgn==1)?0:3;  break;
    case 3: rd=(sgn==1)?2:5;  break;
    case 4: rd=4; break; 
    case 5: rd=(sgn==1)?3:0;  break; 
    default: printf("error y rotd\n"); exit(-17); break;
    }
  }
  else if (xyz==2) {
    switch (d) {
    case 0: rd=(sgn==1)?4:1;  break;
    case 1: rd=(sgn==1)?0:3;  break;
    case 2: rd=2; break;
    case 3: rd=(sgn==1)?1:4;  break; 
    case 4: rd=(sgn==1)?3:0;  break;
    case 5: rd=5; break; 
    default: printf("error z rotd\n"); exit(-17); break;
    }
  }
  else {
    printf("error in rotd\n"); exit(-17); 
  }
  return rd;
}
/****************************************************************************/
/***** MEASUREMENTS *********************************************************/
/****************************************************************************/
void measurements(long mc,int nf,double o[])
{
  int mxs;
  double dim1,dim2;
  extern void cluster_analysis(long mc,int nf,int *mxs,double *dim1,double *dim2);
  if (mc<0) cluster_analysis(-1,nf,&mxs,&dim1,&dim2); 
  o[0]=E; 
  o[1]=E*E;
  cluster_analysis(mc,nf,&mxs,&dim1,&dim2);
  o[2]=mxs;
  o[3]=dim1;
  o[4]=dim2;
  return;
}
/****************************************************************************/
void cluster_analysis(long mc,int nf,int *mxs,double *dim1,double *dim2)
{
  extern void inertia_egn(int s,int list[],double *dim1, double *dim2);
  int i,j,k,m,s,ic,d1c;
  int done[N],list[N];
  int ncd[2];
  short cd[4];
  int cnt1l[L];
  double cb,da,db;
  static double da1[L];
  FILE *fp;

  if (mc<0) {
    da1[1]=0;
    for (i=1;i<L/2;i++) da1[2*i+1]+=da1[2*i-1]+2*i*i;  
    da1[0]=0;
    for (i=1;i<L/2;i++) da1[2*i]+=da1[2*(i-1)]+2*(i-.5)*(i-.5);
    for (i=1;i<L;i++) da1[i]=2*sqrt(da1[i]/i);  
    return; 
  }
  fp=fopen("tmp_cluster","a");
  (*mxs)=0;
  for (i=0;i<N;i++) {done[i]=0;}
  for (i=0;i<L;i++) {cnt1l[i]=0;}

  for (i=0;i<N;i++) {
    if (done[i]==1) continue;
    done[i]=1; 
    list[(ic=0)]=i;
    s=1;
    d1c=d1[i];
    ncd[0]=ncd[1]=0;
    for (k=0;k<4;k++) {
      if (dir[d1c][k]%3==0) {cd[k]=0;}
      else if (dir[d1c][k]%3==2) {cd[k]=1;}
      else {cd[k]=(d1c%3==0) ? 0:1;}
    }
    while (ic<s) {
      m=list[ic];
      for (k=0;k<4;k++) {
        if ((j=occ[nb(m,dir[d1c][k])])<0) continue;
	if ((abs((d1c-d1[j])%3)!=0)) continue;
        ncd[cd[k]]++; 
        if (done[j]==0) {list[s++]=j; done[j]=1;}
      }
      ic++; 
    }
    if (ncd[1]>ncd[0]) {
      j=ncd[1];
      ncd[1]=ncd[0];
      ncd[0]=j;
    }
    ncd[0]/=2; cb=(ncd[1]/=2);
    cb/=s;
    cb=1./(1-cb);
    if (ncd[1]>0) {
      inertia_egn(s,list,&da,&db);
      if (nf==FINF-1) 
        fprintf(fp,"%li %i %i %i %f %f %f %f\n",mc,E,1,s,s/cb,cb,da,db);
    }
    else {
      if (ncd[0]>=L) {printf("too long cluster\n");}
      cnt1l[ncd[0]]++; 
    }    
    if (s>(*mxs)) {
      *mxs=s; 
      if (ncd[1]>0) {
	*dim1=sqrt(3*da*da+1); *dim2=sqrt(3*db*db+1);
      }
      else {
        *dim1=s; *dim2=1;
      }
    }
  }

  for (i=L-1;i>=0;i--) {
    if (nf==FINF-1 && cnt1l[i]>0) fprintf(fp,"%li %i %i %i %f %f %f %f\n",
		    mc,E,cnt1l[i],i+1,i+1.0,1.0,da1[i+1],0.0);
  }
  fclose(fp);
}
/****************************************************************************/
void inertia_egn(int s,int list[],double *da, double *db)
{
  int i,j,d;
  int a[N],b[N],o[L];
  int min0;
  double sa=0,sb=0;
  double ax,bx;
  double iaa=0,ibb=0,iab=0;
  double s3=(1.0/s)/(s*s);

  d=d1[list[0]]%3;  
  for (i=0;i<s;i++) {
    j=site[list[i]];
    if (d==0) { 
      a[i]=(j/L)%L; b[i]=j/LSQ;
    } 
    else if (d==1) {
      a[i]=j%L; b[i]=j/LSQ;
    }
    else {
      a[i]=j%L; b[i]=(j/L)%L;
    }
  }

  for (i=0;i<L;i++) o[i]=0;
  for (i=0;i<s;i++) o[a[i]]=1;
  if (o[0]==1 && o[L-1]==1) {
    min0=0;
    while (o[min0]==1) min0++;
    for (i=0;i<s;i++) {
      if (a[i]<min0) a[i]+=L;
    }
  } 
  for (i=0;i<L;i++) o[i]=0;
  for (i=0;i<s;i++) o[b[i]]=1;
  if (o[0]==1 && o[L-1]==1) {
    min0=0;
    while (o[min0]==1) min0++;
    for (i=0;i<s;i++) {
      if (b[i]<min0) b[i]+=L;
    }
  } 

  for (i=0;i<s;i++) { 
    sa+=a[i]; sb+=b[i]; 
  }
  for (i=0;i<s;i++) { 
    ax=s*a[i]-sa; bx=s*b[i]-sb; 
    iaa+=bx*bx;
    ibb+=ax*ax;
    iab+=ax*bx;
  }
  iaa*=s3; ibb*=s3; iab*=s3;  
  *da=.5*(iaa+ibb)+sqrt(iab*iab+.25*(iaa-ibb)*(iaa-ibb));
  *db=iaa+ibb-(*da);
  if (*db<0) *db=0;
  *da=2*sqrt(*da);
  *db=2*sqrt(*db);
  return;
}  
/****************************************************************************/
/***** INPUT/OUTPUT *********************************************************/
/****************************************************************************/
void printheading(void)
{
  printf("***********************************\n");
  printf("*          fibwl.c                *\n");
  printf("***********************************\n");
  printf("Number of peptides N = %i\n",N);
  printf("Box size L = %i\n",L); 
  printf("Energy parameters e_p = %i   e_ap = %i   e_hh = %i   e_oo = %i   e_ho= %i\n",e_p,e_ap,e_hh,e_oo,e_ho);
  printf("Lower energy cut %i\n",-CUT); 
  printf("MC parameters:");
  printf("  NRT %i  INIF %i  FINF %i  TSW %f seed %li\n",
         NRT,INIF,FINF,TSW,seed);
  printf("-----\n\n");
  fflush(0);
} 
/****************************************************************************/
void runtime(long mc,int nf,double o[]) 
{
  int i;
  FILE *fp;
  fp=fopen("tmp_rt","a");
  fprintf(fp,"%li %i ",mc,nf);
  for (i=0;i<NOBS;i++) {if (i!=1) fprintf(fp,"%f ",o[i]);}
  fprintf(fp,"\n");
  fclose(fp);
}
/****************************************************************************/
/***** INITIALIZATION *******************************************************/
/****************************************************************************/
void init(void)
{
  int i,j;
  char filename[12];
  FILE *fp;

  if (INIF==0) {
    for (i=0;i<CUT;i++) lng[i]=0;
  }
  else {
    sprintf(filename,"lng_%i",INIF-1);
    fp=fopen(filename,"r");
    for (i=0;i<CUT;i++) fscanf(fp,"%i %lf\n",&j,&lng[i]);
    fclose(fp);
  }

  for (i=0;i<V;i++) occ[i]=-1;
  for (i=0;i<N;i++) {
    while (occ[(site[i]=V*ran3n(&seed))]>=0) {}; 
    occ[site[i]]=i; 
    d1[i]=6*ran3n(&seed);
    while (abs((d2[i]=6*ran3n(&seed))-d1[i])%3==0) {};
  }

  /* dir[i][k] lists directions perpendicular to i, k=0,1,2,3 */
  dir[0][0]=dir[3][0]=1;
  dir[0][1]=dir[3][1]=2;
  dir[0][2]=dir[3][2]=4;
  dir[0][3]=dir[3][3]=5;
  dir[1][0]=dir[4][0]=0;
  dir[1][1]=dir[4][1]=2;
  dir[1][2]=dir[4][2]=3;
  dir[1][3]=dir[4][3]=5;
  dir[2][0]=dir[5][0]=0;
  dir[2][1]=dir[5][1]=1;
  dir[2][2]=dir[5][2]=3;
  dir[2][3]=dir[5][3]=4;

  cluster(-1);
  E=energy(-1);
}
/****************************************************************************/
/***** NAVIGATOR ************************************************************/
/****************************************************************************/
# define A (L-1)
# define B (L-1)*L
# define C (L-1)*LSQ
int nb(int n,int dir)
{
  int delta;
  switch (dir) {
  case 0: delta=(site[n]%L==A) ? -A:1; break;
  case 1: delta=((site[n]/L)%L==A) ? -B:L; break;  
  case 2: delta=(site[n]/LSQ==A) ? -C:LSQ; break;  
  case 3: delta=(site[n]%L==0) ? A:-1; break;  
  case 4: delta=((site[n]/L)%L==0) ? B:-L; break;  
  case 5: delta=(site[n]/LSQ==0) ? C:-LSQ; break;
  default: printf("error in nb\n"); exit(-17);
  }
  return site[n]+delta;  
}
# undef A
# undef B
# undef C
/****************************************************************************/

