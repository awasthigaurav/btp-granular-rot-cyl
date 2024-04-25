#include <stdio.h>
#include <math.h>
long i,j,N,M,Nparts,id,iw,Nr;
long ia,im,iq,ir,idum;
double am,rho,Ly,eps,r,dth;
double pi,d,th,sp,di,deld;
double x[9000],z[9000],R;
FILE *inp,*out;

double ran0(void)
{
  long k;
  double res;
  k=idum/iq;
  idum=ia*(idum-k*iq)-ir*k;
  if (idum < 0) idum=idum+im;
  res=2*am*idum-1;
  return res;
}

int main( )
{
  ia = 16807;
  im = 2147483647;
  am = 1.0/im;
  iq = 127773;
  ir = 2836;

  pi=atan(1.0)*4;
  d=0.1;
  deld=0.1*d;
  rho=2.5;
  R=5;
  N=500;
  M=20;
  Ly=M*d/2;
  eps=1e-6;
  idum=31425617;
  
  for (i=1;i<=N;++i){
    th=2*pi/N*i;
    x[i]=(R-d/2)*cos(th);
    z[i]=(R-d/2)*sin(th);
  }
  Nparts=N*M;
  out = fopen("cyl.data","w");
  fprintf(out,"LAMMPS data file via write_data, version 12 Dec 2018, timestep = %i\n",0);
  fprintf(out,"\n");
  fprintf(out,"%li atoms\n",Nparts);
  fprintf(out,"%i atom types\n",2);
  fprintf(out,"\n");
  fprintf(out,"%lf %lf xlo xhi\n",-R,R);
  fprintf(out,"%lf %lf ylo yhi\n",0.0,Ly);
  fprintf(out,"%lf %lf zlo zhi\n",-R,R);
  fprintf(out,"\n");
  fprintf(out,"Atoms # sphere\n");
  fprintf(out,"\n");
  id=0;
  for (j=1;j<=M;++j){
    for (i=1;i<=N;++i){
      id=id+1;
      di=d+deld*ran0();
      fprintf(out,"%li %i %f %f %f %f %f %i %i %i\n", id,1,di,rho,x[i],d/2*j-d/4,z[i],0,0,0);
    }
  }
  fprintf(out,"\n");
  fprintf(out,"Velocities\n");
  fprintf(out,"\n");
  for (i=1;i<=Nparts;++i)
    fprintf(out,"%li %f %f %f %f %f %f\n",i,0.0,0.0,0.0,0.0,0.0,0.0);
  fclose(out);
  
  sp=1.1*d;N=5000;
  r=(R-d/2)-2*sp;dth=sp/r;th=0;
  for (i=1;i<=N;++i){
    th=th+dth;
    if (th>2*pi-2*dth){
      r=r-sp;
      dth=sp/r;
      th=0;
    }
    x[i]=r*cos(th);z[i]=r*sin(th);
  }
  M=trunc(Ly/sp);
  Nparts=N*(M-1);
  out = fopen("pos.data","w");
  fprintf(out,"LAMMPS data file via write_data, version 12 Dec 2018, timestep = %i\n",0);
  fprintf(out,"\n");
  fprintf(out,"%li atoms\n",Nparts);
  fprintf(out,"%i atom types\n",3);
  fprintf(out,"\n");
  fprintf(out,"%lf %lf xlo xhi\n",-R,R);
  fprintf(out,"%lf %lf ylo yhi\n",0.0,Ly);
  fprintf(out,"%lf %lf zlo zhi\n",-R,R);
  fprintf(out,"\n");
  fprintf(out,"Atoms # sphere\n");
  fprintf(out,"\n");
  id=0;
  for (j=1;j<=M-1;++j){
    for (i=1;i<=N;++i){
      id=id+1;
      di=d+deld*ran0();
      fprintf(out,"%li %i %f %f %f %f %f %i %i %i\n", id,2,di,rho,x[i],sp*(j-0.5),z[i],0,0,0);
    }
  }
  fprintf(out,"\n");
  fprintf(out,"Velocities\n");
  fprintf(out,"\n");
  for (i=1;i<=Nparts;++i)
    fprintf(out,"%li %f %f %f %f %f %f\n",i,0.0,0.0,0.0,0.0,0.0,0.0);
  fclose(out);

}


