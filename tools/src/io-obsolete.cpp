#include <config.h>

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.h"
#include "constants.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int section_loadr1(char *name, int nndes, float *buffer,int column)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,n,nitems;
  float dum;
  FILE *in;
  char c;

  in=fopen(name,"r");
  if(in==0)
    TRAP_ERR_RETURN(errno,1,"fopen(\"%s\",\"r\") error (%d %s)\n",name,errno,strerror(errno));

  /* first 3 lines should be comment only */
  for(n=0;n<nndes;n++) {
    for (k=0;k<2+column; k++) nitems=fscanf(in,"%f",&dum);
    nitems=fscanf(in,"%f",&buffer[n]);
    do { c=fgetc(in);   }  while ((c != '\n') && !feof(in));
    }

  fclose(in);
  return (0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_qoddy_inq(const string &path,poc_global_t *global,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  size_t nndes;
  const size_t n=1000;
  char  line[n], *s;
  FILE *in;
  
  in=fopen(path.c_str(),"r");
  if(in==0)
    TRAP_ERR_RETURN(errno,1,"fopen(\""+path+"\",\"r\") error (%d %s)\n",errno,strerror(errno));
  
  const string extension=get_extension(path);
  
  if( extension.size()!=3 or
      strchr("sv",extension[0])==0 or
      extension[1]!='2' or
      strchr("rc",extension[2])==0 ){
    if(verbose>=0) STDERR_BASE_LINE_FUNC("Extension "+extension+" of file "+path+" is not of quoddy format\n");
    return -1;
    }
  
/*------------------------------------------------------------------------------
  first 2 or 3 lines should be comments line only */
  s=fgets(line,n,in);
  s=fgets(line,n,in);
  if(extension[2]=='c')
    s=fgets(line,n,in);
  
/*------------------------------------------------------------------------------
  count nodes */
  nndes=0;
  while(s!=0){
    s=fgets(line,n,in);
    if(s==0)
      break;
    nndes++;
    }
  
  fclose(in);
  
  if(nndes==0){
    TRAP_ERR_RETURN(NC_EDIMSIZE,1,"Not enough lines in \""+path+"\"\n");
    }
  
/*------------------------------------------------------------------------------
  fill header */
  const poc_dim_t N("N",nndes);
  poc_var_t var("",NC_FLOAT);
  var.dimensions<<N;
  
  global->dimensions.clear(1);
  global->attributes.clear();
  
  if(extension=="s2r"){
    global->variables.clear(1);
    var.name="H";
    global->variables<<var;
    }
  else if(extension=="s2c"){
    global->variables.clear(2);
    var.name="Ha";
    global->variables<<var;
    var.name="Hg";
    global->variables<<var;
    }
  else if(extension=="v2r"){
    global->variables.clear(2);
    var.name="U";
    global->variables<<var;
    var.name="V";
    global->variables<<var;
    }
  else if(extension=="v2c"){
    global->variables.clear(4);
    var.name="Ua";
    global->variables<<var;
    var.name="Ug";
    global->variables<<var;
    var.name="Va";
    global->variables<<var;
    var.name="Vg";
    global->variables<<var;
    }
  else{
    TRAP_ERR_EXIT(ENOEXEC,"%s programming error with extension "+extension+"\n",__func__);
    }
  
  return (0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_loadr1(const char *name, int nndes, short *buffer, char *text)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,start=0,dum,nitems;
  char  *line;
  FILE *in;
  
  line=new char[1000];

  in=fopen(name,"r");
  if(in==0)
    TRAP_ERR_RETURN(errno,1,"fopen(\"%s\",\"r\") error (%d %s)\n",name,errno,strerror(errno));

  /* first 3 lines should be comment only */
  n=999;
  fgets(line,n,in);
  fgets(line,n,in);
  
  if(text!=0) strcpy(text,line);

  nitems=sscanf(line,"%d %f",&dum,&buffer[0]);

//  if((nitems ==2) && (dum==1)) start=1;
  while((nitems !=2) || (dum!=1)) {
    fgets(line,n,in);
    nitems=sscanf(line,"%d %d",&dum,&buffer[0]);
    if(feof(in)) goto error;
    }

  start=1;

  for(n=start;n<nndes;n++) {
    nitems=fscanf(in,"%d %d",&dum,&buffer[n]);
    if(nitems !=2) goto error;
    }

  delete[] line;
  fclose(in);
  return (0);

 error:
  delete[] line;
  fclose(in);
  return (-1);

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_loadr1(const char *name, int nndes, float *buffer, char *text)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,start=0,dum,nitems;
  char  *line;
  FILE *in;
  
  line=new char[1000];

  in=fopen(name,"r");
  if(in==0)
    TRAP_ERR_RETURN(errno,1,"fopen(\"%s\",\"r\") error (%d %s)\n",name,errno,strerror(errno));

  /* first 3 lines should be comment only */
  n=999;
  fgets(line,n,in);
  fgets(line,n,in);
  
  if(text!=0) strcpy(text,line);

  nitems=sscanf(line,"%d %f",&dum,&buffer[0]);

//  if((nitems ==2) && (dum==1)) start=1;
  while((nitems !=2) || (dum!=1)) {
    fgets(line,n,in);
    nitems=sscanf(line,"%d %f",&dum,&buffer[0]);
    if(feof(in)) goto error;
    }

  start=1;

  for(n=start;n<nndes;n++) {
    nitems=fscanf(in,"%d %f",&dum,&buffer[n]);
    if(nitems !=2) goto error;
    }

  delete[] line;
  fclose(in);
  return (0);

 error:
  delete[] line;
  fclose(in);
  return (-1);

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_loadr1(const char *name, int nndes, double *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,start=0,dum,nitems,status=0;
  char  line[1000];
  FILE *in;

  in=fopen(name,"r");
  if(in==0)
    TRAP_ERR_RETURN(errno,1,"fopen(\"%s\",\"r\") error (%d %s)\n",name,errno,strerror(errno));

  /* first 3 lines should be comment only */
  n=999;
  fgets(line,n,in);
  fgets(line,n,in);

  nitems=sscanf(line,"%d %lf",&dum,&buffer[0]);

  if((nitems ==2) && (dum==1)) start=1;

  for(n=start;n<nndes;n++) {
    nitems=fscanf(in,"%d %lf",&dum,&buffer[n]);
    if(nitems !=2){
      status=errno;
      break;
      }
    }

  fclose(in);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_loadr1(const char *name, int nndes, short *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=quoddy_loadr1(name, nndes, buffer, (char *) 0);
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_loadr1(const char *name, int nndes, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=quoddy_loadr1(name, nndes, buffer, (char *) 0);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_loadd1(const char *name, int nndes, double *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,start=0,dum,nitems;
  char  line[1000];
  FILE *in;

  in=fopen(name,"r");
  if(in==0)
    TRAP_ERR_RETURN(errno,1,"fopen(\"%s\",\"r\") error (%d %s)\n",name,errno,strerror(errno));

  /* first 3 lines should be comment only */
  n=999;
  fgets(line,n,in);
  fgets(line,n,in);

  nitems=sscanf(line,"%d %lf",&dum,&buffer[0]);

  if((nitems ==2) && (dum==1)) start=1;

  for(n=start;n<nndes;n++) {
    nitems=fscanf(in,"%d %lf",&dum,&buffer[n]);
    if(nitems !=2) goto error;
    }

  fclose(in);
  return (0);

 error:
  fclose(in);
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_loadr2(const char *name, int nndes, float *bufx, float *bufy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,start=0,dum,nitems;
  char  line[1000];
  FILE *in;

  in=fopen(name,"r");
  if(in == NULL) {
    printf("Unable to open %s\n",name);
    return (-1);
    }

  /* first 3 lines should be comment only */
  n=999;
  fgets(line,n,in);
  fgets(line,n,in);

  n=0;
  nitems=sscanf(line,"%d %lf",&dum,&bufx[n],&bufy[n]);
  if((nitems ==2) && (dum==1)) start=1;

  for(n=start;n<nndes;n++) {
    fscanf(in,"%d %f %f",&dum,&bufx[n],&bufy[n]);
    }

  fclose(in);
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_loadr2(const char *name, int nndes, double *bufx, double *bufy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,start=0,dum,nitems;
  char  line[1000];
  FILE *in;

  in=fopen(name,"r");
  if(in == NULL) {
    printf("Unable to open %s\n",name);
    return (-1);
    }

  /* first 3 lines should be comment only */
  n=999;
  fgets(line,n,in);
  fgets(line,n,in);

  n=0;
  nitems=sscanf(line,"%d %lf %lf",&dum,&bufx[n],&bufy[n]);
  if((nitems ==3) && (dum==1)) start=1;

  for(n=start;n<nndes;n++) {
    fscanf(in,"%d %lf %lf",&dum,&bufx[n],&bufy[n]);
    }

  fclose(in);
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int quoddy_loadc1_template(const char *name, int nndes, complex<T> *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, nitems;
  
  const size_t n=1000;
  char  line[n], *s;
  double x,y;
  double a,p;
  FILE *in;

  in=fopen(name,"r");
  if(in==0)
    TRAP_ERR_RETURN(errno,1,"fopen(\"%s\",\"r\") error (%d %s)\n",name,errno,strerror(errno));

/*------------------------------------------------------------------------------
  first 3 lines should be comments line only */
  fgets(line,n,in);
  fgets(line,n,in);
  fgets(line,n,in);
  
  for(i=0;i<nndes;i++) {
    s=fgets(line,n,in);
    if(s==0) break;
    nitems=count_token(line);
    if(nitems==2) {
      nitems=sscanf(line, "%lf %lf",&a,&p);
      }
    else if(nitems==3) {
      nitems=sscanf(line,"%*d %lf %lf",&a,&p);
      }
    else if(nitems==4) {
      nitems=sscanf(line,"%lf %lf %lf %lf",&x,&y,&a,&p);
      }
    else{
      fclose(in);
      TRAP_ERR_RETURN(-1,"%s: error reading %s : %d items found on line %d",__func__,name,nitems,i+3);
      }
//     fscanf(in,"%*d %f %f",&a,&p);
    buffer[i]=polar<T>(a,-p*d2r);
    }

  fclose(in);
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_loadc1(const char *name, int nndes, complex<float> *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=quoddy_loadc1_template(name,nndes,buffer);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_loadc1(const char *name, int nndes, complex<double> *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=quoddy_loadc1_template(name,nndes,buffer);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_loadc2(char *name, int nndes, fcomplex *bufferx, fcomplex *buffery)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,nitems,dum,count=0;
  float a,p;
  FILE *in;
  char c;

  in=fopen(name,"r");
  if(in==0)
    TRAP_ERR_RETURN(errno,1,"fopen(\"%s\",\"r\") error (%d %s)\n",name,errno,strerror(errno));
 /*  warning, should have 2 or 3 lines in header (3=Quoddy) */
  do { c=fgetc(in);   }  while ((c != '\n') && !feof(in));
  do { c=fgetc(in);   }  while ((c != '\n') && !feof(in));

  for(n=0;n<nndes;n++) {
    redo:
    nitems=fscanf(in,"%d %f %f",&dum,&a,&p);
    if(nitems !=3) {
      do { c=fgetc(in);   }  while ((c != '\n') && !feof(in));
      goto redo;
      }
/*
    bufferx[n].r=a*cos(-p*d2r);
    bufferx[n].i=a*sin(-p*d2r);
*/
    bufferx[n]=fcomplex(a*cos(-p*d2r),a*sin(-p*d2r));
    nitems=fscanf(in,"%f %f",&a,&p);
    if(nitems !=2) goto error;
/*
    buffery[n].r=a*cos(-p*d2r);
    buffery[n].i=a*sin(-p*d2r);
*/
    buffery[n]=fcomplex(a*cos(-p*d2r),a*sin(-p*d2r));
    count++;
    }

  if(count!=nndes) goto error;

  fclose(in);
  return (0);

 error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int load_s2r(char *filename, float *buf, int nnde,  int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int   n,dum;
  char  line[1000];
  FILE *in;

  in=fopen(filename,"r");
  if(in == NULL) {
    printf("Unable to open %s\n",filename);
    return (-1);
    }

  n=999;
  fgets(line,n,in);
  fgets(line,n,in);
  for(n=1;n<=nnde;n++) {
    fscanf(in,"%d %f",&dum,&buf[n-1]);
    }

  fclose(in);
  return (0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int load_v2r(char *filename, float *bufx, float *bufy, int nnde,  int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int   n,dum;
  char  line[1000];
  FILE *in;

  in=fopen(filename,"r");

  if(in == NULL) {
    printf("Unable to open %s\n",filename);
    return (-1);
    }

  n=999;
  fgets(line,n,in);
  fgets(line,n,in);
  for(n=1;n<=nnde;n++) {
    fscanf(in,"%d %f %f",&dum,&bufx[n-1],&bufy[n-1]);
    }

  fclose(in);
  return (0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_saver1(const char *name, int nndes, float *buffer, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,nitems;
  FILE *in;

  in=fopen(name,"w");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  nitems=fprintf(in,"%s\n",comment[0]);
  if(nitems ==0) goto error;
  nitems=fprintf(in,"%s\n",comment[1]);
  if(nitems ==0) goto error;

  for(n=0;n<nndes;n++) {
    nitems=fprintf(in,"%d %g\n",n+1,buffer[n]);
    if(nitems ==0) goto error;
    }

  fclose(in);
  return (0);

 error:
  fclose(in);
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_saver1(const char *name, int nndes, double *buffer, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,nitems;
  FILE *in;

  in=fopen(name,"w");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

//   nitems=fprintf(in,"%s\n",comment[0]);
//   if(nitems ==0) goto error;
//   nitems=fprintf(in,"%s\n",comment[1]);
//   if(nitems ==0) goto error;
  if(comment!=0) {
    nitems=fprintf(in,"%s\n",comment[0]);
    if(nitems ==0) goto error;
    nitems=fprintf(in,"%s\n",comment[1]);
    if(nitems ==0) goto error;
    }
  else {
    nitems=fprintf(in,"no comment\n");
    if(nitems ==0) goto error;
    nitems=fprintf(in,"no comment\n\n");
    if(nitems ==0) goto error;
    }

  for(n=0;n<nndes;n++) {
    nitems=fprintf(in,"%d %lf\n",n+1,buffer[n]);
    if(nitems ==0) goto error;
    }

  fclose(in);
  return (0);

 error:
  fclose(in);
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_savec1(const char *name, int nndes, fcomplex *buffer, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,nitems;
  float a,p;
  FILE *in;

  in=fopen(name,"w");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  nitems=fprintf(in,"%s\n",comment[0]);
  if(nitems ==0) goto error;
  nitems=fprintf(in,"%s\n",comment[1]);
  if(nitems ==0) goto error;

  for(n=0;n<nndes;n++) {
/*
    a=sqrt(buffer[n].r*buffer[n].r+buffer[n].i*buffer[n].i);
    p=-atan2(buffer[n].i,buffer[n].r)*r2d;
*/
    a=abs(buffer[n]);
    p=-arg(buffer[n])*r2d;
    nitems=fprintf(in,"%d %f %f \n",n+1,a,p);
    if(nitems ==0) goto error;
    }

  fclose(in);
  return (0);

 error:
  fclose(in);
  return (-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_savec1(const char *name, int nndes, float *a, float *G, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,nitems;
  FILE *in;

  in=fopen(name,"w");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  nitems=fprintf(in,"%s\n",comment[0]);
  if(nitems ==0) goto error;
  nitems=fprintf(in,"%s\n",comment[1]);
  if(nitems ==0) goto error;

  for(n=0;n<nndes;n++) {
    nitems=fprintf(in,"%d %f %f \n",n+1,a[n],G[n]);
    if(nitems ==0) goto error;
    }

  fclose(in);
  return (0);

 error:
  fclose(in);
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_savec1(const char *name, int nndes, complex<double> *buffer, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,nitems;
  float a,p;
  FILE *in;

  in=fopen(name,"w");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  nitems=fprintf(in,"%s\n",comment[0]);
  if(nitems ==0) goto error;
  nitems=fprintf(in,"%s\n",comment[1]);
  if(nitems ==0) goto error;

  for(n=0;n<nndes;n++) {
/*
    a=sqrt(buffer[n].r*buffer[n].r+buffer[n].i*buffer[n].i);
    p=-atan2(buffer[n].i,buffer[n].r)*r2d;
*/
    a=abs(buffer[n]);
    p=-arg(buffer[n])*r2d;
    nitems=fprintf(in,"%d %lf %lf \n",n+1,a,p);
    if(nitems ==0) goto error;
    }

  fclose(in);
  return (0);

 error:
  fclose(in);
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_saver2(const char *name, int nndes, float *bufferx, float *buffery,char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,nitems;
  FILE *in;

  in=fopen(name,"w");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  nitems=fprintf(in,"%s\n",comment[0]);
  if(nitems ==0) goto error;
  nitems=fprintf(in,"%s\n",comment[1]);
  if(nitems ==0) goto error;

  for(n=0;n<nndes;n++) {
    nitems=fprintf(in,"%d %f %f\n",n+1,bufferx[n],buffery[n]);
    if(nitems ==0) goto error;
    }

  fclose(in);
  return (0);

 error:
  fclose(in);
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int quoddy_saver2(const char *name, int nndes, double *bufferx, double *buffery,char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,nitems;
  FILE *in;

  in=fopen(name,"w");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  nitems=fprintf(in,"%s\n",comment[0]);
  if(nitems ==0) goto error;
  nitems=fprintf(in,"%s\n",comment[1]);
  if(nitems ==0) goto error;

  for(n=0;n<nndes;n++) {
    nitems=fprintf(in,"%d %lf %lf\n",n+1,bufferx[n],buffery[n]);
    if(nitems ==0) goto error;
    }

  fclose(in);
  return (0);

 error:
  fclose(in);
  return (-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int quoddy_savec2_template(const char *name, int nndes, complex<T> *bufferx, complex<T> *buffery, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,nitems;
  T a,p;
  FILE *in;

  in=fopen(name,"w");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  nitems=fprintf(in,"%s\n",comment[0]);
  if(nitems ==0) goto error;
  nitems=fprintf(in,"%s\n",comment[1]);
  if(nitems ==0) goto error;

  for(n=0;n<nndes;n++) {
/*     nitems=fprintf(in,"%d %f %f %f %f \n",n,bufferx[n].r,bufferx[n].i,buffery[n].r,buffery[n].i); */
/*
    a=sqrt(bufferx[n].r*bufferx[n].r+bufferx[n].i*bufferx[n].i);
    p=-atan2(bufferx[n].i,bufferx[n].r)*r2d;
*/
    a=abs(bufferx[n]);
    p=-arg(bufferx[n])*r2d;
    nitems=fprintf(in,"%d %f %f ",n+1,a,p);
    if(nitems ==0) goto error;
/*
    a=sqrt(buffery[n].r*buffery[n].r+buffery[n].i*buffery[n].i);
    p=-atan2(buffery[n].i,buffery[n].r)*r2d;
*/
    a=abs(buffery[n]);
    p=-arg(buffery[n])*r2d;
    nitems=fprintf(in,"%f %f \n",a,p);
    if(nitems ==0) goto error;
    }

  fclose(in);
  return (0);

 error:
  fclose(in);
  return (-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int quoddy_savec2_template(const char *name, int nndes, T *a_u, T *G_u, T *a_v, T *G_v, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,nitems;
  FILE *in;

  in=fopen(name,"w");
  if(in == NULL) {
    printf("unable to open %s\n",name);
    return (-1);
    }

  nitems=fprintf(in,"%s\n",comment[0]);
  if(nitems ==0) goto error;
  nitems=fprintf(in,"%s\n",comment[1]);
  if(nitems ==0) goto error;

  for(n=0;n<nndes;n++) {
    nitems=fprintf(in,"%d %f %f ",n+1,a_u[n],G_u[n]);
    if(nitems ==0) goto error;
    nitems=fprintf(in,"%f %f \n",n+1,a_v[n],G_v[n]);
    if(nitems ==0) goto error;
    }

  fclose(in);
  return (0);

 error:
  fclose(in);
  return (-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int quoddy_savec2(const char *name, int nndes, fcomplex *bufferx, fcomplex *buffery, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=quoddy_savec2_template(name,nndes,bufferx,buffery,comment);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int quoddy_savec2(const char *name, int nndes, dcomplex *bufferx, dcomplex *buffery, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=quoddy_savec2_template(name,nndes,bufferx,buffery,comment);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int quoddy_savec2(const char *name, int nndes, float *a_u, float *G_u, float *a_v, float *G_v, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=quoddy_savec2_template(name,nndes,a_u,G_u,a_v,G_v,comment);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int quoddy_savec2(const char *name, int nndes, double *a_u, double *G_u, double *a_v, double *G_v, char **comment)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=quoddy_savec2_template(name,nndes,a_u,G_u,a_v,G_v,comment);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cefmo_loadc1(char *name, int nndes, fcomplex *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,nitems;
  char  line[64], *s;
//   float a,p;
  FILE *in;
  double z1,z2;

  in=fopen(name,"r");
  if(in==0)
    TRAP_ERR_RETURN(errno,1,"fopen(\"%s\",\"r\") error (%d %s)\n",name,errno,strerror(errno));

  do {
    fgets(line,1024,in);
    if(feof(in)) break;
    s=strstr(line,"#FILE HEADER END");
    } while(s==0);

  for(n=0;n<nndes;n++) {
    fgets(line,1024,in);
    do {
      s=strstr(line,"(");
      if(s!=0) s[0]=' ';
      } while(s!=0);
    do {
      s=strstr(line,",");
      if(s!=0) s[0]=' ';
      } while(s!=0);
    do {
      s=strstr(line,")");
      if(s!=0) s[0]=' ';
      } while(s!=0);
    do {
      s=strstr(line,"D");
      if(s!=0) s[0]='e';
      } while(s!=0);
    nitems=sscanf(line,"%lf %lf",&z1,&z2);
    buffer[n]=fcomplex((float) z1, (float) z2);
//    buffer[n]=polar<float>(a,-p*d2r);
    }

  fclose(in);
  return (0);
}
