/*compile with $g++ -o homo-lumo homo-lumo-x.y.c */
/*to use:   $homo-lumo   or homo-lumo x       where x is the minimum occupancy of the lowest band (defaul is 0) */
/*or $homo-lumo 0 OUTCAR*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
int main(int argc, char *argv[])
{
FILE *outcar;
int i,j,k,s,nef,nbands,nkpt,band,kh,kl;
char str[200]="",ch;
float tol,gapup,gapdown;


if( argc > 1 )tol=atof(argv[1]);
else {tol=0;}

if( argc > 2 )strcat(str,argv[2]);
else {strcat(str,"OUTCAR");}


outcar = fopen(str,"r"); /* Arquivo ASCII, para leitura */
if(!outcar)
{
printf("Erro na abertura do arquivo %s\n",str);
exit(0);
}

nef=0;j=0;
while (fscanf(outcar,"%s",str) != EOF){
if(strcmp(str,"E-fermi")==0)nef++;                      /*verifica se o arquivo outcar esta completo e conta qtos e-fermi*/
if(strcmp(str,"Voluntary")==0)j++;
}

if(j == 0){printf("\n\nIncomplete OUTCAR!! bye! bye! \n\n");exit(0);} /*verifica se o calculo terminou corretamente*/

rewind(outcar);
do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra ISPIN*/
while(strcmp(str,"ISPIN")!=0);
fscanf(outcar,"%s",str);        /*pula o =*/
fscanf(outcar,"%d",&s);        

do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra NKPTS*/
while(strcmp(str,"NKPTS")!=0);
fscanf(outcar,"%s",str);        /*pula o =*/
fscanf(outcar,"%d",&nkpt);      /*lê o numero de kpts*/

float a,b[2*nkpt],c[2*nkpt],homo,lumo,occupb[2*nkpt],occupc[2*nkpt],kpoint[3],directgapup[2*nkpt],directgapdown[2*nkpt];

do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra NBANDS=*/
while(strcmp(str,"NBANDS=")!=0);

fscanf(outcar,"%d",&nbands);      /*lê o numero de bandas*/

/*-----------------------------spin----------------------------------------------------------------------------------------------------*/
if(s==2){

int sum1=0,sum2=0,bndh[2*nkpt],bndl[2*nkpt];

for(i=0;i<nef;i++){
do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra E-fermi*/
while(strcmp(str,"E-fermi")!=0);
}

for(i=0;i<5;i++){
do
ch = getc(outcar);           /*vai pra proxima linha*/
while(ch!='\n');
}



c[0]=0;
printf("           NKPT   HOMO         LUMO       LUMO-HOMO     bandHOMO    occupHOMO     bandLUMO    occupLUMO   | NKPT                   KPOINT\n");

for(k=0;k<nkpt;k++){

for(i=0;i<3;i++)fscanf(outcar,"%s",str);
for(i=0;i<3;i++)fscanf(outcar,"%f",&kpoint[i]); 

for(i=0;i<2;i++){
do
ch = getc(outcar);           /*vai pra proxima linha*/
while(ch!='\n');
}

bndh[k]=bndl[k]=i=0;
do{
b[k]=c[k];occupb[k]=occupc[k];bndh[k]=bndl[k];
fscanf(outcar,"%f",&a);
fscanf(outcar,"%f",&c[k]);
fscanf(outcar,"%f",&a);occupc[k]=a;
bndl[k]++;i++;
}while(a>tol);

directgapup[k]=c[k]-b[k];

printf("spin-up   %3d   %f    %9f     %9f        %d         %.5f         %d         %.5f    | %3d     %7.4f       %7.4f       %7.4f\n",k+1,b[k],c[k],directgapup[k],bndh[k],occupb[k],bndl[k],occupc[k],k+1,kpoint[0],kpoint[1],kpoint[2]);



for(j=0;j<nbands-i+2;j++){
do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');
}
}

for(i=0;i<2;i++){
do
ch = getc(outcar);           /*vai pra proxima linha*/
while(ch!='\n');
}


for(k=0;k<nkpt;k++){

for(i=0;i<3;i++)fscanf(outcar,"%s",str);
for(i=0;i<3;i++)fscanf(outcar,"%f",&kpoint[i]);

for(i=0;i<2;i++){
do
ch = getc(outcar);           /*vai pra proxima linha*/
while(ch!='\n');
}

bndh[k+nkpt]=bndl[k+nkpt]=i=0;
do{
b[k+nkpt]=c[k+nkpt];occupb[k]=occupc[k];bndh[k+nkpt]=bndl[k+nkpt];
fscanf(outcar,"%f",&a);
fscanf(outcar,"%f",&c[k+nkpt]);
fscanf(outcar,"%f",&a);occupc[k]=a;
bndl[k+nkpt]++;i++;
}while(a>tol);

directgapdown[k+nkpt]=c[k+nkpt]-b[k+nkpt];

printf("spin-down %3d   %f    %9f     %9f        %d         %.5f         %d         %.5f    | %3d     %7.4f       %7.4f       %7.4f\n",k+nkpt+1,b[k+nkpt],c[k+nkpt],directgapdown[k+nkpt],bndh[k+nkpt],occupb[k],bndl[k+nkpt],occupc[k],k+nkpt+1,kpoint[0],kpoint[1],kpoint[2]);



for(j=0;j<nbands-i+2;j++){
do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');
}
}



printf("\n"); 
printf("   HOMO(k-point)           LUMO(k-point)         homoBAND       lumoBAND     direct gap up (k-point)   direct gap(k-point)\n");
homo=b[0];
lumo=c[0];
kh=kl=0;

for(k=0;k<2*nkpt;k++){
if(homo<b[k]){homo=b[k];gapdown=c[k]-b[k];kh=k;}
if(lumo>c[k]){lumo=c[k];kl=k;}
}

int gup;
gapup=20;
for(j=0;j<nkpt;j++){
if(directgapup[j]<gapup){gapup=directgapup[j];gup=j;}
}

int gdown;
gapdown=20;
for(j=nkpt;j<2*nkpt;j++){
if(directgapdown[j]<gapdown){gapdown=directgapdown[j];gdown=j;}
}


if(kh<=nkpt & kl<=nkpt)printf("%f (%d)up        %9f (%d)up          %d             %d           %9f (%d)           %9f (%d)\n\n",homo,kh+1,lumo,kl+1,bndh[kh],bndl[kl],gapup,gup+1,gapdown,gdown+1);
if(kh>nkpt & kl>nkpt)printf("%f (%d)down      %9f (%d)down        %d             %d           %9f (%d)           %9f (%d)\n\n",homo,kh+1,lumo,kl+1,bndh[kh],bndl[kl],gapup,gup+1,gapdown,gdown+1);
if(kh<=nkpt & kl>nkpt)printf("%f (%d)up        %9f (%d)down        %d             %d           %9f (%d)           %9f (%d)\n\n",homo,kh+1,lumo,kl+1,bndh[kh],bndl[kl],gapup,gup+1,gapdown,gdown+1);
if(kh>nkpt & kl<=nkpt)printf("%f (%d)down      %9f (%d)up          %d             %d           %9f (%d)           %9f (%d)\n\n",homo,kh+1,lumo,kl+1,bndh[kh],bndl[kl],gapup,gup+1,gapdown,gdown+1);


printf("LUMO-HOMO = %f (%d-%d)\n\n",lumo-homo,kl+1,kh+1);


for(k=0;k<nkpt-1;k++){if(bndh[k]!=bndh[k+1])sum1++;}
if(sum1==0)printf("SEMICONDUCTOR SPIN-UP!\n");
else {printf("METAL SPIN-UP!\n");}

for(k=nkpt;k<2*nkpt-1;k++){if(bndh[k]!=bndh[k+1])sum2++;}
if(sum2==0)printf("SEMICONDUCTOR SPIN-DOWN!\n");
else {printf("METAL SPIN-DOWN!\n");}

}

/*------------------------no-spin-------------------------------------------------------------------------------------------------------*/
else{

int sum=0,bndh[nkpt],bndl[nkpt];

for(i=0;i<nef;i++){
do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra E-fermi*/
while(strcmp(str,"E-fermi")!=0);
}

for(i=0;i<3;i++){
do
ch = getc(outcar);           /*vai pra proxima linha*/
while(ch!='\n');
}

c[0]=0;

printf("NKPT     HOMO         LUMO        LUMO-HOMO     bandHOMO      occupHOMO     bandLUMO     occupLUMO     | NKPT                   KPOINT\n");
for(k=0;k<nkpt;k++){

for(i=0;i<3;i++)fscanf(outcar,"%s",str);
for(i=0;i<3;i++)fscanf(outcar,"%f",&kpoint[i]);

for(i=0;i<2;i++){
do
ch = getc(outcar);           /*vai pra proxima linha*/
while(ch!='\n');
}

bndh[k]=bndl[k]=i=0;
do{
b[k]=c[k];occupb[k]=occupc[k];bndh[k]=bndl[k];
fscanf(outcar,"%f",&a);
fscanf(outcar,"%f",&c[k]);
fscanf(outcar,"%f",&a);occupc[k]=a;
bndl[k]++;i++;
}while(a>tol);

directgapup[k]=c[k]-b[k];

printf("%3d   %9f    %9f     %9f        %d         %.5f         %d         %.5f      | %3d     %7.4f       %7.4f       %7.4f\n",k+1,b[k],c[k],directgapup[k],bndh[k],occupb[k],bndl[k],occupc[k],k+1,kpoint[0],kpoint[1],kpoint[2]);

for(j=0;j<nbands-i+2;j++){
do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');
}
}

printf("\n");
printf("HOMO(k-point)        LUMO(k-point)        homoBAND     lumoBAND     direct gap(k-point)\n");
homo=b[0];
lumo=c[0];
kh=kl=0;

for(k=0;k<nkpt;k++){
if(homo<b[k]){homo=b[k];kh=k;}
if(lumo>c[k]){lumo=c[k];kl=k;}
}

int g;
gapup=20;
for(j=0;j<nkpt;j++){
if(directgapup[j]<gapup){gapup=directgapup[j];g=j;}
}





printf("%f (%d)        %9f (%d)           %d            %d           %9f (%d)\n\n",homo,kh+1,lumo,kl+1,bndh[kl],bndl[kh],gapup,g+1);

printf("LUMO-HOMO = %f (%d-%d)\n\n",lumo-homo,kl+1,kh+1);


for(k=0;k<nkpt-1;k++){if(bndh[k]!=bndh[k+1])sum++;}
if(sum==0)printf("SEMICONDUCTOR!\n");
else {printf("METAL!\n");}




}

fclose(outcar);

printf("\n\n------------------------------------------------------\n");
printf("Homo-Lumo energy for VASP results\n");
printf("by F. Matusalem - 03/2013 - filipematus@gmail.com\n");
printf("Version 5.6 10/2016\n\n");
   
}
