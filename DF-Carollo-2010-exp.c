// Método das Diferenças Finitas Aplicado à problemas de Transferência de Calor

// Informações sobre os dados experimentais e parâmetros retirados de:

// Carollo, L.F.S. Estimação simultânea de propriedades termofísicas de materiais metálicos.
// Dissertação (Mestrado em Engenharia Mecânica). Universidade Federal de Itajubá.
// Itajubá, MG, 2010.

// Material: Aço Inox AISI 304 (Seção 5.1).
// Dados experimentais de temperatura retirados da Figura 5.4 (Y2).
// Dados para delta_x, L, delta_t, tempo total de experimento retirados da Tabela 3.1.
// Número de nós espaciais igual 10 e fluxo de calor retirados do texto da seção 5.1 e da Tabela 3.1.

// Configuração do Método Gauss-Seidel retirados da seção 3.2.

// Temperatura inicial retirada do experimento 26, da Tabela 5.2.
// Valores dos parâmetros de interesse retirados do experimento 26, da Tabela 5.3.

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define NX 10 //Número de nós da malha espacial
#define NT 801 //número de nós da malha temporal

///Termo fonte de acordo com o tempo de simulação
int termof (double tempo){
    int f=0;
    if(tempo >= 0 & tempo <= 20)
        f = 2640;
    else if(tempo > 20 & tempo <= 140)
        f = 640;
    else
        f = 0;
    return f;
    }

void main(){


//Dados experimentais do tempo
double t_exp[] = {0.00, 4.17, 6.67, 9.17, 11.25, 13.33, 15.00, 17.92, 19.58,
  20.83, 23.33, 24.58, 27.08, 30.83, 33.75, 37.92, 41.67, 46.67, 51.67, 57.08,
  62.50, 67.08, 72.50, 77.50, 82.50, 88.75, 93.75, 98.75, 102.50, 107.92,
  112.92, 118.33, 123.33, 127.92, 133.33, 139.17, 142.92, 148.33, 153.75,
  159.17};

//Dados experimentais da temperatura
double float_exp[] = {18.84, 18.85, 18.91, 19.02, 19.15, 19.26, 19.38, 19.49,
  19.63 , 19.76, 19.88, 19.98, 20.08, 20.16, 20.23, 20.29, 20.35, 20.43, 20.51,
  20.59, 20.65, 20.73, 20.80, 20.89, 20.98, 21.05, 21.13, 21.18, 21.25, 21.33,
  21.41, 21.49, 21.54, 21.63, 21.70, 21.79, 21.86, 21.90, 21.88, 21.88};

 //Dados/parâmetros do modelo matemático
 double lambda = 14.611;
 double pc = 3.9073*pow(10,6);
 double alpha = pc/lambda;
 double T_inicial = 18.84; //Temperatura inicial retirada do experimento 26 (Tabela 5.2)
 double xf = 10.88*pow(10,-3); //Comprimento (altura) da placa dada em metros
 double tf = 160;         //Tempo do experimento (tempo final)

 int f;
 double t[NT];

 double T_num[NT] = {0.0};
 double T[NX] = {0.0};
 double Tf[NX] = {0.0};
 double B[NX] = {0.0};
 double M[NX][NX] = {0.0};

 double dx = xf/(NX-1); //Delta x
 double dt = tf/(NT-1);   //Delta t
 int i, j;

 FILE *arquivo;
 FILE *arquivoa;
 FILE *arquivob;

 arquivo=fopen("dados-carollo.txt","w"); //Abrindo o arquivo para gravar os resultados
 arquivoa=fopen("todos-explicito.txt","w"); //Abrindo o arquivo para gravar os resultados
 arquivob=fopen("explicito-10-800.txt","w");

 //Gerando arquivo com os resultados obtidos por Carollo
  for(int i=0; i<40; i++){
    fprintf(arquivo,"%16.10f, %16.10f\n", t_exp[i], float_exp[i]);
  }

 //Vetor com a temperatura inicial do problema
 for(i=0; i<NX;i++){
    T[i] = T_inicial;
  }//for

  //Construção da matriz tri-diagonal de diferenças finitas
  for(i=0; i<NX; i++){
      for(j=0; j<NX; j++){
          if(i==j){
              M[i][j-1] = dt/(alpha*dx*dx);
              M[i][j] =1-((2*dt)/(alpha*dx*dx));
              M[i][j+1] = dt/(alpha*dx*dx);
          }//if
      }//for
  }//for
  M[0][0] = 1-((2*dt)/(alpha*dx*dx));
  M[0][1] = (2*dt)/(alpha*dx*dx);
  M[NX-1][NX-2] = (2*dt)/(alpha*dx*dx);
  M[NX-1][NX-1] = 1-((2*dt)/(alpha*dx*dx));

  //imprimindo matriz M
  printf("\n\n\n");
  printf("===========================\n");
  printf("Imprimindo matriz tri-diagonal M:\n");
  printf("===========================\n");

   for (i=0; i<NX; i++){
     for (j=0; j<NX; j++){
       printf("%7.5f\t", M[i][j]);
     }
     printf("\n");
   }
    t[0]=0;
    T_num[0] = T_inicial;
   //funcao geral
   for(i=1;i<NT;i++){
     t[i]=i*dt;
///Também necessitaremos calcular o termo independente B. Mas apenas o primeiro, pois sobre a influência do termo fonte.
     B[0]=((2*dt)/(alpha*dx*lambda))*termof((i-1)*dt);
     Tf[0]=(1-((2*dt)/(alpha*dx*dx)))*T[0]+((2*dt)/(alpha*dx*dx))*T[1]+((2*dt)/(alpha*dx*lambda))*termof(i*dt);//(Eq. (2.21).
     for(j=1;j<NX-1;j++){ //equacao 2.19 da dissertacao
        Tf[j]=(dt/(alpha*dx*dx))*T[j-1]+(1-((2*dt)/(alpha*dx*dx)))*T[j]+(dt/(alpha*dx*dx))*T[j+1];
     }//for
     Tf[NX-1]=((2*dt)/(alpha*dx*dx))*T[NX-2]+(1-((2*dt)/(alpha*dx*dx)))*T[NX-1];//Eq. (2.23).

     for(j=0;j<NX;j++){
        T[j] = Tf[j];
     }
     T_num[i] = Tf[NX-1];
   }//for

   //imprimindo matriz B
   printf("\n\n\n");
   printf("===========================\n");
   printf("Imprimindo vetor B:\n");
   printf("===========================\n");

    for (i=0; i<NX; i++){
        printf("%7.5f\n", B[i]);
    }

    ///Temperaturas no contorno em x = L, onde são coletados os dados experimentais.
    printf("\n\n\n");
    printf("===========================\n");
    printf("Imprimindo vetor T_num:\n");
    printf("===========================\n");

     for (i=0; i<NT; i++){
         printf("(linha %d)= %f\n ",i,T_num[i]);
         fprintf(arquivoa,"T_num[%d] = %16.10f\n", i, T_num[i]);
         fprintf(arquivob,"%f, %16.10f\n", t[i], T_num[i]);
      }

      //imprimindo matriz Tf
      printf("\n\n\n");
      printf("===========================\n");
      printf("Imprimindo vetor Tf:\n");
      printf("===========================\n");

       for (i=0; i<NX; i++){
           printf("%7.5f\n", Tf[i]);
       }

       fclose(arquivo);
       fclose(arquivoa);
       fclose(arquivob);

       FILE *arquivoc;
       arquivoc=fopen("grafico-explicito.txt","w"); //Abrindo o arquivo para gravar os resultados
       fprintf(arquivoc,"reset\n");
       fprintf(arquivoc,"set size 0.5\n");
       fprintf(arquivoc,"set xlabel 'tempo(s)'\n");
       fprintf(arquivoc,"set ylabel 'Temperatura(°C)'\n");
       fprintf(arquivoc,"set grid\n");
       fprintf(arquivoc,"set key right top\n");
       //fprintf(arquivoc,"set key left bottom\n");
       fprintf(arquivoc,"set title \"Formulação explícita com  " "{/Symbol D}x = 0,001209\n");
       //fprintf(arquivoc,"plot 'explicito.txt' using ($1):($2) t\"metodo explicito\" with linespoints pointinterval 10\n");
    //   fprintf(arquivoc,"plot 'explicito-10-1600.txt' using ($1):($2) t\" {/Symbol D}t = 0,1\" with linespoints pointinterval 100\n");
       fprintf(arquivoc,"plot 'explicito-10-800.txt' using ($1):($2) t\" {/Symbol D}t = 0,2\" with linespoints pointinterval 100\n");
    //   fprintf(arquivoc,"replot 'explicito-10-320.txt' using ($1):($2) t\"{/Symbol D}t = 0,5\" with linespoints pointinterval 100\n");
    //   fprintf(arquivoc,"replot 'explicito-10-160.txt' using ($1):($2) t\"{/Symbol D}t = 1\" with linespoints pointinterval 100\n");
    //   fprintf(arquivoc,"replot 'explicito-10-80.txt' using ($1):($2) t\"{/Symbol D}t = 2\" with linespoints pointinterval 100\n");
    //   fprintf(arquivoc,"replot 'explicito-10-32.txt' using ($1):($2) t\"{/Symbol D}t = 5\" with linespoints pointinterval 100\n");
    //   fprintf(arquivoc,"replot 'explicito-10-16.txt' using ($1):($2) t\"{/Symbol D}t = 10\" with linespoints pointinterval 100\n");
       fprintf(arquivo, "replot 'dados-carollo.txt' using ($1):($2) t\"Carollo\" with linespoints pointinterval 10\n");
       fprintf(arquivoc,"replot\n;pause -1");
       fclose(arquivoc);

       system("gnuplot grafico-explicito.txt");

}//main
