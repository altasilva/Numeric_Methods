/****************************************************PROBLEMAS INVERSOS EM ENGENHARIA 2015-01**************************************************/
/* Implementacao do algoritmo de Luus-Jaakola                                                                                                 */
/* Autor: Wagner Rambaldi Telles                                                                                                              */
/* Data: 30 de abril 2015                                                                                                                     */
/**********************************************************************************************************************************************/

//Bibliotecas
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
//Variaveis Globais
#define QTDVAR 2   //Numero de Vari�veis do problema
///=============================================================
//Manipular
#define NOUT 2000  //Numero de iteracoes externas
#define NINT 1000  //Numero de iteracoes internas
#define EPS 0.20  //Fator de reducao do intervalo de busca
///=============================================================
#define PI 3.1415 //Estimativa para pi

//Funcoes utilizadas
double f(double []);
double randomico();


//============================================================Funcao Objetivo==================================================================//
//Defindo a funcao objetivo a ser usada no algoritmo
double f(double x[]){
   double f, S1, S2, S11, S12, S21, S22;

   //2.3 - Rotated hyper-ellipsoid function
   // NOUT=2000     NINT=100    EPS=0.20  intervalo [-500 a 500]
   //f=pow(x[0], 2) + (pow(x[0], 2) + pow(x[1], 2));

   //2.4 - Rosenbrock's valley / banana function / second function of De Jong
   // NOUT=2000     NINT=1000    EPS=0.20  intervalo [-5 a 10]
   //f=100*(pow(x[1] - pow(x[0], 2), 2)) + pow(1-x[0], 2);

   //2.6 - Schwefel's function  nout=2000   nint=1000    eps=0.20  intervalo [0-500]
   //f = - x[0] * sin(sqrt(fabs(x[0]))) - x[1] * sin(sqrt(fabs(x[1])));

   //2.11- Michalewicz's function  nout=2000   nint=1000    eps=0.20  intervalo [0-PI] para m=10
   //double a = (pow(1*x[0],2)/PI);
   //double b = (pow(2*x[1],2)/PI);
   //f = - (sin(x[0]) * (pow(sin(a),20)) + sin(x[1]) * (pow(sin(b),20)));

   //2.12- Branins's function  nout=20   nint=10    eps=0.20  intervalo [-5 a 10] e [0 a 15]
   f = pow((x[1]-5.1/(4*pow(PI,2))*pow(x[0],2)+(5/PI)*x[0]-6),2)+(10*(1-(1/(8*PI)))*cos(x[0]))+10;

//FIM DAS FUNÇÕES PARA O SEMINÁRIO
//**********************************************************************************

   //funcao Jong  -->  NOUT=10   NINT=5  EPS=0.20
  // f = pow(x[0], 2) + pow(x[1], 2);

   //funcao Restrigin   -->  NOUT=10   NINT=30  EPS=0.18
   //f = 20 + (pow(x[0], 2) - 10 * cos(2 * PI * x[0])) + (pow(x[1], 2) - 10 * cos (2 * PI * x[1]));

//   f = pow((x1 - 1), 2) + pow((x2 - 1), 2) - x1 * x2;                                               //Exerc�cio F1
//   f = pow(x1, 2) + 3 * pow(x2, 2) - 2 * x1 + 2 * x2;                                               //Exerc�cio F2
//   f = pow(x[0], 4) + pow(x[1], 4) - 3;                                                            //Trabalho F1
//   f = 20 + (pow(x[0], 2) - 10 * cos(2 * PI * x[0])) + (pow(x[1], 2) - 10 * cos (2 * PI * x[1]));  //Trabalho F2
//   S1 = (4 - (2.1 * (pow(x[0], 2))) + ((pow(x[0], 4)) / 3)) * (pow(x[0], 2));                       //Trabalho F3
//   S2 = (4 * (pow(x[1], 2)) - 4) * (pow(x[1], 2));                                                  //Trabalho F3
//   f = S1 + (x[0] * x[1]) + S2;                                                                     //Trabalho F3
//   S11 = pow((x1 + x2 + 1), 2);                                                                     //Trabalho F4
//   S12 = (19 - 14 * x1 + 3 * (pow(x1, 2)) - 14 * x2 + 6 * x1 * x2 + 3 * (pow(x2, 2)));              //Trabalho F4
//   S1 = 1 + S11 * S12;                                                                              //Trabalho F4
//   S21 = pow((2 * x1 - 3 * x2), 2);                                                                 //Trabalho F4
//   S22 = (18 - 32 * x1 + 12 * (pow(x1, 2)) + 48 * x2 - 36 * x1 * x2 + 27 * (pow(x2, 2)));           //Trabalho F4
//   S2 = 30 + S21 * S22;                                                                             //Trabalho F4
//   f = S1 * S2;                                                                                     //Trabalho F4
//    f = -((sin(x[0])*pow(sin((1*x[0]*x[0])/M_PI), 20)) + (sin(x[1])*pow(sin((2.0*x[1]*x[1])/M_PI), 20)));

   return (f);
}
//=============================================================================================================================================//

//===========================================================Funcao Randomico==================================================================//
//Funcao que gera um numero aleatorio entre [0, 1]
double randomico(){
   double aleatorio;

   aleatorio = ((double)rand() / (double)RAND_MAX);

   return (aleatorio);
}
//=============================================================================================================================================//

//======================================================Programa Principal=====================================================================//
int main(){
   int iteracao, i, j, k, cont=0, corrida=0, tot_corrida=100, aux=0, cont_max=200000;
   double r[QTDVAR] = {0}, min[QTDVAR] = {0}, max[QTDVAR] = {0}, xi[QTDVAR] = {0}, xf[QTDVAR] = {0},\
   result_x1[tot_corrida], result_x2[tot_corrida];
   double soma1=0, soma2=0, DP1, DP2, d, tol=1E-4;
   double x1=-PI, x2=12.275, media1=0.0, media2=0.0;

   int semente_inicial = 31; ///INICIALIZEI A SEMENTE
///   srand(time(NULL));


///=============================================================
///Manipular
//--------------------------------------------------------Intervalo de busca das funcoes-------------------------------------------------------//
   min[0] = -5; //0.0; //Valor minimo do intervalo de busca de x1.
   max[0] = 10; //M_PI; //Valor maximo do intervalo de busca de x1.
   min[1] = 0; //0.0; //Valor minimo do intervalo de busca de x2.
   max[1] = 15; //M_PI; //Valor maximo do intervalo de busca de x2.
//---------------------------------------------------------------------Fim---------------------------------------------------------------------//
///============================================================

   clock_t tInicio, tFim, tDecorrido;    //Calculando o tempo de execu��o do programa
   tInicio=clock();                      //Iniciando o tempo de execu��o do programa

   FILE *arquivo;
   arquivo=fopen("LuusJaakola.txt","w"); //Abrindo o arquivo para gravar os resultados
   fprintf(arquivo,"\t x inicial\t\t\t x final\n");             //Variaveis que serao
   fprintf(arquivo,"x1\t\t x2\t\t x1\t\t x2\t\t F(x1,x2)\t Iteracoes\t\n"); //impressas no arquivo-texto.
   FILE *arquivo1;
   arquivo1=fopen("fx_Melhor.txt","w"); //Abrindo o arquivo para gravar os resultados
   FILE *arquivo2;
   arquivo2=fopen("fx_Todos.txt","w"); //Abrindo o arquivo para gravar os resultados

//---------------------------------------------------------Gerando a estimativa inicial--------------------------------------------------------//
///PASSEI O INÍCIO DAS CORRIDAS AQUI.
///O PROGRAMA VAI REDUZINDO O INTERVALO DE BUSCA (R). ENTÃO AO INICIARMOS UMA NOVA CORRIDA, O VALOR DE R DEVE SER -50CALCULADO
///NOVAMENTE, ASSIM COMO INICIALIZAR O CONTADOR E TODAS AS OUTRAS VARIÁVE500.
  for (corrida=0; corrida<tot_corrida; corrida++){

   unsigned semente = semente_inicial;///PROCEDIMENTO PARA GERAR UM NÚMERO ALEATORIO BASEADO NA SEMENTE QUE DEFINI LOGO NO INICIO.
   srand(semente);
   cont = 0; ///INICIEI O CONTADOR DE NOVO.

   for(i = 0; i < QTDVAR; i++){
      r[i] = max[i] - min[i];                //Intervalo de busca da solucao
      xf[i] = min[i] + (r[i] * randomico()); //Estimativa Inicial
      fprintf(arquivo,"%f\t", xf[i]);
   }
  // fprintf(arquivo,"%16.10f\t %16.10f\t %9.10f\n", xf[0], xf[1], f(xf));
   fprintf(arquivo1,"%d\t%.5f\n", cont, f(xf));
   fprintf(arquivo2,"%d\t%.5f\n", cont, f(xf));
//---------------------------------------------------------------------Fim---------------------------------------------------------------------//
//*********************************************************Executando o algoritmo**************************************************************/

     for(i = 1; i <= NOUT; i++){                             //Numero de iteracoes externas do programa
        for(j = 1; j <= NINT; j++){                          //Numero de iteracoes internas
           for(k = 0; k < QTDVAR; k++){                      //Definindo o novo vetor solucao
              xi[k] = xf[k] + (r[k] * (-0.5 + randomico())); //Nova solucao
              if(xi[k] < min[k])                             //Definindo o limite minimo como nova solucao caso a mesma seja menor que min[k]
                 xi[k] = min[k];
              if(xi[k] > max[k])                             //Definindo o limite maximo como nova solucao caso a mesma seja maior que mmax[k]
                 xi[k] = max[k];
           }

           if(f(xi) < f(xf)){                                 //Se a nova estimativa for melhor que a anterior, atualiza os par�metros
             for(k = 0; k < QTDVAR; k++)
               xf[k] = xi[k];
           }

           fprintf(arquivo1,"%d\t%.5f\n", cont, f(xf));
           fprintf(arquivo2,"%d\t%.5f\n", cont, f(xi));
           cont++;


           if((cont>=cont_max) || ((fabs(xf[0]-x1) <= tol) && (fabs(xf[1]-x2) <= tol))){
             goto finaliza;
           }

        }//fim do for NINT

        for(k = 0; k < QTDVAR; k++){  //Reduzindo o intervalo de busca
          r[k] *= (1 - EPS);        //Novo intervalo de busca.
        }
     }//fim do for NOUT
//************************************************************Fim do algoritmo*****************************************************************//
   finaliza:
   result_x1[corrida]=xf[0]; //Salvando as solucoes de x1 em um vetor.
   result_x2[corrida]=xf[1]; //Salvando as solucoes de x2 em um vetor.
   media1+=result_x1[corrida]/tot_corrida; //Calculando a media das solucoes de x1.
   media2+=result_x2[corrida]/tot_corrida; //Calculando a media das solucoes de x2.
   fprintf(arquivo,"%f\t %f\t %f\t %d\n", result_x1[corrida], result_x2[corrida], \
   f(xf), cont); //Imprime os valores da melhor solucao.
  // fprintf(arquivo," x1 = %16.5f\n x2 = %16.5f\n f(x1, x2) = %9.5f\n", xf[0], xf[1], f(xf));
   semente_inicial++; 

 }//fim for corrida

   for (i=0; i<tot_corrida; i++){ //Gerando a media das solucoes.
       soma1+=pow((result_x1[i]-media1),2); //Somando as diferencas entre a media e as solucoes de x1.
       soma2+=pow((result_x2[i]-media2),2); //Somando as diferencas entre a media e as solucoes de x2.
   } //Fim (Gerando a media das solucoes).
   fprintf(arquivo,"\nMedia: %f\t %f\n", media1, media2);
   DP1=sqrt(soma1/(tot_corrida-1)); //Calculando o desvio-padrao de x1.
   DP2=sqrt(soma2/(tot_corrida-1)); //Calculando o desvio-padrao de x2.
   fprintf(arquivo,"\nDesvio Padrao: %f\t %f\n", DP1, DP2);
   d=sqrt(((pow(1+(DP1/media1),2))+(pow(1+(DP2/media2),2)))/2)-1;
   fprintf(arquivo,"\nDispersao: %f\n", d);

   tFim = clock();                                                                       //Fechando o tempo de execu��o do progrma e imprimindo
   tDecorrido = ((tFim - tInicio) / (CLOCKS_PER_SEC/1000));
   fprintf(arquivo,"\nTempo de execucao: %ld milisegundos.\n", tDecorrido);

   fclose(arquivo);
   fclose(arquivo1);
   fclose(arquivo2);

   FILE *arquivo3;
   arquivo3=fopen("graficoerro.txt","w"); //Abrindo o arquivo para gravar os resultados
   fprintf(arquivo3,"reset\n");
   fprintf(arquivo3,"set xlabel 'Numero de Avaliacoes'\n");
   fprintf(arquivo3,"set ylabel 'f(x1, x2)'\n");
   fprintf(arquivo3,"set grid\n");
   fprintf(arquivo3,"plot 'fx_Todos.txt' using ($1):($2) t\"Todos\" with linespoints pointinterval 3\n");
   //fprintf(arquivo3,"plot 'fx_Melhor.txt' using ($1):($2) t\"Melhor\" with linespoints pointinterval 3\n");
   fprintf(arquivo3,"replot 'fx_Melhor.txt' using ($1):($2) t\"Melhor\" with linespoints pointinterval 10\n");
   fprintf(arquivo3,"replot\n;pause -1");
   fclose(arquivo3);

   system("gnuplot graficoerro.txt");
//   system ("pause");
   return (0);
} //fecha main
//=======================================================Fim do Programa Principal=============================================================//
