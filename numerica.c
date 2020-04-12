#include <stdio.h>
#include <math.h>
#include<stdlib.h>

#define MAX 100

//////////PROGRAMA PARA RESOLUÇÃO DE EQUAÇÕES DIFERENCIAVEIS//////////////////////////////////////////////////////////

//////////////////////////////////////////FUNÇÕES PROFESSOR///////////////////////////////////////////////////////////

int ler_ficheiro(float *polinomio)
{
    FILE *fp;
    float coeff;
    int grau, i;

    fp = fopen("polinomio.txt", "r");
    if (fp == NULL)
    {
        printf("Erro ao abrir o ficheiro!\n");
        return -1;
    }

    fscanf(fp, "%d", &grau);

    for(i=0; i<=grau; i++)
    {
        fscanf(fp, "%f", &coeff);
        polinomio[i] = coeff;
    }
    fclose(fp);
    return grau;
}

float calcula_polinomio(float *polinomio, int grau, float x)
{
    int i;
    float resultado = 0;

    for(i=0; i<=grau; i++)
    {
        resultado += polinomio[i] * pow(x, grau-i);
    }
    return resultado;
}

//////////////Funções para calcular derivadas//////////////////////////////////////////

float derivada_polinomio(float *polinomio, int grau, float x)
{
    int i;
    float result = 0;
    for (i = 0; i <= grau; i++)
        result += (grau-i)*polinomio[i]*pow(x, (grau-i)-1);
    return result;
}

//////////////////////////////////FUNÇÃO PARA APLICAR O MÉTODO DE NEWTON RAPHSON/////////////////////////////////////

float newton_raphson(float *polinomio, int grau, float x, float e)
{
    float f, df, result;
    f=calcula_polinomio(polinomio, grau, x);
    df=derivada_polinomio(polinomio, grau, x);
    result = x - (f/df);
    while(fabsf(result - x)>=e)
    {
        x=result;
        f=calcula_polinomio(polinomio, grau, x);
        df=derivada_polinomio(polinomio, grau, x);
        result = x - (f/df);
    }
    return result;
}


//////////////////////////PROGRAMA PARA INTERPOLAÇÃO POLINOMIAL//////////////////////////////////////////////////


int read_file(float *x, float *y)
{
    FILE *f;
    int size;
    float val;
    f = fopen("values.txt", "r");
    if (f == NULL)
    {
        printf("Erro ao abrir o ficheiro, verifique o ficheiro e volte a correr o programa\n");
        return 0;
    }
    fscanf(f, "%d", &size);

    for (int i = 0; i < size; i++)
    {
        fscanf(f, "%f", &val);
        x[i] = val;
    }

    for (int i = 0; i < size; i++)
    {
        fscanf(f, "%f", &val);
        y[i] = val;
    }
    fclose(f);
    return size;
}

//Verifica o menor valor no intervalo x
float min(float *x, int n)
{
    float result = x[0];
    for (int i = 0; i < n; i++)
        if (result > x[i])
            result = x[i];
    return result;
}

//Verifica o maior valor no intervalo x
float max(float *x, int n)
{
    float result = x[0];
    for (int i = 1; i < n; i++)
        if (result < x[i])
            result = x[i];
    return result;
}

float lagrange(float *y, float *x, int n, float xp)
{
    float P = 0;
    for (int i=0;i<n;i++)
    {
        float L=1;
        for (int j=0;j<n;j++)
        {
            if (i != j)
            {
                L = L * (xp-x[j])/(x[i]-x[j]);
            }
        }
        P = P + y[i]*L;
    }
    return P;
}

void dif(float *x, float *y, int n, float *np)
{
    float d[n][n];
    for (int i = 0; i < n; i++)
    {
        d[i][0] = y[i];
    }
    for (int j = 1; j < n; j++)
    {
        for (int i = j; i < n; i++)
        {
            d[i][j] = (d[i][j - 1] - d[i - 1][j - 1]) / (x[i] - x[i - j]);
            if(i==j)
            {
                np[i] = d[i][j];//Coloca os valores das diferencas divididas em um array
            }
        }
    }
}

float erro(float *x, float *d, int n, float xp)
{
    float result=0;
    float e=1;
    for (int i = 0; i < n; i++)
    {
        e*=(xp-x[i]);
    }
    result=fabs(e*d[n-2]);
    return result;
}

////////////////////////////////////////APLICAÇÃO DOS MENUS////////////////////////////////////////////////////////


void print_menu_eq()
{
    float *polinomio, x, nr, e;
    polinomio=(float *) malloc(MAX * sizeof(float));
    int grau;
    grau=ler_ficheiro(polinomio);
    system("cls");
    printf("O polinomio em analise e: \n");
    for (int i = 0; i <= grau; i++)
    {
        printf("+(%f).x^%d", polinomio[i], grau-i);
    }

    printf("\nDefina um ponto x: ");
    scanf("%f", &x);
    printf("\nDefina o erro: ");
    scanf("%f", &e);
    nr = newton_raphson(polinomio, grau, x, e);
    printf("\nO valor de uma das raizes e %f\n\nClique R para recomecar\n\nClique S para sair", nr);
}

void menu_equacao()
{
    char op;
    system("cls");
    printf("\n1 - Metodo de Newton Raphson\n\nS - Saida");
    op=toupper(getch());
        if (op == '1')
        {
            system("cls");
            do {
                print_menu_eq();
                op=toupper(getch());
                if(op=='S')
                    menu_equacao();
                if(op=='R')
                {
                 system("cls");
                 print_menu_eq();
                 op=toupper(getch());
                }
            }while(op != 'S');
        }
        if (op == 'S')
            menu();
}

void menu_inter()
{
    float *x, *y, xp, np[MAX], e, polagran;
    x=(float *) malloc(MAX* sizeof(float));
    y=(float *) malloc(MAX* sizeof(float));
    int size=read_file(x,y);
    float d[MAX];
    system("cls");
    printf("\t\t\tINTERPOLACAO POLINOMIAL\n\nTabela com valores de x e y respectivamente:  \n");
    for (int i = 0; i < size; i++)
    {
            printf("%f    %f\n", x[i], y[i]);
    }
    printf("Defina um valor um novo x: \n");
    scanf("%f", &xp);
        if(xp > max(x,size) || xp < min(x,size))
        {
            system("cls");
            printf("O valor de x tem que esta dentro do intervalo, escolha um novo valor para x\n\nAperte qualquer tecla para voltar");
            getch();
            menu_inter();
        }
        if(xp < max(x, size) && xp > min(x, size))
        {
            char op;
            polagran = lagrange(y, x, size, xp);
            printf("O polinomio interpolador de Lagrange e: \n");
            for (int i = 0; i < size; i++)
            {
                printf("(%f)*L%d(x) +", y[i], i);
            }
            printf("\nO resultado para x=[%f] aproxima-se de %f\n", xp, polagran);

            dif(x, y, size, np);
            for (int i = 1; i < size; i++) {//Passa os valores para um novo array iniciando no 0 e o ultimo passará a ser 0
                d[i - 1] = np[i];
            }
            for (int i = 0; i < size - 1; i++) {
                printf("%f ", d[i]);
            }
            e = erro(x, d, size, xp);
            printf("\nO erro de aproximacao pelas diferencas divididas e: %f", e);

            printf("\n\nClique R para recomecar\n\nClique S para voltar");
                op = toupper(getch());
                if (op == 'R')
                    menu_inter();
                if (op == 'S')
                    menu();
        }
}

void menu()
{
    char op;
    do{
    system("cls");
    printf("\t\tTRABALHO ANALISE NUMERICA\n\n1 - Equacoes nao lineares\n\n2 - Interpolacao polinomial\n\n3 - Aproximacao de funcoes\n\n4 - Equacoes diferenciais\n\nI - Instrucoes de uso\n\nS - Sair");
    op=toupper(getch());

        if (op == '1')
            menu_equacao();
        if (op == '2')
            menu_inter();
        if (op == '3') {
            system("cls");
            printf("To be continued...\n\nAperte qualquer tecla para sair");
            getch();
            menu();
        }
        if (op == '4') {
            system("cls");
            printf("To be continued...\n\nAperte qualquer tecla para sair");
            getch();
            menu();
        }
        if (op == 'I') {
            system("cls");
            printf("Para o '2' devera ter um ficheiro de texto com o nome de values.txt no mesmo diretorio do programa on o ficheiro funcionara de tal forma:\n");
            printf("O primeiro valor será um inteiro indicando o numero de x's e y's na tabela, a seguir devera ser expresso os valores de x's e depois os valores de y's");
            printf("\nO programa indicara o resto do caminho a seguir\n\nAperte qqr tecla para sair.");
            getch();
            menu();
        }
        if (op == 'S')
            exit(0);
    }while(op!='s' && op!='S');
}

int main()
{
    menu();
    return 0;
}