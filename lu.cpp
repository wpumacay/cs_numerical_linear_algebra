#include<iostream>
#include<iomanip>
using namespace std;


class Matriz
{
 
 public:
        double A[100][100];
        void Ingresar(int N);
        void Mostrar(int N);
        void Llenar_L(int N);
        void Llenar_U(int N);
        
}A,L,U;
 
class Vector
{
 
 public:
        double B[100];
        
        void Ingresar(int N)
        {
         for ( int i = 1 ; i <= N ; i++ )
             {
              cout<<"\t\t - ["<<i<<"] : " ;
              cin>>B[i];
             } 
        }
        
        void Mostrar(int N)
        {
          
          for( int i = 1 ; i <= N ; i++ )
              cout<<"\t\t º "<<setw(6)<<B[i]<<setw(5)<<"º"<<endl;;
        } 
}B,X,Y;


class Menu
{
 private :
         int opcion;
 public :
        int mostrar()
        {
         cout<<endl;
         cout<<"\t\t ÉÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ»\n";
         cout<<"\t\t º        Factorizacion LU               º\n";
         cout<<"\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n";
         cout<<"\t\t º     1. Factorizar en LU               º\n";
         cout<<"\t\t º     2. Calcular Determinate co LU     º\n"; 
         cout<<"\t\t º     3. Resolvemos ecuacion con LU     º\n";
         cout<<"\t\t º     4. Salir.                         º\n";
         cout<<"\t\t ÈÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍÍ¼\n";
 
         cout<<endl<<"\t Ingrese Opcion  "; 
         cout<<endl<<"\t ÄÄÄÄÄÄÄÄÄÄÄÄÄÄ : "; 
         cin>>opcion;
         return opcion;
       }
       
       
};

class L_U
{

 
 public:
        L_U Factorizar(int N)
        {
        double suma,suma1;
  
        L.Llenar_L(N);
        U.Llenar_U(N);
        
        for ( int i = 1 ; i <= N ; i++)
            {
            for ( int j = i ; j <= N ; j++)
                {
                suma = 0;
                for ( int k = 1 ; k <= i - 1 ; k++)
                     suma = L.A[i][k]*U.A[k][j] + suma;
                 
                 U.A[i][j] = A.A[i][j] - suma;
                }
            
            for ( int j = i + 1; j <= N ; j++)    
                {
                 suma = 0;     
                 for ( int k = 1 ; k <= i - 1 ; k++)
                      suma = L.A[j][k]*U.A[k][i] + suma;
             
                 L.A[j][i] = (A.A[j][i] - suma)/U.A[i][i];
                }
             L.A[i][i] = 1;
           }          


        }
        
        void mostrar_LU(int N)
        {
         cout<<endl<<endl<<"\t L = ";
         L.Mostrar(N);
                     
         cout<<endl<<endl<<"\t U = ";
         U.Mostrar(N);
        }
        
        
        L_U Ecuacion(int N)
        {
         Factorizar(N);
         mostrar_LU(N);
                  
         Vector X;
         
         double sum = 0;
         
         
         
         for( int i = 1 ; i <= N ; i++ )
            {
             sum = 0;
             
             for( int j = 1 ; j <= i - 1 ; j++)
                 sum = sum + L.A[i][j]*Y.B[j];
                 
             Y.B[i] = B.B[i]-sum/L.A[i][i];   
            }  
          
         cout<<endl<<endl<<"\t Y = "<<endl<<endl;      
         Y.Mostrar(N);
         
         X.B[N] = Y.B[N]/U.A[N][N];
         
         for( int i = N - 1 ; i >= 1 ; i-- )
            {
             sum = 0;
             
             for( int j = i + 1 ; j <= N ; j++)
                 sum = sum + U.A[i][j]*X.B[j];
             
                 
             X.B[i] = (Y.B[i]-sum)/U.A[i][i];  
            }
            
            cout<<endl<<endl<<"\t X = "<<endl<<endl;
            X.Mostrar(N);
            

        } 
        
        L_U Determinante(int N)
        {
         
         L_U Lu;
         Lu.Factorizar(N);
          
         double mult = 1;
         
         for(int i = 1 ; i <= N ; i++)
              mult = mult*U.A[i][i];
             
         cout<<endl<<"\t\t La Determinate Es : "<<mult<<endl;
         
        }
        
         
}LU;



void Matriz :: Ingresar(int N)
{
 cout<<endl;
 for ( int i = 1 ; i <= N ; i++ )
     for ( int j = 1 ; j <= N ; j++ )
         {
          cout<<"\t\t - ["<<i<<","<<j<<"] : " ;
          cin>>A[i][j];
         }  
}




















void Matriz :: Mostrar(int N)
{
 
 cout<<endl;
 for( int i = 1 ; i <= N ; i++ )
    {
     cout<<"\t\t º";
     for( int j = 1 ; j <= N ; j++ )
if((A[i][j]<0.00001&&A[i][j]>0)||(A[i][j]>-0.00001&&A[i][j]<0))
             cout<<setw(6)<<0<<setw(6);
          else 
             cout<<setw(6)<<setprecision(3)<<A[i][j]<<setw(6);
        cout<<"º"<<endl;
    }
    cout<<endl;
}


void Matriz :: Llenar_L(int N)
{
 
 for ( int i = 1 ; i <= N ; i++ )
      for ( int j = 1 ; j <= N ; j++ )
          {
           if ( i == j )
              A[i][j] = 1;
           else
              A[i][j] = 0;   
          }
}

void Matriz :: Llenar_U(int N)
{
 
 for ( int i = 1 ; i <= N ; i++ )
      for ( int j = 1 ; j <= N ; j++ )
          A[i][j] = 0;
}



int main()
{
 
 Menu op;
 Matriz L;

 int opcion;
 int N;
 char rpt = 's';
 
 bool salir=false;
 
  
 system("color 2f");
 
 while(rpt=='s'&&!salir)
      {
      opcion = op.mostrar();
 
      switch(opcion)
            {
             case 1 : 
                     cout<<"\n\n\t > Dimension de la Matriz : ";
                     cin>>N;
   
                     A.Ingresar(N);
                     cout<<endl<<endl<<"\t A = ";
                     A.Mostrar(N);

                     LU.Factorizar(N);
                     LU.mostrar_LU(N);
                     
                    break;
                    
             case 2 : 
                     cout<<"\n\n\t > Dimension de la Matriz : ";
                     cin>>N;
                     A.Ingresar(N);
                     cout<<endl<<endl<<"\t A = ";
                     A.Mostrar(N);
                    
                     LU.Determinante(N);

                                          
                    break;
                                        
             case 3 :
     
                     cout<<"\n\n\t > Dimension de la Matriz : ";
                     cin>>N;
                     
                     A.Ingresar(N);
                     cout<<endl<<endl<<"\t A = ";
                     A.Mostrar(N);
                    
                     cout<<"\n\n\t > Ingrese Vector : "<<endl<<endl;
                     B.Ingresar(N);
                     
                     LU.Ecuacion(N);
                     break;
                
             case 4 :
                     salir = true;
             }
             
       if(!salir)
         {
          cout<<"\n\n\t\t\t Desea continuar ?? : "; 
          cin>>rpt;
         }
         
         system("cls");      
                         
       }
       
 
 system("pause>null");
}
