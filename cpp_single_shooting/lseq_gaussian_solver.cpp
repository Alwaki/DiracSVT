#include<math.h>
#include <vector>

std::vector<double> gaussianElimination(std::vector<std::vector<double>> mat)
{
    int i,j,k;
    const int n = mat.size();
    std::vector<double> res(n, 0);
    for(i=0;i<n;i++) 
    {                   
        for(j=i+1;j<n;j++)
        {
            if(abs(mat[i][i]) < abs(mat[j][i]))
            {
                for(k=0;k<n+1;k++)
                {
                    /* swapping mat[i][k] and mat[j][k] */
                    mat[i][k]=mat[i][k]+mat[j][k];
                    mat[j][k]=mat[i][k]-mat[j][k];
                    mat[i][k]=mat[i][k]-mat[j][k];
                }
            }
      }
    }
   
     /* performing Gaussian elimination */
    for(i=0;i<n-1;i++)
    {
        for(j=i+1;j<n;j++)
        {
            float f=mat[j][i]/mat[i][i];
            for(k=0;k<n+1;k++)
            {
              mat[j][k]=mat[j][k]-f*mat[i][k];
      }
        }
    }
    /* Backward substitution for discovering values of unknowns */
    for(i=n-1;i>=0;i--)          
    {                     
        res[i]=mat[i][n];
                    
        for(j=i+1;j<n;j++)
        {
          if(i!=j)
          {
              res[i]=res[i]-mat[i][j]*res[j];
    }          
  }
  res[i]=res[i]/mat[i][i];  
    }
return res;
}