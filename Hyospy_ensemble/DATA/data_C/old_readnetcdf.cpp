#include <iostream>
#include <netcdfcpp.h>
#include <string>  
#include <sstream>
#include <stdio.h>

using namespace std;


static const int NX=data_length;
static const int NY=2;
static const int NC_ERR=2;
static const int nsize=particle_num;
static const int row=row_num;
static const int col=column_num;
static const float west=west_point;
static const float east=east_point;
static const float south=south_point;
static const float north=north_point;

int main(){
    
    std::ostringstream fn;
    NcFile dataFile("GNOME_file", NcFile::ReadOnly);
    int nstep=NX/nsize;
    float dataIn1[NX];
    float dataIn2[NX];
    float location[NX][NY];
    float nloc[nstep][nsize][NY];
    float dy, dx;
    bool find;

    if (!dataFile.is_valid()){
        cout << "Couldn't open file!\n";
        return NC_ERR;
        }
        
    NcVar *data1=dataFile.get_var("latitude");
    data1->get(&dataIn1[0], NX);

    NcVar *data2=dataFile.get_var("longitude");
    data2->get(&dataIn2[0], NX);
    
    cout.precision(15);
        
    for (int i=0; i<=NX-1; i++){
            location[i][0]=dataIn1[i];
            location[i][1]=dataIn2[i];
    }
    
    //**regroup the locations, every 200 points at one time step is in one group
    for (int i=0; i<=nstep-1; i++){
         for (int j=0; j<=nsize-1; j++){
             nloc[i][j][0]=location[i*nsize+j][0];
             nloc[i][j][1]=location[i*nsize+j][1];
         }
    }
    
    /*
    step 2 -- find out the 4 vertex points for the 'map'
    */    
    //west=-95.78; east=-94.45;
    //south=28.3;  north=29.9; 
    /*
    step 3 -- subdivide the rectangle into small cells (row*col) and initiate count to be 0
    */
    dy=(north-south)/row; 
    dx=(east-west)/col;
    int count[row][col];
      for (int i=0; i<=row-1; i++){
        for (int j = 0; j <= col-1; j++)
            count[i][j]=0;
      }
    
    /*
    step 4 -- write the output into netcdf file
    */  
    
    int accumulator=0;
    for (int i=0; i<=nstep-1; i++){
      for (int j=0; j<=nsize-1; j++){
            for (int p = 0; p <= row-1; p++){
                find= false;
                for (int k = 0; k <= col-1; k++){
                    if ((south+dy*p) <= nloc[i][j][0] && nloc[i][j][0] <= (south+dy*(p+1)) && (west+dx*k) <=nloc[i][j][1] &&  nloc[i][j][1]<=(west+dx*(k+1))){
                       count[p][k]+=1;  
                       find=true;
                       accumulator++;
                       break;}       
                }  
                if (find==true)
                    break;
             }
       }
       
       char fn[32];
       snprintf(fn, sizeof(char) * 32, "out_count%i.nc", i);
       NcFile dataFile2(fn, NcFile::Replace); 
       if (!dataFile2.is_valid())
           {
              cout << "Couldn't open file!\n";
              return NC_ERR;
           }

       NcDim* xDim = dataFile2.add_dim("x", row);
       NcDim* yDim = dataFile2.add_dim("y", col); 
 
       NcVar *data = dataFile2.add_var("data", ncInt, xDim, yDim);
       data->put(&count[0][0], row, col);    
    }
    cout << "*** SUCCESS writing example file out_count.nc!" << endl; 
       
    return 0;
}
