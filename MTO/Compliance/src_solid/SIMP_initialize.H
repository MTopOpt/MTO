// 
    std::remove("Voluse.txt");
    std::remove("Time.txt");
    std::remove("Compliance.txt"); 

    scalar C(0);
    scalar V(0);
    scalar scalef(0);
    scalar eps(1e-8);
    int i=0,j=0,k=0,jk,allcells,nallcells,mycells=0,opt=1;
    scalar allcell=gamma.size();//allcell:length of gamma in the i-th core
    allcells=static_cast<int>(allcell); 
    List<scalar> cellsize(Pstream::nProcs());
    for(i=0;i<Pstream::nProcs();i++)
    {
       cellsize[i]=0;
    }
    cellsize[Pstream::myProcNo()]=allcells;
    reduce(cellsize, sumOp<List<scalar>>()); 
    reduce(allcell, sumOp<scalar>()); 
    nallcells=static_cast<int>(allcell);//nallcells:total length of gamma
    for(i=0;i<Pstream::myProcNo();i++)
    {
       mycells=mycells+cellsize[i];
    }     

    double add1=0,add2=0,len=std::pow(mesh.V()[0],1.0/3);
    double eta0,eta1,eta5=0,del=0.1,y0,y1,y5; 

    double *Faa =new double[nallcells];
    double *Gbb =new double[nallcells];
    double *Gcc =new double[nallcells];
    double *Gdd =new double[nallcells];
    double *gsens =new double[nallcells];
    double *gamma_densityfilter=new double[nallcells];
    double *drho=new double[nallcells];

 // Heaviside density filter (only works on uniform mesh)
    List<vectorField> AllcellList(Pstream::nProcs());
    vectorField localcell(mesh.C());
    AllcellList[Pstream::myProcNo()] = localcell;
    Pstream::gatherList(AllcellList);
    Pstream::scatterList(AllcellList);
    vectorField Allcellcenter 
    (
        ListListOps::combine<vectorField>
        (
            AllcellList,
            accessOp<vectorField>()
        )
    );  
    int FR=filter_Radius;
    int Nnum[nallcells];
    int (*Nlabel)[125]=new int [nallcells][125];
    double (*Nweight)[125]=new double [nallcells][125];
    double geo_xmax=-1000,geo_ymax=-1000,geo_zmax=-1000,geo_xmin=1000,geo_ymin=1000,geo_zmin=1000;
    for(i=0;i<nallcells;i++)
    {
       geo_xmax=std::max(geo_xmax,Allcellcenter[i][0]);
       geo_ymax=std::max(geo_ymax,Allcellcenter[i][1]);
       geo_zmax=std::max(geo_zmax,Allcellcenter[i][2]);
       geo_xmin=std::min(geo_xmin,Allcellcenter[i][0]);
       geo_ymin=std::min(geo_ymin,Allcellcenter[i][1]);
       geo_zmin=std::min(geo_zmin,Allcellcenter[i][2]);
    }
    int nx,ny,nz,Nx,Ny,Nz;
    nx=static_cast<int>((geo_xmax-geo_xmin+len/4.0)/len)+FR*2+1;
    ny=static_cast<int>((geo_ymax-geo_ymin+len/4.0)/len)+FR*2+1;
    nz=static_cast<int>((geo_zmax-geo_zmin+len/4.0)/len)+FR*2+1;

    int ***ALLlabel;
    ALLlabel = new int **[nx]; 
    for(i=0; i<nx; i++) 
    { 
       ALLlabel[i]=new int *[ny]; 
       for(j=0; j<ny; j++)
       { 
          ALLlabel[i][j]=new int[nz]; 
       }
    }
    for(i=0;i<nx;i++)
    {
       for(j=0;j<ny;j++)
       {
          for(k=0;k<nz;k++)
          {
             ALLlabel[i][j][k]=-1;
          }
       }
    }
   for(jk=0;jk<nallcells;jk++)
   {
      nx=static_cast<int>((Allcellcenter[jk][0]-geo_xmin+len/4.0)/len);
      ny=static_cast<int>((Allcellcenter[jk][1]-geo_ymin+len/4.0)/len);
      nz=static_cast<int>((Allcellcenter[jk][2]-geo_zmin+len/4.0)/len);
      ALLlabel[nx+FR][ny+FR][nz+FR]=jk;
   }

   for(jk=0;jk<nallcells;jk++)
   {
      Nnum[jk]=0;
      nx=static_cast<int>((Allcellcenter[jk][0]-geo_xmin+len/4.0)/len);
      ny=static_cast<int>((Allcellcenter[jk][1]-geo_ymin+len/4.0)/len);
      nz=static_cast<int>((Allcellcenter[jk][2]-geo_zmin+len/4.0)/len);
      for(i=-FR;i<=FR;i++)
      {
         for(j=-FR;j<=FR;j++)
         {
            for(k=-FR;k<=FR;k++)
            {
               Nx=nx+i;
               Ny=ny+j;
               Nz=nz+k; 
               if(ALLlabel[Nx+FR][Ny+FR][Nz+FR]>=0)
               {
                  Nlabel[jk][Nnum[jk]]=ALLlabel[Nx+FR][Ny+FR][Nz+FR];
                  Nweight[jk][Nnum[jk]]=(3.8-std::sqrt(i*i+j*j+k*k))*len;
                  Nnum[jk]++;
               }
            }
         }
      }
   }
//***********************************//Petsc
PetscInt n_dof,N_dof,m=1;//m:number of constraints
n_dof=allcells;
N_dof=nallcells;

PetscScalar Xmin=0; // Min. value of design variables
PetscScalar Xmax=1; // Max. value of design variables

Vec xmma,xold,xmin,xmax;
PetscScalar *xmma_array;
VecCreate(PETSC_COMM_WORLD,&xmma);
VecSetSizes(xmma,n_dof,N_dof);
VecSetFromOptions(xmma);
VecGetArray(xmma,&xmma_array);
for(i=0;i<allcells;i++)
{
  xmma_array[i]=gamma[i];
}
VecRestoreArray(xmma,&xmma_array);

VecDuplicate(xmma,&xold); 
VecDuplicate(xmma,&xmin);
VecDuplicate(xmma,&xmax); 

PetscScalar aMMA[m];
PetscScalar cMMA[m];
PetscScalar dMMA[m];

for (PetscInt i=0;i<m;i++)
{
  aMMA[i]=0.0;
  dMMA[i]=1.0;
  cMMA[i]=1000.0;
}

MMA *mma;
mma=new MMA(N_dof,m,xmma, aMMA, cMMA, dMMA,mma_init,mma_dec,mma_inc);

//**********************************************//dfdx
Vec dfdx; // Sensitivities of objective
VecDuplicate(xmma,&dfdx); 
PetscScalar *dfdx_array;
VecGetArray(dfdx,&dfdx_array);
VecRestoreArray(dfdx,&dfdx_array); 
//*****************************************************//gx
PetscScalar *gx= new PetscScalar[m]; // Array with constraint values

//*****************************************************//dgdx
Vec *dgdx; // Sensitivities of constraints (vector array)
VecDuplicateVecs(xmma,m, &dgdx); 
PetscScalar **dgdx_array;
VecGetArrays(dgdx,m,&dgdx_array);
VecRestoreArrays(dgdx,m,&dgdx_array);

gamma.correctBoundaryConditions();
mu=(E*eps+gamma*gamma*gamma*(1-eps)*E)/(2.0*(1.0 + Po));
lambda=(E*eps+gamma*gamma*gamma*(1-eps)*E)*Po/((1.0 + Po)*(1.0 - 2.0*Po));
if (planeStress)
{
   lambda = (E*eps+gamma*gamma*gamma*(1-eps)*E)*Po/((1.0 + Po)*(1.0 - Po));
}
