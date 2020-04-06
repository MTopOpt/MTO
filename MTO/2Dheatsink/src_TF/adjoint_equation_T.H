
    for(i=0;i<2;i++)
    {
       #include "adjoint_heat_T.H" 
    }

    if(opt==1)
    { 
      for(i=0;i<200;i++)
      {
          #include "adjoint_flow_T.H"  
      }
    }
    else
    {
      for(i=0;i<60;i++)
      {
         #include "adjoint_flow_T.H"  
      }
    }
