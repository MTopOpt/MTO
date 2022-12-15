//*************flow problem******************************//
   if(opt==1)
   { 
     for(i=0;i<200;i++)
     {
        #include "primal_flow.H"
     }
   }
   else
   {
     for(i=0;i<40;i++)
     {
        #include "primal_flow.H"
     }
   }
//*************heat transfer*******************************//
  for(i=0;i<2;i++)
  {
     #include "primal_heat.H"  
  }

