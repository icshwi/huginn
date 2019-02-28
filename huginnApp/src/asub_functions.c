#include <aSubRecord.h>
#include <registryFunction.h>
#include <epicsExport.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

 
int getparam1(char *p, char *patern, float *a, int nr);
//===============================================================================================================
// This function implement linear interpolation 
// For giving arrays values of CalibTempa (aCT) and tbase 
// the fuction interpolate a value from arrays DeltaUp (aDU) amd DeltaDown (aDD)
//
// aCT MUST BE MONOTONIC from max value at the first element to min value at the last element
// as it is done in .sav file, tabY does not have to be monotonic.
// Fuction asumes that number of elements in both tab are the same.
 static long asub_limits(aSubRecord *precord) {

    float *aCT, *aDU, *aDD;
    unsigned long nCT, nDU, nDD;
    float tbase;
    float setPoint;
    int delayOff;
    int counter;
    float safetyFactor;

    float dx, dy;
    float up, down;
    int isSetPointOK;

    int i; 
 
    aCT 	 = (float *)precord->a;		// CalibTemp
    nCT 	 = precord->nea;   	 	// nr of elements in CalibTemp
     
    aDU 	 = (float *)precord->b;	 	// DeltaUp
    nDU 	 = precord->neb;   		// nr of elements in DeltaUp
 
    aDD 	 = (float *)precord->c;		// DeltaDown
    nDD		 = precord->nec;		// nr of elements in DeltaDown
 
    tbase	 = *(float*)precord->d; 	// tbase from subcriostat 
    setPoint	 = *(float*)precord->e; 	// subCryostat setpoint 
    safetyFactor = *(float*)precord->f; 	// safety factor
    delayOff 	 = *(short*)precord->g; 	// delay to switch off sub-cryo 
    counter	 = *(short*)precord->h; 	// counter to count seconds

//////////////////////////////////////////////////////////////////////////////////////
/*  
    printf("\n------ [DEBUG]:nCT=%ld ", nCT);
    for (i = 0; i < nCT; i++) {
       printf("%f ", aCT[i]);
    }
    printf("\n------ [DEBUG]:nDU=%ld ", nDU);
    for (i = 0; i < nDU; i++) {
       printf("%f ", aDU[i]);
    }
    printf("\n------ [DEBUG]:nDD=%ld ", nDD);
    for (i = 0; i < nDD; i++) {
       printf("%f ", aDD[i]);
    }
    printf("\n------ [DEBUG]:tbase=%f\n", tbase);
*/
//////////////////////////////////////////////////////////////////////////////////////
// first for DU

//    printf("------ [DEBUG]:tbase = %f\n", tbase);
    if (tbase > aCT[0]) { 		
       up = aDU[0];
//       printf("------ [DEBUG]:up = %f\n", aDU[0]);
    }
 
    else if (tbase < aCT[nCT-1]) {	
       up = aDU[nDU-1];
//       printf("------ [DEBUG]:up =  %f\n", aDU[nDU-1]);
    }
    else {                     		
       for (i = 0; i < nCT-1; i++) {
           if (aCT[i] > tbase) {
               break;
           }
       }
       dx = aCT[i+1] - aCT[i];
       dy = aDU[i+1] - aDU[i];
       up = aDU[i] + (tbase - aCT[i]) * dy / dx; 
//       printf("------ [DEBUG]: up interpolate = %f \n", up );
   }

////////////////////////////////////////////////////////////////////////////////////////////////
// then for DD
    if (tbase > aCT[0]) {
       down = aDD[0];
//       printf("------ [DEBUG]:down = %f\n", aDD[0]);
    }
 
    else if (tbase < aCT[nCT-1]) {
       down = aDD[nDD-1];
//       printf("------ [DEBUG]:down =  %f\n", aDD[nDD-1]);
    }
    else {                     		
       for (i = 0; i < nCT-1; i++) {
           if (aCT[i] > tbase) {
               break;
           }
       }
       dx = aCT[i+1] - aCT[i];
       dy = aDD[i+1] - aDD[i];
       down = aDD[i] + (tbase - aCT[i]) * dy / dx; 
//       printf("------ [DEBUG]: down interpolate = %f \n",down );
   }



   *(float *)precord->vala = tbase + up;
   *(float *)precord->valb = tbase - down;
//   printf("limits are up:%f, down:%f \n", tbase + up, tbase - down);  
    
    //printf("------ [DEBUG]:\tsetPoint = %f\n", setPoint );
    if( (setPoint > (tbase - down*safetyFactor)) && (setPoint < (tbase + up*safetyFactor))  )
      isSetPointOK = 1;
    else
      isSetPointOK = 0;
    //printf("------ [DEBUG]:\tisSetPointOK = %d\n", isSetPointOK );


    if(isSetPointOK == 1){
       *(short *)precord->vald = 1;
       *(short *)precord->valc = counter = delayOff;
       //printf("------ [DEBUG]:(1)\tisSetPointOK = %d, counter = %d\n", isSetPointOK, counter );
    }
    else if (isSetPointOK == 0 && counter > 0){
       counter--;
       *(short *)precord->vald = 1;
       *(short *)precord->valc = counter;
       //printf("------ [DEBUG]:(2)\tisSetPointOK = %d, counter = %d\n", isSetPointOK, counter );
    }
    else {  
       *(short *)precord->vald = 0;
       *(short *)precord->valc = 0;
       //printf("------ [DEBUG]:(3)\tisSetPointOK = %d, counter = %d\n", isSetPointOK, counter );
    }


return 0;
}
 
//===============================================================================================================

static long asub_inrange(aSubRecord *precord) {

    float tp1		= *(float*)precord->a; 
    float tp2		= *(float*)precord->b; 
    float setpoint	= *(float*)precord->c; 
    float tolerance	= *(float*)precord->d; 
    int   toleranceTime	= *(short*)precord->e;
    int   counter	= *(short*)precord->f;
 
    int inrange;

    if( (abs(tp1 - setpoint) < tolerance ) && (abs(tp2 - setpoint) < tolerance) )
        inrange = 1;
    else 
        inrange = 0;


//   printf("------ [DEBUG]: inrange = %d\n", inrange );

    if(inrange && counter > 0){
       counter--;
       *(short *)precord->vala = 0;
       *(short *)precord->valb = counter;
//       printf("------ [DEBUG]:\tstart counting, count = %d\n", count );
    }
    else if (inrange && counter <= 0){
       *(short *)precord->vala = 1;
       *(short *)precord->valb = 0;
//       printf("------ [DEBUG]:\tready to measure\n" );
    }
    else {  
       *(short *)precord->vala = 0;
       *(short *)precord->valb = toleranceTime;
//       printf("------ [DEBUG]:\thug:Counts = %d\n", toleranceTime );
    }

return 0;
}

//===============================================================================================================
static long asub_init(aSubRecord *precord) {

    printf("=================================\n");
    printf("===    asub, initialization     =\n");



    char *s = (char*)precord->a;
    printf("=== precord->a = %s\n", s);
 
    int chn = 0, a=0, b=0, c=0;
    sscanf(s,"%d , %d , %d , %d",&chn, &a, &b, &c); 

    printf("=== aSub conversion: a=%d, b=%d, c=%d, s=%s\n",a,b,c,s);
    *(short *)precord->vala = (short)a;
    *(short *)precord->valb = (short)b;
    *(short *)precord->valc = (short)c;

    printf("=================================\n");
return 0;
}

//===============================================================================================================

static long asub_test(aSubRecord *p) {
return 0;
}


//===============================================================================================================
int getparam1(char *p, char *patern, float *a, int nr){
char *p_const = p;
char *c;

        while(isspace(*p))
                p++;

        if(*p == '\0')
                return 0;
        if(*p == '#')
                return 0;

        if( (c = strstr(p, patern)) ){
           sscanf(p_const,"%*s %f , %f , %f , %f",&a[0], &a[1], &a[2], &a[3]);
           return 4;
        }


return 0;
}

 
epicsRegisterFunction(asub_init);
epicsRegisterFunction(asub_test);
epicsRegisterFunction(asub_inrange);
epicsRegisterFunction(asub_limits);
