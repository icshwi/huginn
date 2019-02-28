program huginn


%%#include <stdio.h>
%%#include <stdlib.h>
%%#include <errlog.h>

int	huginnMode;
assign	huginnMode 	to "{P}:Mode";
monitor	huginnMode;

int	currentMode;		
int	currentIsSetPointOK;	
int	currentSubSetPoint;	
int	currentMainSetPoint;	
int	currentMainOn;		
int	currentSubOn;		

float	subSetPoint;
assign	subSetPoint	to "{P}:Sub_SETP";
monitor	subSetPoint;

float	mainSetPoint;
assign	mainSetPoint	to "{P}:MainCryo_SETP";
monitor	mainSetPoint;

int	subOn;
assign	subOn		to "{P}:Sub_ON";
monitor	subOn;

int	mainOn;
assign	mainOn		to "{P}:MainON";
monitor	mainOn;

int	subEnable;
assign	subEnable	to "{P}:SubEnable";
monitor	subEnable;

int	mainEnable;
assign	mainEnable	to "{P}:MainEnable";
monitor	mainEnable;

int	subSetPointEnable;
assign	subSetPointEnable to "{P}:SubSetPointEnable";
monitor	subSetPointEnable;

int	mainSetPointEnable;
assign	mainSetPointEnable to "{P}:SubSetPointEnable"; 
monitor	mainSetPointEnable;

int	isSetPointOK;
assign	isSetPointOK	to "{P}:IsSetPointOK";
monitor isSetPointOK;

string	message;
assign	message		to "{P}:Message";

string	warning;
assign	warning		to "{P}:Warning";

//option +r;

ss MODE {
//--------------------------------------
   state init {

      when(){
       printf("MESSAGE: init\n");
       strcpy(message,"Start Huginn");
       pvPut(message, SYNC);
      } state mode

   }

//--------------------------------------
   state mode {

      when(huginnMode == 0) {
          printf("--- mode, 0\n");
      } state manual

      when(huginnMode == 1) {
          printf("--- mode, 1\n");
      } state autoMainExtended

      when(huginnMode == 2) {
          printf("--- mode, 2\n");
      } state autoSub


   }
 
//--------------------------------------
    state waitForChangeState{
      
      entry{
            }  

      when(currentMode != huginnMode) {
          printf("wait for changing state, Mode\n");
      } state mode

      when(currentIsSetPointOK != isSetPointOK) {
          printf("wait for changing state, isSetPointOK\n");
      } state mode

      when(currentSubSetPoint != subSetPoint)  {
          printf("wait for changing state, subSetPoint\n");
      } state mode

      when(currentMainSetPoint != mainSetPoint) {
          printf("wait for changing state, mainSetPoint\n");
      } state mode

      when(currentMainOn != mainOn) {
          printf("wait for changing state, main ON/OFF\n");
      } state mode

      when(currentSubOn != subOn) {
          printf("wait for changing state, sub ON/OFF\n");
      } state mode

    }

//--------------------------------------
   state manual{

      entry{ 
	currentMode 		= huginnMode;
        currentSubSetPoint	= subSetPoint;
        currentIsSetPointOK	= isSetPointOK;
        currentMainSetPoint	= mainSetPoint;
        currentMainOn		= mainOn;
        currentSubOn		= subOn; 

        strcpy(message,"Init Manual Mode");
        pvPut(message, SYNC);

//	subOn = 0;
//        pvPut(subOn, SYNC);

//        mainOn = 0;
//        pvPut(mainOn, SYNC);

	mainEnable = 1;
	pvPut(mainEnable);
      }
    
      when(isSetPointOK){
	subEnable = 1;
	pvPut(subEnable);

        strcpy(warning,"");
        pvPut(warning, SYNC);
        printf("MESSAGE Manual Mode, isSetPointOk = 1\n");
     } state waitForChangeState

      when(!isSetPointOK){
	subEnable = 0;
        pvPut(subEnable, SYNC);
 
        if(subOn == 1){
	  subOn = 0;
          pvPut(subOn);
        }

        strcpy(warning,"SetPoint is outside the limits!");
        pvPut(warning, SYNC);
        printf("MESSAGE Manual Mode, isSetPointOk = 0\n");
        } state waitForChangeState

   }

//--------------------------------------
   state autoMainExtended {
      entry{ 
        currentMode = huginnMode; 
        currentSubSetPoint	= subSetPoint;
        currentIsSetPointOK	= isSetPointOK;
        currentMainSetPoint	= mainSetPoint;
        currentMainOn		= mainOn;
        currentSubOn		= subOn; 

        strcpy(message,"Init Auto-Extended Mode");
        pvPut(message, SYNC);

	subEnable = 0;
	pvPut(subEnable, SYNC);

	mainEnable = 1;
	pvPut(mainEnable, SYNC);

        subSetPointEnable = 1;
        pvPut(subSetPointEnable, SYNC);

        mainSetPointEnable = 1;
        pvPut(mainSetPointEnable, SYNC);
      }

      when(isSetPointOK){
          //subOn = 1;
          //pvPut(subOn, SYNC);
          strcpy(warning,"");
          pvPut(warning, SYNC);
          printf("MESSAGE Auto-Extend Mode, isSetPointOk = 1\n");
      } state waitForChangeState

      when(!isSetPointOK){
          //subOn = 0;
          //pvPut(subOn, SYNC);
          strcpy(warning,"SetPoint is outside the limits!");
          pvPut(warning, SYNC);
          printf("MESSAGE Auto-Extend Mode, isSetPointOk = 0\n");
      } state waitForChangeState

//        printf("MESSAGE: Mode Auto-Extend\n");
   
   }

//--------------------------------------
   state autoSub{
      entry{ 
	currentMode = huginnMode;
        currentSubSetPoint	= subSetPoint;
        currentIsSetPointOK	= isSetPointOK;
        currentMainSetPoint	= mainSetPoint;
        currentMainOn		= mainOn;
        currentSubOn		= subOn; 
        strcpy(message,"Init Auto Mode");
        pvPut(message, SYNC);

	subEnable = 0;
	pvPut(subEnable, SYNC);

	mainEnable = 0;
	pvPut(mainEnable, SYNC);

        mainSetPoint = subSetPoint;
        pvPut(mainSetPoint, SYNC); 
        mainSetPointEnable = 0;
        pvPut(mainSetPointEnable);

        subSetPointEnable = 1;
        pvPut(subSetPointEnable, SYNC);

        printf("MESSAGE Auto Mode\n");
      }

      when(isSetPointOK){
          //subOn = 1;
          //pvPut(subOn, SYNC);
          strcpy(warning,"");
          pvPut(warning, SYNC);
          printf("MESSAGE Auto Mode, isSetPointOk = 1\n");
      } state waitForChangeState

      when(!isSetPointOK){
          //subOn = 0;
          //pvPut(subOn, SYNC);
          strcpy(warning,"SetPoint is outside the limits!");
          pvPut(warning, SYNC);
          printf("MESSAGE Auto Mode, isSetPointOk = 0\n");
      } state waitForChangeState
   }
}

