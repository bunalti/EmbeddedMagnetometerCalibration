/****************************************************************
 * @brief            ICM 20948 Embedded Magnetometer Calibration 
 * 
 * 
 * 
 * Initialize elliptical fit using matrix inverse and eigenvalue
 * decomposition. 
 * 
 * 
 * 
 * 
 * 
 * 
 * @functions        MagCal_Run(void):  Checks if enough data is  
 *                   collected and applies possible elliptic fit
 *                   methods accounting the number of data.
 *                   (Check magcal.c lines 44-46)
 *                   
 *                   
 *                   Mag_Collect(icm_20948_DMP_data_t data):
 *                   Check if valid data is available on
 *                   FIFO of ICM bus. Writes raw magnetometer
 *                   data to appropiate calibration samples
 *                   buffers
 *                
 *                   
 *                   
 * 
 * 
 * Kerem Tecirlioğlu, April 25th, 2021
 * Distributed as-is; no warranty is given.
 ***************************************************************/

#include <ICM_20948.h>
#include <magcal.h>


#define SERIAL_PORT Serial

#define WIRE_PORT Wire // Your desired Wire port.   
#define AD0_VAL 0      // The value of the last bit of the I2C address.                

ICM_20948_I2C myICM;          // Create an ICM_20948_I2C object
icm_20948_DMP_data_t data;    // Create DMP data storage object

Point_t accelerometer;
Point_t gyroscope;
Point_t compass;


double roll;
double pitch;
double yaw;
double roll_deg;
double pitch_deg;
double yaw_deg;




int Mag_Collect(icm_20948_DMP_data_t data,uint16_t samples){

  uint16_t dex = 0;

  
  while(true){
  
    // Read any DMP data waiting in the FIFO
    myICM.readDMPdataFromFIFO(&data);
    
   
    if ((myICM.status == ICM_20948_Stat_Ok) || (myICM.status == ICM_20948_Stat_FIFOMoreDataAvail)) // Was valid data available?
    {
  
      if ((data.header & DMP_header_bitmap_Compass) > 0) // Check for Compass
      {

        magcal.BpFast[0][dex] = (float)data.Compass.Data.X; // Extract the compass data
        magcal.BpFast[1][dex] = (float)data.Compass.Data.Y;
        magcal.BpFast[2][dex] = (float)data.Compass.Data.Z;
        magcal.valid[dex] = 1;
        dex++;
  
        
        if(dex > samples){
          
          dex = 0;
          Serial.println(F("Data collection complete"));
          
          MagCal_Run();


          SERIAL_PORT.println(F("------------------------------------"));
          SERIAL_PORT.print(F("Hard Iron Offset:      X:"));
          SERIAL_PORT.print(magcal.trV[0]);
          SERIAL_PORT.print(F(" Y:"));
          SERIAL_PORT.print(magcal.trV[1]);
          SERIAL_PORT.print(F(" Z:"));
          SERIAL_PORT.println(magcal.trV[2]);
          SERIAL_PORT.println("");
    
          SERIAL_PORT.print(F("Soft Iron Offset:      XX:"));
          SERIAL_PORT.print(magcal.trinvW[0][0]);
          SERIAL_PORT.print(F(" XY:"));
          SERIAL_PORT.print(magcal.trinvW[0][1]);
          SERIAL_PORT.print(F(" XZ:"));
          SERIAL_PORT.println(magcal.trinvW[0][2]);
    
    
          SERIAL_PORT.print(F("Soft Iron Offset:      YX:"));
          SERIAL_PORT.print(magcal.trinvW[1][0]);
          SERIAL_PORT.print(F(" YY:"));
          SERIAL_PORT.print(magcal.trinvW[1][1]);
          SERIAL_PORT.print(F(" YZ:"));
          SERIAL_PORT.println(magcal.trinvW[1][2]);
    
          SERIAL_PORT.print(F("Soft Iron Offset:      ZX:"));
          SERIAL_PORT.print(magcal.trinvW[2][0]);
          SERIAL_PORT.print(F(" ZY:"));
          SERIAL_PORT.print(magcal.trinvW[2][1]);
          SERIAL_PORT.print(F(" ZZ:"));
          SERIAL_PORT.println(magcal.trinvW[2][2]);

          SERIAL_PORT.println(F("------------------------------------"));   
        while(1);
        
        }
    
      }
    
    } 
  
  }
  
  
  
  
  
  
}


void setup()
{

  SERIAL_PORT.begin(115200); // Start the serial console
  SERIAL_PORT.println(F("ICM 20948 Embedded Magnetometer Calibration"));
  
  delay(100);

  while (SERIAL_PORT.available()) // Make sure the serial RX buffer is empty
    SERIAL_PORT.read();


  
  WIRE_PORT.begin(21,22);
  WIRE_PORT.setClock(400000);

  bool initialized = false;
  while (!initialized)
  {

    // Initialize the ICM-20948
    // If the DMP is enabled, .begin performs a minimal startup. We need to configure the sample mode etc. manually.
    myICM.begin(WIRE_PORT, AD0_VAL);


    SERIAL_PORT.print(F("Initialization of the sensor returned: "));
    SERIAL_PORT.println(myICM.statusString());
    if (myICM.status != ICM_20948_Stat_Ok)
    {
      SERIAL_PORT.println(F("Trying again..."));
      delay(500);
    }
    else
    {
      initialized = true;
    }
  }

  SERIAL_PORT.println(F("Device connected!"));

  bool success = true; // Use success to show if the DMP configuration was successful

  // Initialize the DMP. initializeDMP is a weak function. You can overwrite it if you want to e.g. to change the sample rate
  success &= (myICM.initializeDMP() == ICM_20948_Stat_Ok);

  // Enable sensors / features
  success &= (myICM.enableDMPSensor(INV_ICM20948_SENSOR_MAGNETIC_FIELD_UNCALIBRATED) == ICM_20948_Stat_Ok);


  // Configuring DMP to output data at multiple ODRs:
  // DMP is capable of outputting multiple sensor data at different rates to FIFO.
  // Setting value can be calculated as follows:
  // Value = (DMP running rate / ODR ) - 1
  success &= (myICM.setDMPODRrate(DMP_ODR_Reg_Cpass, 0) == ICM_20948_Stat_Ok);        //


  // Enable the FIFO
  success &= (myICM.enableFIFO() == ICM_20948_Stat_Ok);

  // Enable the DMP
  success &= (myICM.enableDMP() == ICM_20948_Stat_Ok);

  // Reset DMP
  success &= (myICM.resetDMP() == ICM_20948_Stat_Ok);

  // Reset FIFO
  success &= (myICM.resetFIFO() == ICM_20948_Stat_Ok);

  // Check success
  if (success)
  {
    SERIAL_PORT.println(F("DMP enabled!"));
  }
  else
  {
    SERIAL_PORT.println(F("Enable DMP failed!"));
    SERIAL_PORT.println(F("Please check that you have uncommented line 29 (#define ICM_20948_USE_DMP) in ICM_20948_C.h..."));
    while (1)
      ; // Do nothing more
  }
}

void loop()
{

  Mag_Collect(data,400);

  if (myICM.status != ICM_20948_Stat_FIFOMoreDataAvail) // If more data is available then we should read it right away - and not delay
  {
    delay(10);
  }
}
