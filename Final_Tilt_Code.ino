#include <Arduino_LSM9DS1.h>
#include <MadgwickAHRS.h>
#include <Wire.h>
// Madgwick
Madgwick filter;

int ADXL345 = 0x53; // The ADXL345 sensor I2C address

// sensor's sample rate is fixed at 119 Hz:
const float sensorRate = 119;

// Measured inclination based on accelerometer X axis
float thetaM;
// Measured inclination based on accelerometer y axis
float phiM;
// Measured inclination based on gyroscope x axis
float thetaG = 0;
// Measured inclination based on gyroscope y axis
float phiG = 0;

// Final tilt angle values
float roll;
float pitch;
float yaw;

// Time Management
float dt;
unsigned long millisOld;
int timeInterval = 100;

void setup() {
  IMU.begin();
  Serial.begin(9600);

  // Caibration of Magnetometer - Replaceable
  IMU.setMagnetFS(0);
  IMU.setMagnetODR(8);
  IMU.setMagnetOffset(19.246216, 18.193359, 3.260498);
  IMU.setMagnetSlope (1.859274, 1.662833, 1.755710);
  // End of Calibration Values
  IMU.magnetUnit = MICROTESLA;

  Wire.begin(); // initiate Wire library
  Wire.beginTransmission(ADXL345); // Start communicating with the device
  Wire.write(0x2D); // Access/ talk to POWER_CTL Register - 0x2D
  // Enable measurement
  Wire.write(8); // Bit D3 High for measuring enable (8dec -> 0000 1000 binary)
  Wire.endTransmission();
  delay(10);

  //Off-set Calibration
  //X-axis
  Wire.beginTransmission(ADXL345);
  Wire.write(0x1E);
  Wire.write(1);
  Wire.endTransmission();
  delay(10);
  //Y-axis
  Wire.beginTransmission(ADXL345);
  Wire.write(0x1F);
  Wire.write(-2);
  Wire.endTransmission();
  delay(10);

  //Z-axis
  Wire.beginTransmission(ADXL345);
  Wire.write(0x20);
  Wire.write(-9);
  Wire.endTransmission();
  delay(10);

  filter.begin(119);
  delay(1000);
  millisOld = millis();
}

void loop() {

  Wire.beginTransmission(ADXL345);
  Wire.write(0x32); // Start with register 0x32 (ACCEL_XOUT_H)
  Wire.endTransmission(false);

  float accelx, accely, accelz;
  float gyrox, gyroy, gyroz;
  float magx, magy, magz;

  // Read Accelerometer Data
  IMU.readAcceleration(accelx, accely, accelz);
  // Read Gyroscope Data
  IMU.readGyroscope(gyrox, gyroy, gyroz);
  // Read Magnetomter Data
  IMU.readMagneticField(magx, magy, magz);

  // Calculating the inclination based on accelerometer data x axis
  thetaM = atan2(accelx, accelz) * (180 / PI);
  // Calculating the inclination based on accelerometer data y axis
  phiM = atan2(accely, sqrt(pow(accelx, 2) + pow(accelz, 2))) * (180 / PI);

  dt = (millis() - millisOld) / 1000;
  // Resetting the time with current time
  millisOld = millis();

  // Calculate the new Theta angle based on gyroscope
  thetaG = thetaG + gyroy * dt;
  // Calculate the new Phi angle based on gyroscope
  phiG = phiG - gyrox * dt;

  // Apply complementary filter, gyro high pass, accel low pass to account for gyro drift long term
  // Estimated values with low lag and high resistance to noise due to acceleration
  pitch = ((pitch + gyroy * dt) * .94) + (thetaM * .06);
  roll = ((roll - gyrox * dt) * .94) + (phiM * .06);

  /* Not that great
    // Madwick Sensor Fusion for Yaw
    filter.update(gyrox, gyroy, gyroz, accelx, accely, accelz, magx, magy, magz);
    yaw = filter.getYaw();
  */

  // Yaw Calculation from Library - Calibrated
  doNMeasurements (50, magx, magy, magz);
  float yaw = atan2(magy, magx) * 180 / PI;

  /* Use for Processing Animation, for actual data, in working orientation the Roll angle is useless
   *  Also, if using processing, the roll and pitch are inverse so -roll and -pitch would be what the processing software needs to accuratly protray motion
  Serial.print(roll);
  Serial.print("/");
  */
  Serial.print(pitch-90);
  Serial.print("/");
  Serial.println(yaw);
  delay(timeInterval);
}

//Average measurements to reduce noise
void doNMeasurements(unsigned int N, float& averX, float& averY, float& averZ)
{ float x, y, z;
  averX = 0; averY = 0; averZ = 0;
  for (int i = 1; i <= N; i++)
  { while (!IMU.magnetAvailable());
    IMU.readMagnet(x, y, z);
    averX += x; averY += y;  averZ += z;
  }
  averX /= N;    averY /= N;  averZ /= N;
}
