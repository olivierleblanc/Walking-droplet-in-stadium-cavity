/*
  ADXL3xx

  Reads an Analog Devices ADXL3xx accelerometer and communicates the
  acceleration to the computer. 

  The circuit:
  - analog 0: accelerometer self test
  - analog 1: z-axis
  - analog 2: y-axis
  - analog 3: x-axis
  - analog 4: ground
  - analog 5: vcc

  http://www.arduino.cc/en/Tutorial/ADXL3xx
*/

// these constants describe the pins. They won't change:
const int groundpin = 18;             // analog input pin 4 -- ground
const int powerpin = 19;              // analog input pin 5 -- voltage
const int xInput = A3;                  // x-axis of the accelerometer
const int yInput = A2;                  // y-axis
const int zInput = A1;                  // z-axis (only on 3-axis models)
const int buttonPin = 2;
const int orangePin = 5;
const int yellowPin = 4;
const double gAcceleration = 1000.0; // [mg]

// Raw Ranges:
// initialize to mid-range and allow calibration to
// find the minimum and maximum for each axis
int xRawMin = 491;
int xRawMax = 529;
int yRawMin = 494;
int yRawMax = 530;
int zRawMin = 373;
int zRawMax = 569;
// NullAcc corresponds to the middle range in tension, then 512 (between 0 and 1023). 

// Take multiple samples to reduce noise
const double sampleSize = 1.0;


void setup() {
  analogReference(EXTERNAL);
  // initialize the serial communications:
  Serial.begin(250000);

  // Provide ground and power by using the analog inputs as normal digital pins.
  // This makes it possible to directly connect the breakout board to the
  // Arduino. If you use the normal 5V and GND pins on the Arduino,
  // you can remove these lines.
  pinMode(groundpin, OUTPUT);
  pinMode(powerpin, OUTPUT);
  pinMode(orangePin, OUTPUT);
  pinMode(yellowPin, OUTPUT);
  pinMode(buttonPin, INPUT);
  digitalWrite(orangePin, HIGH);
  digitalWrite(yellowPin, LOW);
  digitalWrite(groundpin, LOW);
  digitalWrite(powerpin, HIGH);
}

void loop() {

  /* Accelerometer reading */

  double tstart = micros();
  double xRaw = ReadAxis(xInput);
  double yRaw = ReadAxis(yInput);
  double zRaw = ReadAxis(zInput);

//    Serial.println();
//    Serial.print(xRaw);
//    Serial.print(", ");
//    Serial.print(yRaw);
//    Serial.print(", ");
//    Serial.print(zRaw);
    // Convert raw values to 'milli-Gs"
    long xScaled = map(xRaw, xRawMin, xRawMax, -gAcceleration, gAcceleration);
    long yScaled = map(yRaw, yRawMin, yRawMax, -gAcceleration, gAcceleration);
    long zScaled = map(zRaw, zRawMin, zRawMax, -gAcceleration, gAcceleration);
 
 //   Serial.print(" AccX =  ");
    Serial.println("");
    Serial.print(-xScaled);
    Serial.print(" ");
   // Serial.print("mg, ");
   // Serial.print(" AccY =  ");
    Serial.print(-yScaled);
    Serial.print(" ");
   // Serial.print("mg, ");
   // Serial.print(" AccZ =  ");
    Serial.println(-zScaled);
    Serial.print(" ");
//   // Serial.println("mg");


    /* DOD generator button pressing */
    if (digitalRead(buttonPin)==HIGH) {
      Serial.println("button has been pressed");
      digitalWrite(orangePin, LOW);
      digitalWrite(yellowPin, HIGH);
      //delay(5);
      delayMicroseconds(1000);
      digitalWrite(orangePin, HIGH);
      digitalWrite(yellowPin, LOW);
      delay(500);
    }
    double tend = micros();
//    Serial.println(tend-tstart);
// delay(200);
}

//
// Read "sampleSize" samples and report the average
//
double ReadAxis(int axisPin)
{
  long reading = 0;
  analogRead(axisPin);
  for (int i = 0; i < sampleSize; i++)
  {
    reading += analogRead(axisPin);
  }
  return reading/sampleSize;
}
