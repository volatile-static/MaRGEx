// Code to control the autotuning circuit with python.

int cPins[] = {24, 32, 30, 28, 26, 31, 29, 27, 25, 23, 33, 41, 39, 37, 35};
int nPins = 15;

void setup() {
  // Start the Serial communication
  Serial.begin(115200);
  Serial.setTimeout(100);
  // Set the pins for tuning, matching and series capacitors
  for (int c=0; c<15; c++) {
    pinMode(cPins[c], OUTPUT);
  }
}

void loop() {
  // Wait until there are any data available into the serial port
  if (Serial.available()>0) {
    delay(10);
    String state = Serial.readString();
    for (int c=0; c<nPins; c++) {
      digitalWrite(cPins[c], String(state[c]).toInt());
    }
    Serial.write("Ready!");
  }
}
