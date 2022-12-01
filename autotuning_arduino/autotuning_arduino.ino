// Code to control the autotuning circuit with python.

int cPins[] = {24, 26, 28, 30, 32, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41};
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
    String state = Serial.readString();
    for (int c=0; c<nPins; c++) {
      digitalWrite(cPins[c], String(state[c]).toInt());
    }
    Serial.write("Ready!");
  }
}
