// Code to control the autotuning circuit with python.

int ttl = 22;
String tmState;
String serialState;
int nStates = 1024;
int cPins[] = {30, 32, 34, 36, 38, 13, 11, 9, 7, 5, 12, 10, 8, 6, 4};
String bestTmState;
String bestSState;
int oldVoltage;

void setup() {
  // Start the Serial communication
  Serial.begin(115200);
  Serial.setTimeout(100);
  // Set the ttl port as input
  pinMode(ttl, INPUT);
  // Set the pins for tuning, matching and series capacitors
  for (int c=0; c<15; c++) {
    pinMode(cPins[c], OUTPUT);
  }
  // Set series capacitor 0 to HIGH
  digitalWrite(cPins[0], 1);
}

void loop() {
  // Wait until there are any data available into the serial port
  if (Serial.available()>0) {
    String waste = Serial.readString();
    oldVoltage = 1023;
    // Wait until the TTL is high
    while (digitalRead(ttl)==0) {}
    if (digitalRead(ttl)==1) {
      delay(100);                           // To ensure the big relay switches on
      // Set all capacitors to False
      for (int c=0; c<15; c++) {
          digitalWrite(cPins[c], LOW);
      }
      // Test the 5 series capacitors
      for (int sCap=0; sCap<5; sCap++) {
          // Setup the series capacitor
          serialState = String(int(round(pow(2,sCap))), BIN);
          int nPins = serialState.length();
          for (int c=0; c<nPins; c++) {
              int k = nPins-c-1;
              digitalWrite(cPins[c], String(serialState[k]).toInt());
          }
          // For current series capacitor, sweep all tuning/matching states
          for (int n=0; n<nStates; n++) {
            // Setup the tuning/matching capacitors
            tmState = String(n, BIN);
            nPins = tmState.length();
            for (int c=0; c<nPins; c++) {
              int k = nPins-c-1;
              digitalWrite(cPins[5+c], String(tmState[k]).toInt());
            }
            delay(10);                         // Wait until relays are switched
            int ttlVal = 1;
            // If ttl goes down while waiting the 10 milliseconds, wait until ttl goes up again and main relay switch on
            while (digitalRead(ttl)==0) {
              ttlVal = 0;
            }
            if (ttlVal==0) {
              delay(100);
            }
            // Read voltage and check if it is lower than previous voltages
            int voltage = analogRead(A0);
            if (voltage<oldVoltage) {
              oldVoltage = voltage;
              bestTmState = tmState;
              bestSState = serialState;
            }
          }
      }

      // Set all pins to false again
      for (int c=0; c<15; c++) {
        digitalWrite(cPins[c], LOW);
      }

      // Set the best serial capacitor
      int nPins = bestSState.length();
      for (int c=0; c<nPins; c++) {
          int k = nPins-c-1;
          digitalWrite(cPins[c], String(bestSState[k]).toInt());
      }

      // Set the best tuning/matching capacitors
      nPins = bestTmState.length();
      for (int c=0; c<nPins; c++) {
          int k = nPins-c-1;
          digitalWrite(cPins[5+c], String(bestTmState[k]).toInt());
      }

      // Print result to serial port
      Serial.println(bestSState+","+bestTmState+","+oldVoltage);
    }
  }
}
