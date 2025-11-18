unsigned long previousMillis1 = 0;
unsigned long previousMillis2 = 0;
unsigned long previousMillis3 = 0;

const long interval1 = 500;  // LED 1 blink interval (500ms)
const long interval2 = 250;  // LED 2 blink interval (250ms)
const long interval3 = 1000; // LED 3 blink interval (1000ms)

void setup()
{
  pinMode(12, OUTPUT);  // LED 1
  pinMode(13, OUTPUT);  // LED 2
  pinMode(8, OUTPUT);   // LED 3
}

void loop()
{
  unsigned long currentMillis = millis();

  // LED 1 (pin 12) toggles every 200ms (5 Hz)
  if (currentMillis - previousMillis1 >= interval1) {
    previousMillis1 = currentMillis;
    digitalWrite(12, !digitalRead(12));  // Toggle LED 1
  }

  // LED 2 (pin 13) toggles every 250ms (4 Hz)
  if (currentMillis - previousMillis2 >= interval2) {
    previousMillis2 = currentMillis;
    digitalWrite(13, !digitalRead(13));  // Toggle LED 2
  }

  // LED 3 (pin 8) toggles every 300ms (3.33 Hz)
  if (currentMillis - previousMillis3 >= interval3) {
    previousMillis3 = currentMillis;
    digitalWrite(8, !digitalRead(8));    // Toggle LED 3
  }
}