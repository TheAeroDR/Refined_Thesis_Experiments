unsigned long previousMillis1 = 0;
unsigned long previousMillis2 = 0;
unsigned long previousMillis3 = 0;

const int bottom_LED = 13;
const int middle_LED = 12;
const int top_LED = 8;

const double bottom_F = 1100.0;
const double middle_F = 1700.0;
const double top_F = 2500.0;

const double bottom_time = 1000.0 / (2.0 * bottom_F);
const double top_time = 1000.0 / (2.0 * top_F);
const double middle_time = 1000.0 / (2.0 * middle_F);

void setup()
{
  pinMode(middle_LED, OUTPUT);
  pinMode(bottom_LED, OUTPUT); 
  pinMode(top_LED, OUTPUT);
}

void loop()
{
  unsigned long currentMillis = millis();

  static bool middle_state = LOW;
  static bool bottom_state = LOW;
  static bool top_state = LOW;

  if (currentMillis - previousMillis1 >= middle_time) {
    previousMillis1 = currentMillis;
    middle_state = !middle_state;
    digitalWrite(middle_LED, middle_state);
  }

  if (currentMillis - previousMillis2 >= bottom_time) {
    previousMillis2 = currentMillis;
    bottom_state = !bottom_state;
    digitalWrite(bottom_LED, bottom_state);
  }

  if (currentMillis - previousMillis3 >= top_time) {
    previousMillis3 = currentMillis;
    top_state = !top_state;
    digitalWrite(top_LED, top_state);
  }
}