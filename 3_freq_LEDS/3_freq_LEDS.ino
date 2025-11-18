#include <avr/io.h>
#include <avr/interrupt.h>

const int bottom_LED = 6;
const int middle_LED = 11;
const int top_LED = 9;

const double bottom_F = 100.0;
const double middle_F = 200.0;
const double top_F = 300.0;

void setup()
{
  cli();

  pinMode(middle_LED, OUTPUT);
  pinMode(bottom_LED, OUTPUT); 
  pinMode(top_LED, OUTPUT);

  TCCR1A = 0;
  TCCR1B = 0;

  TCCR1A = (1<<COM1A0);
  TCCR1B = (1<<WGM12) | (1<<CS10);

  OCR1A = (uint16_t)(16000000/(2*top_F)-1);

  TCNT1 = 0;

  TCCR2A = 0;
  TCCR2B = 0;

  TCCR2A = (1<<COM2A0) | (1<<WGM21);
  TCCR2B = (1<<CS22) | (1<<CS21) | (1<<CS20);

  OCR2A = (uint8_t)(16000000/(2*1024*middle_F)-1);
  TCNT2 = 0;

  TCCR0A = 0;
  TCCR0B = 0;

  TCCR0A = (1<<COM0A0) | (1<<WGM01);
  TCCR0B = (1<<CS02) | (1<<CS10);

  OCR0A = (uint8_t)(16000000/(2*1024*bottom_F)-1);

  TCNT0 = 0;

  sei();
}

void loop()
{

}