/*
 * Copyright 2016 <Admobilize>
 * All rights reserved.
 */
#include <fftw3.h>
#include <gflags/gflags.h>
#include <stdint.h>
#include <string.h>
#include <wiringPi.h>

#include <fstream>
#include <iostream>
#include <string>
#include <valarray>
#include <ctime>

#include "../cpp/driver/direction_of_arrival_music.h"
#include "../cpp/driver/everloop.h"
#include "../cpp/driver/everloop_image.h"
#include "../cpp/driver/matrixio_bus.h"
#include "../cpp/driver/microphone_array.h"
#include "../cpp/driver/microphone_core.h"

DEFINE_bool(big_menu, true, "Include 'advanced' options in the menu listing");
DEFINE_int32(sampling_frequency, 16000, "Sampling Frequency");

namespace hal = matrix_hal;

int led_offset[] = {23, 27, 32, 1, 6, 10, 14, 19};
int lut[] = {1, 2, 10, 200, 10, 2, 1};

typedef std::pair<float,float> pair;

int main(int argc, char *agrv[]) {
  google::ParseCommandLineFlags(&argc, &agrv, true);
  int duration = 10;
  if (argc > 1) {
      duration = atoi(agrv[1]); 
  }
  hal::MatrixIOBus bus;
  if (!bus.Init()) return false;

  if (!bus.IsDirectBus()) {
    std::cerr << "Kernel Modules has been loaded. Use ALSA examples "
              << std::endl;
  }

  int sampling_rate = FLAGS_sampling_frequency;

  hal::Everloop everloop;
  everloop.Setup(&bus);

  hal::EverloopImage image1d(bus.MatrixLeds());

  hal::MicrophoneArray mics(false);
  mics.Setup(&bus);
  mics.SetSamplingRate(sampling_rate);
  mics.ShowConfiguration();

  hal::MicrophoneCore mic_core(mics);
  mic_core.Setup(&bus);

  hal::DirectionOfArrivalMusic doam(mics);
  doam.Init();

  float azimutal_angle;
  float polar_angle;
  int mic;
  
  time_t start = time(NULL); //get current time
  time_t end = start; 

  while (true && (end-start) < duration) {
    mics.Read(); /* Reading 8-mics buffer from de FPGA */

    doam.Calculate();

    azimutal_angle = doam.GetAzimutalAngle() * 180 / M_PI;
    polar_angle = doam.GetPolarAngle() * 180 / M_PI;
    mic = doam.GetNearestMicrophone();

    std::cout << azimutal_angle << " " << polar_angle << " " << mic
              << std::endl;

    for (hal::LedValue &led : image1d.leds) {
      led.blue = 0;
    }

    for (int i = led_offset[mic] - 3, j = 0; i < led_offset[mic] + 3;
         ++i, ++j) {
      if (i < 0) {
        image1d.leds[image1d.leds.size() + i].blue = lut[j];
      } else {
        image1d.leds[i % image1d.leds.size()].blue = lut[j];
      }

      everloop.Write(&image1d);
    }
    end = time(NULL);
  }

  return 0;
}
